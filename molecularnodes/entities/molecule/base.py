"""The Universe-backed ``Molecule`` entity.

Provides the ``Molecule`` class for loading and visualizing molecular structures and
molecular dynamics trajectories in Blender, using an MDAnalysis ``Universe`` as the
underlying data model.
"""

import functools
import inspect
import io
import logging
import re
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Dict, Sequence
import bpy
import databpy as db
import MDAnalysis as mda
import numpy as np
from MDAnalysis.core.groups import AtomGroup
from nodebpy.builder import MaterialBuilder
from nodebpy.nodes.geometry import NamedAttribute
from ... import download
from ...assets import data
from ...blender import coll, path_resolve, set_obj_active
from ...blender import utils as blender_utils
from ...converters import universe_from_atoms
from ...material import PresetMaterial, append_material
from ...nodes import handlers as node_handlers
from ...nodes._assets import link_asset_groups
from ...nodes._utils import STYLE_LITERALS, STYLE_NODE_MAPPING
from ...nodes.geometry import AssemblyInstance
from ...utils import (
    _UNSET,
    count_value_changes,
    is_periodic_box,
    temp_override_property,
)
from ..base import EntityType, MolecularEntity, MolecularTree
from ..utilities import (
    BoolObjectMNProperty,
    IntObjectMNProperty,
    StringObjectMNProperty,
    _validate_non_negative,
)
from .annotations import MoleculeAnnotationManager
from .calculations import CalculationManager
from .helpers import (
    FrameManager,
    _ag_to_bool,
    _isin,
)
from .selections import SelectionManager

if TYPE_CHECKING:
    from ..density import Density

logger = logging.getLogger(__name__)

#: Colour keywords understood by ``Molecule.add_style(color=...)``. Anything else must
#: name an existing attribute, be an RGBA sequence, or be a callable building the nodes.
COLOR_KEYWORDS = ("common", "default", "plddt")


# only Cycles ray-traces point clouds; every other engine draws them as a crude
# polyhedron, so a sphere style left on its "Point" default silently renders wrong
_POINT_CLOUD_ENGINES = {"CYCLES"}


def _sphere_for_engine(style_node: type, engine: str) -> str | None:
    """
    The ``sphere`` to use for a style, when the default won't render.

    Returns ``None`` when the default is fine - either the engine can draw point
    clouds, or this style doesn't default to them in the first place.
    ``StyleBallAndStick`` already defaults to ``"Instance"``, so only
    ``StyleSpheres`` is affected.
    """
    if engine in _POINT_CLOUD_ENGINES:
        return None
    parameter = inspect.signature(style_node).parameters.get("sphere")
    if parameter is None or parameter.default != "Point":
        return None
    return "Instance"


def add_style_to_tree(
    node_tree: bpy.types.GeometryNodeTree,
    style: STYLE_LITERALS = "spheres",
    material: bpy.types.Material | PresetMaterial | MaterialBuilder | str = "Default",
    color: Sequence[float] | None = None,
    selection_attribute: str | None = None,
    name: str | None = None,
) -> None:
    """
    Add a style branch directly to a node tree, without a session entity.

    The tree-level counterpart of [](`~mn.Molecule.add_style`), used by the UI
    when an object is not linked to a session entity. Everything that needs the
    entity's universe is unavailable here: selections must name an existing
    boolean attribute (MDAnalysis phrases can't be evaluated), color is limited
    to a uniform RGBA (schemes need chain information), and assembly instancing
    is not offered.

    Parameters
    ----------
    node_tree : bpy.types.GeometryNodeTree
        The tree to add the style branch to, typically the tree of an object's
        "Molecular Nodes" modifier.
    style : str, default "spheres"
        Name of a predefined style (see ``STYLE_NODE_MAPPING``).
    material : bpy.types.Material | PresetMaterial | MaterialBuilder | str, default "Default"
        Material for the styled geometry; a string is appended from the asset
        file.
    color : Sequence[float] | None, optional
        Uniform RGBA color applied via a ``Set Color`` node. ``None`` (default)
        leaves the baked ``Color`` attribute in use.
    selection_attribute : str | None, optional
        Name of an existing boolean attribute to restrict the style to.
    name : str | None, optional
        Label for the added style node.
    """
    from ...nodes import geometry as g

    if style not in STYLE_NODE_MAPPING:
        raise ValueError(
            f"Invalid style '{style}'. Supported styles are "
            f"{sorted(STYLE_NODE_MAPPING)}"
        )

    kwargs = {}
    # spheres default to point clouds, which only Cycles can draw
    geometry = _sphere_for_engine(
        STYLE_NODE_MAPPING[style], bpy.context.scene.render.engine
    )
    if geometry is not None:
        kwargs["sphere"] = geometry

    if isinstance(material, (PresetMaterial, MaterialBuilder)):
        material = material.material
    material = append_material(material) if isinstance(material, str) else material

    with MolecularTree(None, node_tree) as tree:
        selection_input = (
            NamedAttribute.boolean(selection_attribute) if selection_attribute else None
        )
        style_node = STYLE_NODE_MAPPING[style](
            selection=selection_input, material=material, **kwargs
        )
        if name:
            style_node.node.label = name
        (
            tree.atoms
            >> (g.SetColor(color=tuple(color)) if color is not None else None)
            >> style_node
            >> tree.join
        )


class Molecule(MolecularEntity):
    """Universe-backed molecular entity for Blender visualization.

    Complete interface for loading, visualizing, and manipulating molecular structures
    and molecular dynamics trajectories in Blender using an MDAnalysis ``Universe``.

    Features: structure/trajectory loading, attribute computation, selection management,
    visual styling, frame interpolation/averaging, periodic boundary handling,
    Blender animation integration.

    Attributes
    ----------
    universe : mda.Universe
        MDAnalysis Universe with topology and (optionally) trajectory frames
    frame_manager : FrameManager
        Position caching and frame updates
    selections : SelectionManager
        Dynamic atom selections
    calculations : CalculationManager
        Calculations that write attributes from the universe, such as DSSP's
        `sec_struct`
    annotations : MoleculeAnnotationManager
        Molecule annotations
    world_scale : float
        Scale factor from Angstroms to Blender units
    frame : int
        Current animation frame (synced with Blender)
    subframes : int
        Interpolation steps between frames
    offset : int
        Frame offset for playback
    average : int
        Number of frames to average (smoothing)
    correct_periodic : bool
        Apply periodic boundary corrections
    interpolate : bool
        Enable position interpolation

    Examples
    --------
    ```{python}
    #| warning: false
    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PSF, DCD
    import molecularnodes as mn
    canvas = mn.Canvas()
    u = mda.Universe(PSF, DCD)
    traj = mn.Molecule(u)
    traj.add_style("spheres", sphere="Instance", selection="resname LYS")
    canvas.look_at(traj)
    canvas.snapshot()
    ```
    """

    # Blender property descriptors with validation
    frame = IntObjectMNProperty("frame")
    subframes = IntObjectMNProperty("subframes", validate_fn=_validate_non_negative)
    offset = IntObjectMNProperty("offset")
    average = IntObjectMNProperty("average", validate_fn=_validate_non_negative)
    correct_periodic = BoolObjectMNProperty("correct_periodic")
    interpolate = BoolObjectMNProperty("interpolate")

    _mn_frame = BoolObjectMNProperty("frame_hidden")
    _mn_entity_type = StringObjectMNProperty("entity_type")
    _mn_filepath_topology = StringObjectMNProperty("filepath_topology")
    _mn_filepath_trajectory = StringObjectMNProperty("filepath_trajectory")
    _mn_n_frames = IntObjectMNProperty("n_frames", _validate_non_negative)
    _entity_type: EntityType = EntityType.MOLECULE
    # extra `mda.Universe()` arguments for subclasses whose files can't be
    # recognised from their extension, used when restoring from a saved session
    _universe_kwargs: dict = {}

    def __init__(
        self,
        universe: mda.Universe,
        name: str = "NewUniverseObject",
        world_scale: float = 0.1,
        create_object: bool = True,
    ):
        """Initialize Molecule from an MDAnalysis Universe.

        Parameters
        ----------
        universe : mda.Universe
            MDAnalysis Universe with topology and trajectory
        name : str, default="NewUniverseObject"
            Name for the Blender object
        world_scale : float, default=0.1
            Scale factor from nanometers to Blender units
        create_object : bool, default=True
            Whether to immediately create the Blender object

        Notes
        -----
        Default world_scale of 0.1 converts nanometers to Blender units.
        """
        super().__init__()
        self.universe: mda.Universe = universe
        self.selections: SelectionManager = SelectionManager(self)
        self.calculations = CalculationManager(self)
        self.world_scale = world_scale
        self._updating_in_progress = False
        self.annotations = MoleculeAnnotationManager(self)
        self.frame_manager = FrameManager(self)
        if create_object:
            self.create_object(name=name)

    @property
    def _has_periodic_box(self) -> bool:
        """Check if the current frame has a usable periodic box.

        Periodic boundary correction works for orthorhombic and triclinic
        boxes, and is skipped for frames without one.

        Returns
        -------
        bool
            True if the unit cell is a valid periodic box
        """
        return is_periodic_box(self.universe.dimensions)

    @property
    def atoms(self) -> mda.AtomGroup:
        """All atoms as MDAnalysis AtomGroup."""
        if self.universe.atoms is None:
            raise ValueError(f"Universe {self.universe} has no atoms")
        return self.universe.atoms

    @property
    def _scaled_position(self) -> np.ndarray:
        """Current atom positions in Blender world units.

        Returns
        -------
        np.ndarray
            Shape (n_atoms, 3) in Blender units
        """
        return self.atoms.positions * self.world_scale

    @property
    def uframe(self) -> int:
        """Current frame number in MDAnalysis Universe.

        Differs from `frame` property (Blender scene frame). Query actual trajectory frame.

        Returns
        -------
        int
            Current Universe frame (0-indexed)
        """
        return self.universe.trajectory.frame

    @uframe.setter
    def uframe(self, value: int) -> None:
        """Set current frame in MDAnalysis Universe.

        Only updates if changed to avoid redundant operations.

        Parameters
        ----------
        value : int
            Target frame number (0-indexed)
        """
        if self.universe.trajectory.frame != value:
            self.universe.trajectory[value]

    @functools.cached_property
    def _elements(self) -> np.ndarray:
        """
        Cached computation of element symbols for all atoms.
        This is computed once and reused by multiple attribute computations.
        """
        if hasattr(self.atoms, "elements"):
            return self.atoms.elements

        try:
            default_guesser = mda.guesser.default_guesser.DefaultGuesser(None)  # type: ignore
            guessed_elements = [
                x
                if x in data.elements.keys()
                else default_guesser.guess_atom_element(x)
                for x in self.atoms.names
            ]
            return np.array(guessed_elements)
        except Exception as e:
            logger.warning(f"Failed to compute elements, using placeholder 'X': {e}")
            return np.repeat("X", len(self))

    @property
    def _titled_elements(self) -> np.ndarray:
        # title-case the element symbols so lookups match the data tables (e.g. "FE" ->
        # "Fe"); this mirrors the biotite reader so both backends agree. Cast to a string
        # dtype first since MDAnalysis may hand back an object array.
        return np.char.title(np.asarray(self._elements, dtype=str))

    def _compute_atomic_number(self) -> np.ndarray:
        return np.array(
            [
                data.elements.get(element, {}).get("atomic_number", 0)
                for element in self._titled_elements
            ],
            dtype=int,
        )

    def _compute_vdw_radii(self) -> np.ndarray:
        return (
            np.array(
                [
                    data.elements.get(element, {}).get("vdw_radii", 100)
                    for element in self._titled_elements
                ]
            )
            * 0.01  # pm to Angstrom
            * self.world_scale  # Angstrom to world scale
        )

    def _compute_mass(self) -> np.ndarray:
        # units: daltons
        if hasattr(self.atoms, "masses"):
            return np.asarray(self.atoms.masses)
        else:
            masses = [
                data.elements.get(element, {"standard_mass": 0}).get("standard_mass")
                for element in self._elements
            ]
            return np.array(masses)

    def _compute_res_name(self) -> np.ndarray:
        return np.array(list(map(lambda x: x[0:3], self.atoms.resnames)))

    def _compute_res_name_int(self) -> np.ndarray:
        res_name = self._compute_res_name()
        return np.array(
            [
                data.residues.get(name, data.residues["UNK"]).get("res_name_num")
                for name in res_name
            ],
            dtype=int,
        )

    def _compute_b_factor(self) -> np.ndarray:
        if hasattr(self.atoms, "tempfactors"):
            return self.atoms.tempfactors
        else:
            return np.zeros(len(self))

    def _compute_occupancy(self) -> np.ndarray:
        # corresponds to Occupancies topology attr
        if hasattr(self.atoms, "occupancies"):
            return self.atoms.occupancies
        else:
            return np.zeros(len(self))

    def _compute_charge(self) -> np.ndarray:
        # corresponds to Charges topology attr. Only stored when the topology (or
        # the file, carried across by the biotite converter) provides charges; an
        # all-zero array carries no information. Skipped like any other missing
        # topology attribute; no residue / atom-name lookup is done, that will
        # return as a node instead (#1228).
        charges = self.atoms.charges  # raises mda.NoDataError when absent
        if not np.any(charges != 0):
            raise mda.NoDataError("topology provides no charges")
        return charges

    def _compute_res_id(self) -> np.ndarray:
        return self.atoms.resids

    def _compute_ures_id(self) -> np.ndarray:
        return count_value_changes(self.atoms.resids, self._compute_chain_id_int())

    def _compute_atom_id(self) -> np.ndarray:
        return self.atoms.ids

    def _compute_segindices(self) -> np.ndarray:
        segs = []
        for seg in self.atoms.segments:
            segs.append(seg.atoms[0].segid)

        else:
            try:
                self.props.segments = segs
            except db.LinkedObjectError:
                logger.warning("Failed to store segments metadata on object")

        return self.atoms.segindices

    def _compute_chain_id_int(self) -> np.ndarray:
        chain_ids, chain_id_index = np.unique(self.atoms.chainIDs, return_inverse=True)

        try:
            self.props.chain_ids = chain_ids.astype(str).tolist()
        except db.LinkedObjectError:
            logger.warning("Failed to store chain_ids metadata on object")

        return chain_id_index

    def _compute_atom_type_int(self) -> np.ndarray:
        atom_type_unique, atom_type_index = np.unique(
            self.atoms.types, return_inverse=True
        )

        try:
            self.props.atom_type_unique = atom_type_unique.tolist()
        except db.LinkedObjectError:
            logger.warning("Failed to store atom_type_unique metadata on object")

        return atom_type_index

    def _compute_atom_name_int(self) -> np.ndarray:
        if hasattr(self.atoms, "names"):
            return np.array(
                [data.atom_names.get(x, -1) for x in self.atoms.names],
                dtype=int,
            )
        else:
            return np.repeat(int(-1), len(self))

    def _compute_is_lipid(self) -> np.ndarray:
        return _isin(self.atoms.resnames, data.RESNAMES_LIPID)

    def _compute_is_solvent(self) -> np.ndarray:
        resname_is_solvent = _isin(self.atoms.resnames, data.RESNAMES_SOLVENT)
        name_is_solvent = _isin(self.atoms.names, data.NAMES_SOLVENT)
        return np.logical_or(resname_is_solvent, name_is_solvent)

    # atom names that make up the peptide and nucleic acid backbones. Matches the
    # biotite reader's definition (ReaderBase._compute_is_backbone) so that the two
    # backends agree, while "BB" preserves coarse-grained (Martini) backbone support.
    _BACKBONE_ATOM_NAMES = (
        # Peptide backbone atoms
        "N",
        "C",
        "CA",
        "H",
        "HA",
        "O",
        # Continuous nucleic backbone atoms
        "P",
        "O5'",
        "C5'",
        "C4'",
        "C3'",
        "O3'",
        # Alternative names for phosphate O's
        "O1P",
        "OP1",
        "O2P",
        "OP2",
        # Remaining ribose atoms
        "O4'",
        "C1'",
        "C2'",
        "O2'",
        # Coarse-grained backbone bead
        "BB",
    )

    def _compute_is_peptide(self) -> np.ndarray:
        # MDAnalysis' selection covers force field residue names (HIE, HSD, ...) and
        # the CCD names cover modified residues (HIC, MSE, ...)
        from biotite.structure.info import amino_acid_names

        is_ccd_amino_acid = _isin(self.atoms.resnames, list(amino_acid_names()))
        return np.logical_or(is_ccd_amino_acid, self.select("protein or (name BB SC*)"))

    def _compute_is_polymer(self) -> np.ndarray:
        # amino acid and nucleotide residues, so that ligands and ions such as calcium
        # (atom name CA) are excluded without dropping modified residues from the chain
        from biotite.structure.info import nucleotide_names

        is_ccd_nucleotide = _isin(self.atoms.resnames, list(nucleotide_names()))
        return np.logical_or.reduce(
            [
                self._compute_is_peptide(),
                is_ccd_nucleotide,
                self.select("nucleic"),
            ]
        )

    def _compute_is_alpha_carbon(self) -> np.ndarray:
        is_alpha_carbon = _isin(self.atoms.names, ("CA", "BB"))
        return np.logical_and(is_alpha_carbon, self._compute_is_polymer())

    def _compute_is_backbone(self) -> np.ndarray:
        is_backbone_atom = _isin(self.atoms.names, self._BACKBONE_ATOM_NAMES)
        return np.logical_and(is_backbone_atom, self._compute_is_polymer())

    def _compute_is_side_chain(self) -> np.ndarray:
        # side chain = polymer atoms that are not backbone, but the alpha carbon
        # (or CG backbone bead) is counted as side chain. Mirrors the biotite reader.
        backbone = self._compute_is_backbone()
        is_alpha_carbon = _isin(self.atoms.names, ("CA", "BB"))
        is_polymer = self._compute_is_polymer()
        return np.logical_and(np.logical_or(~backbone, is_alpha_carbon), is_polymer)

    def _compute_color(self) -> np.ndarray:
        from ... import color

        atomic_numbers = self._compute_atomic_number()
        # Colour by chain when chain information is available, otherwise fall back
        # to element-only colouring so that a valid `Color` attribute is *always*
        # produced. Some universes (e.g. trajectory-only formats) carry no
        # `chainIDs` topology attribute, which would otherwise raise and cause the
        # `Color` attribute to be skipped entirely.
        try:
            chain_ids = self.atoms.chainIDs
        except (mda.NoDataError, AttributeError):
            return color.colors_from_elements(atomic_numbers) / 255

        return color.color_chains(atomic_numbers, chain_ids)

    def _compute_is_hetero(self) -> np.ndarray:
        # carried across from biotite by the converter as a custom topology attribute
        return self.atoms.heteros

    def _compute_is_carb(self) -> np.ndarray:
        # carried across from biotite (filter_carbohydrates has no MDAnalysis equivalent)
        return self.atoms.is_carbs

    def _compute_entity_id(self) -> np.ndarray:
        # carried across from the structure file by the converter
        return self.atoms.entity_ids

    def _compute_sec_struct(self) -> np.ndarray:
        try:
            # carried across from the structure file by the converter
            return self.atoms.sec_structs
        except (mda.NoDataError, AttributeError):
            pass
        # otherwise compute it from the current frame so the cartoon has
        # something to show
        if len(self.universe.select_atoms("protein")) == 0:
            raise mda.NoDataError("No protein atoms to compute secondary structure")
        if self.frame_manager.n_frames != 1:
            # a trajectory gets DSSP as a calculation that isn't updating, which
            # runs on this frame and can be set to update on every frame. It is
            # set up first so this frame reuses its analysis. Setting it up
            # writes properties whose update callbacks would otherwise update
            # this molecule to the scene frame mid-import
            self._updating_in_progress = True
            try:
                self.calculations._dssp._setup("protein", updating=False, refresh=False)
            finally:
                self._updating_in_progress = False
        # raises where DSSP can't run, such as for coarse-grained beads
        try:
            return self.calculations._dssp.compute_frame()
        except Exception:
            if "sec_struct" in self.calculations:
                del self.calculations["sec_struct"]
            raise

    def _save_filepaths_on_object(self) -> None:
        """Save file paths to the Blender object for reference"""
        if isinstance(self.universe.filename, (str, Path)):
            self._mn_filepath_topology = str(path_resolve(self.universe.filename))
        if isinstance(self.universe.trajectory.filename, (str, Path)):
            self._mn_filepath_trajectory = str(
                path_resolve(self.universe.trajectory.filename)
            )

    def reset_playback(self) -> None:
        """Set the playback settings to their default values"""
        self.subframes = 0
        self.offset = 0
        self.average = 0
        self.correct_periodic = False
        self.interpolate = False

    @property
    def _blender_attributes(self) -> Dict[str, Callable | str]:
        """Registry of default attributes for Blender object.

        Defines standard molecular attributes computed at trajectory creation.
        Maps attribute names to compute functions or MDAnalysis selection strings.

        Returns
        -------
        dict
            Attribute names to compute functions or selection strings
        """
        return {
            "atomic_number": self._compute_atomic_number,
            "vdw_radii": self._compute_vdw_radii,
            "mass": self._compute_mass,
            "res_id": self._compute_res_id,
            "ures_id": self._compute_ures_id,
            "segid": self._compute_segindices,
            "res_name": self._compute_res_name_int,
            "atom_id": self._compute_atom_id,
            "b_factor": self._compute_b_factor,
            "occupancy": self._compute_occupancy,
            "charge": self._compute_charge,
            "chain_id": self._compute_chain_id_int,
            "atom_types": self._compute_atom_type_int,
            "atom_name": self._compute_atom_name_int,
            "Color": self._compute_color,
            "is_alpha_carbon": self._compute_is_alpha_carbon,
            "is_backbone": self._compute_is_backbone,
            "is_side_chain": self._compute_is_side_chain,
            "is_solvent": self._compute_is_solvent,
            "is_nucleic": "nucleic",
            "is_lipid": self._compute_is_lipid,
            "is_peptide": self._compute_is_peptide,
            "is_hetero": self._compute_is_hetero,
            "is_carb": self._compute_is_carb,
            "entity_id": self._compute_entity_id,
            "sec_struct": self._compute_sec_struct,
        }

    def _store_default_attributes(self) -> None:
        """Store default attributes"""

        for name, item in self._blender_attributes.items():
            try:
                if isinstance(item, str):
                    data = self.select(item)
                elif callable(item):
                    data = item()
                else:
                    raise ValueError("Unable to convert to attribute for storage")
                # "Color" is an (n, 4) RGBA array which would otherwise be guessed as a
                # generic 4-vector; store it as a colour attribute explicitly.
                atype = db.AttributeTypes.FLOAT_COLOR if name == "Color" else None
                self.store_named_attribute(
                    data=data,
                    name=name,
                    atype=atype,
                )
            except (mda.NoDataError, AttributeError) as e:
                logger.debug(f"Skipping attribute '{name}': {e}")
            except Exception as e:
                logger.warning(f"Failed to compute attribute '{name}': {e}")

    def _store_extra_attributes(self) -> None:
        # TODO: enable adding of arbitrary mda.Universe attirbutes not currently applied
        values = self.frame_manager.read_and_cache(self.uframe)
        self._sync_frame_attributes(values)
        values.pop("position")
        self._store_frame_attributes(values, "TRAJECTORY")

    @property
    def frame_attributes(self) -> dict[str, str]:
        """The attributes written on frame changes, and where each comes from.

        Attributes come from the trajectory (positions, and velocities, forces
        and other per-atom values a frame has), from selections and from
        `calculations`. The record is kept on the object as
        `obj.mn.frame_attributes`, and listed in the Attributes panel.

        Returns
        -------
        dict[str, str]
            Source of each attribute by name: `"TRAJECTORY"`, `"SELECTION"` or
            `"CALCULATION"`
        """
        return {item.name: item.source for item in self.object.mn.frame_attributes}

    def _sync_frame_attributes(self, trajectory_names=()) -> None:
        """Bring the object's record of attributes written on frame changes up to date.

        Trajectory attributes are only added: an attribute that a frame lacks
        keeps its last values, and a TRR file can gain velocities on later
        frames. Selections and calculations are listed as they currently are.

        Parameters
        ----------
        trajectory_names : Iterable[str]
            Names of attributes read from the trajectory on this frame
        """
        items = self.object.mn.frame_attributes
        current = [(item.name, item.source) for item in items]
        trajectory = [name for name, source in current if source == "TRAJECTORY"]
        # a single-frame structure has no trajectory to read attributes from, but
        # can still have selections and calculations
        if self.frame_manager.n_frames != 1:
            trajectory += [n for n in trajectory_names if n not in trajectory]
        wanted = [(name, "TRAJECTORY") for name in trajectory]
        wanted += [
            (item.name, "SELECTION") for item in self.selections.ui_items if item.name
        ]
        wanted += [(name, "CALCULATION") for name in self.calculations]
        # only write when something changed, as writing redraws the panel
        if wanted == current:
            return
        items.clear()
        for name, source in wanted:
            item = items.add()
            item.name = name
            item.source = source

    def _store_frame_attributes(
        self, values: dict[str, np.ndarray], source: str
    ) -> None:
        """Store attributes written on frame changes, other than positions.

        Only the attributes given are stored. One that a frame lacks is not
        given, so an existing attribute keeps the values it was last given.

        Parameters
        ----------
        values : dict[str, np.ndarray]
            Values by attribute name
        source : str
            Where the values come from: `"TRAJECTORY"`, `"SELECTION"` or
            `"CALCULATION"`
        """
        # selections are always boolean, other existing attributes keep their
        # type and new ones get the type of their values
        atype = db.AttributeTypes.BOOLEAN if source == "SELECTION" else None
        for name, value in values.items():
            try:
                self.store_named_attribute(data=value, name=name, atype=atype)
            except Exception as e:
                # values a trajectory reader has but Blender can't store are
                # expected, a calculation or selection that can't be stored is not
                expected = source == "TRAJECTORY" and isinstance(e, ValueError)
                level = logging.DEBUG if expected else logging.ERROR
                logger.log(level, f"Skipping attribute '{name}': {e}")

    def _create_object(self, name: str = "NewUniverseObject") -> None:
        """Create Blender mesh object (internal).

        Creates mesh with positions and bonds, stores attributes, sets up modifiers.

        Parameters
        ----------
        name : str
            Name for the Blender object
        """
        self.object = db.create_object(
            name=name,
            collection=coll.mn(),
            vertices=self._scaled_position,
            edges=self.atoms.bonds.indices if hasattr(self.atoms, "bonds") else None,
        )

        # carry bond order/type onto the edge domain (set by the biotite converter,
        # aligned to the universe's bond ordering used for the edges above)
        bond_types = getattr(self.universe, "_mn_bond_types", None)
        if bond_types is not None:
            self.store_named_attribute(
                data=bond_types,
                name="bond_type",
                domain=db.AttributeDomains.EDGE,
                atype=db.AttributeTypes.INT,
            )

        self._mn_entity_type = self._entity_type.value
        try:
            self._mn_n_frames = self.universe.trajectory.n_frames
        except RuntimeError:
            pass

    def create_object(self, name: str = "NewUniverseObject") -> bpy.types.Object:
        """Create and initialize Blender object for trajectory.

        Creates mesh, computes attributes, sets up modifiers, registers with MolecularNodes.

        Parameters
        ----------
        name : str, default="NewUniverseObject"
            Name for the Blender object

        Returns
        -------
        bpy.types.Object
            Created Blender object
        """
        self._create_object(name=name)
        self._store_default_attributes()
        self._store_extra_attributes()
        self._setup_modifiers()
        self._save_filepaths_on_object()
        set_obj_active(self.object)
        return self.object

    @classmethod
    def load(
        cls,
        topology: Path | str,
        coordinates: Path | str | None = None,
        style: STYLE_LITERALS | None = None,
        selection: str | None = None,
        create_object: bool = True,
        name: str | None = None,
        data_block: str | None = None,
        **kwargs,
    ) -> "Molecule":
        """Load a single structure file, or an MD topology + trajectory.

        With only ``topology`` given, it is treated as a single structure file
        (``.pdb``/``.cif``/``.bcif``/``.sdf``/``.mol``/``.cube``) and routed through the biotite
        reader and converter (see :meth:`from_file`). When ``coordinates`` is also
        given, the two are read as an MD topology and trajectory into an
        MDAnalysis ``Universe``.

        Parameters
        ----------
        topology : Path | str
            Structure file, or MD topology file when ``coordinates`` is given.
        coordinates : Path | str | None, optional
            MD trajectory/coordinates file. If omitted, ``topology`` is loaded as a
            single structure file.
        name : str | None, optional
            Name for the Blender object. Defaults to the topology file name
            (extension included).
        style : str | None, optional
            If given, the visual style to apply to the loaded entity. If None (the
            default) no style is added, leaving the node tree empty for manual setup.
        selection : str | None, optional
            Atom selection to restrict the style to, passed to :meth:`add_style`.
        create_object : bool, optional
            Whether to create the Blender object immediately (MD route only).
        data_block : str | None, optional
            For a ``.cif``/``.bcif`` with several structures in separate data blocks,
            the name of the block to load (see :meth:`from_file`).
        kwargs : dict, optional
            Additional keyword arguments to pass to the `MDAnalysis.Universe()` constructor.

        Returns
        -------
        Molecule
            The loaded entity.
        """
        if coordinates is None:
            entity = cls.from_file(topology, name=name, data_block=data_block)
        else:
            entity = cls(
                universe=mda.Universe(topology, coordinates, **kwargs),
                name=name or Path(topology).name,
                create_object=create_object,
            )

        if style is not None and create_object:
            entity.add_style(style=style, selection=selection)

        return entity

    @classmethod
    def from_file(
        cls,
        file_path: str | Path | io.BytesIO,
        name: str | None = None,
        data_block: str | None = None,
    ) -> "Molecule":
        """Load a single structure file into a Universe-backed entity.

        The file (``.pdb``/``.cif``/``.bcif``/``.sdf``/``.mol``/``.cube``) is parsed by the biotite
        readers and converted into an MDAnalysis ``Universe`` via
        :func:`~molecularnodes.converters.universe_from_atoms`. Multi-model files become
        multi-frame universes. Biological assembly and entity/chain metadata parsed from
        the file are stored on the Blender object.

        Small-molecule crystallographic ``.cif`` files (e.g. from the ICSD or COD)
        are instead read by :mod:`~molecularnodes.entities.molecule.corecif`, which
        expands the asymmetric unit by the file's symmetry operators to fill one
        unit cell. ``.xyz`` files are read natively by MDAnalysis, with elements
        taken from the atom names.

        Parameters
        ----------
        file_path : str | Path | io.BytesIO
            Path to the structure file (or an in-memory ``bcif`` buffer).
        name : str | None, optional
            Name for the Blender object. Defaults to the file name (extension included).
        data_block : str | None, optional
            For a ``.cif``/``.bcif`` with several structures in separate data blocks,
            the name of the block to load. Blocks without atoms (such as restraint
            definitions) are skipped, so this is only needed when more than one
            block holds a structure; otherwise a
            :class:`~molecularnodes.entities.molecule.pdbx.MultipleStructuresError`
            listing the blocks is raised.

        Returns
        -------
        Molecule
            The Universe-backed entity representing the structure.
        """
        from ..molecule.reader import read_structure
        from . import corecif

        # small-molecule crystallographic CIFs (ICSD, COD) are a different
        # dialect that biotite cannot parse; route them through the dedicated
        # MDAnalysis parser, which expands the asymmetric unit to a unit cell
        if (
            not isinstance(file_path, io.BytesIO)
            and Path(file_path).suffix == ".cif"
            and corecif.is_core_cif(file_path)
        ):
            universe = mda.Universe(
                str(file_path),
                topology_format=corecif.CoreCIFParser,
                format=corecif.CoreCIFReader,
            )
            return cls(universe, name=name or Path(file_path).name)

        # .xyz files are read natively by MDAnalysis; the atom names are the
        # element symbols, which the element-derived attributes need
        if not isinstance(file_path, io.BytesIO) and Path(file_path).suffix == ".xyz":
            universe = mda.Universe(str(file_path))
            if not hasattr(universe.atoms, "elements"):
                universe.add_TopologyAttr(
                    "elements",
                    [
                        re.sub(r"[^A-Za-z]", "", atom_name)[:2].capitalize()
                        for atom_name in universe.atoms.names
                    ],
                )
            return cls(universe, name=name or Path(file_path).name)

        reader = read_structure(file_path, data_block=data_block)
        universe = universe_from_atoms(reader.array, frame_data=reader.model_data())
        if name is None:
            name = Path(file_path).name if not isinstance(file_path, io.BytesIO) else ""
        entity = cls(universe, name=name)
        entity._store_structure_metadata(reader, file_path)
        return entity

    @classmethod
    def fetch(
        cls,
        code: str,
        format: str = ".bcif",
        cache: Path | str | None = download.CACHE_DIR,
        database: str = "rcsb",
    ) -> "Molecule":
        """Fetch a structure from an online database into a Universe-backed entity.

        Parameters
        ----------
        code : str
            The database accession code (e.g. a PDB id), or for ``"pubchem"`` a
            compound CID (e.g. ``"2244"``) or name (e.g. ``"aspirin"``).
        format : str, optional
            File format to download, by default ``".bcif"``. Ignored for
            ``"pubchem"``, which always downloads the 3D conformer as SDF.
        cache : Path | str | None, optional
            Directory to cache downloads in.
        database : str, optional
            The database to fetch from: ``"rcsb"`` (default), ``"alphafold"`` or
            ``"pubchem"``.

        Returns
        -------
        Molecule
            The Universe-backed entity representing the fetched structure.
        """
        file_path = download.StructureDownloader(cache=cache).download(
            code=code, format=format, database=database
        )
        entity = cls.from_file(file_path, name=code)
        # record the source so the entity can be re-fetched into a fresh session
        entity.props.code = code
        entity.props.database = database
        return entity

    def _store_structure_metadata(
        self, reader, file_path: str | Path | io.BytesIO
    ) -> None:
        """Store file-parsed assembly/entity/chain metadata on the Blender object."""
        self.props.entity_type = EntityType.MOLECULE
        self.props.entity_ids = reader.entity_ids()
        self.props.chain_ids = reader.chain_ids()
        self.props.biological_assemblies = reader.assemblies(as_json_string=True)
        # record the source so the entity can be reloaded into a fresh session
        if not isinstance(file_path, io.BytesIO):
            self.props.filepath = str(file_path)
        self.props.data_block = getattr(reader, "_data_block", None) or ""

    def assemblies(self, as_array: bool = False):
        """Get the biological assemblies parsed from the source file.

        Parameters
        ----------
        as_array : bool, optional
            Return the assemblies as a structured array of per-chain 4x4 transforms
            rather than a dict.

        Returns
        -------
        dict | np.ndarray | None
            The biological assemblies as transformation matrices, or ``None`` when
            the structure has no assembly data.
        """
        from ... import utils

        assemblies_info = self.props.biological_assemblies
        if not assemblies_info:
            return None
        if as_array:
            try:
                return utils.array_transforms_from_dict(assemblies_info)
            except (ValueError, TypeError) as e:
                logger.warning(f"Failed to convert biological assemblies to array: {e}")
                return None
        return assemblies_info

    def align_to(
        self,
        fixed: "Molecule",
        method: str = "sequence",
        selection: str | None = None,
        fixed_selection: str | None = None,
    ):
        """Superimpose this molecule onto ``fixed`` by setting its object transform.

        Residues are paired by sequence (``method="sequence"``, the default) or by
        structure (``method="structure"``, proteins only), so the two molecules do
        not need identical sequences. The coordinates are not changed. See
        :func:`~molecularnodes.entities.molecule.align.align` for details.

        Parameters
        ----------
        fixed : Molecule
            The reference molecule, which stays where it is.
        method : {"sequence", "structure"}, optional
            How corresponding residues are found.
        selection : str | None, optional
            MDAnalysis selection of this molecule's atoms to fit on.
        fixed_selection : str | None, optional
            MDAnalysis selection of the fixed molecule's atoms to fit on.

        Returns
        -------
        AlignResult
            The fitted transform, RMSD and number of anchor atoms.
        """
        from .align import align

        return align(
            self,
            fixed,
            method=method,
            selection=selection,
            fixed_selection=fixed_selection,
        )

    def create_data_object(self) -> bpy.types.Object:
        """Create the data object holding the biological assembly transforms."""
        from ...blender import mesh

        data_obj_name = f".data_{self.name}_assemblies"
        data_obj = bpy.data.objects.get(data_obj_name)
        if not data_obj:
            transforms = self.assemblies(as_array=True)
            data_obj = mesh.create_data_object(array=transforms, name=data_obj_name)

        return data_obj

    def _update_attributes(self, frame: int) -> None:
        """Write the attributes of every source for the scene frame.

        The trajectory goes first: it moves the universe to the frame and sets
        positions, which selections and calculations are then computed from.

        Parameters
        ----------
        frame : int
            Scene frame number
        """
        trajectory = self._update_trajectory_positions(frame)
        self._store_frame_attributes(trajectory, "TRAJECTORY")
        self._store_frame_attributes(self.selections._frame_values(), "SELECTION")
        # DSSP and hydrogen bonds are rebuilt from their saved settings after a
        # file is reopened
        self.calculations._dssp._restore()
        self.calculations._hbonds._restore()
        self._store_frame_attributes(self.calculations._frame_values(), "CALCULATION")
        self._sync_frame_attributes(["position", *trajectory])

    def set_frame(self, frame: int) -> None:
        """Update trajectory state for scene frame.

        Main entry point called by Blender's animation system. Updates positions,
        selections, and calculations with recursion prevention.

        Parameters
        ----------
        frame : int
            Scene frame number (mapping applied to get Universe frame)

        Notes
        -----
        Typically called automatically by frame change handlers, not user code.
        """
        # single-frame structures have no trajectory to update; skip so that
        # manually modified positions aren't reset on scene frame changes.
        # n_frames is None for streaming trajectories, which must still update.
        n_frames = self.frame_manager.n_frames
        if n_frames is not None and n_frames <= 1:
            return
        if self._updating_in_progress:
            logger.debug("Update already in progress, skipping nested update")
            return

        try:
            self._updating_in_progress = True
            self._update_attributes(frame)
            # update annotation object
            self.annotations._update_annotation_object()
            # update periodic box nodes
            self._update_box()
        finally:
            self._updating_in_progress = False

    def frames_to_collection(
        self,
        start: int = 0,
        stop: int | None = None,
        step: int = 1,
        name: str | None = None,
    ) -> bpy.types.Collection:
        """Bake a range of trajectory frames into a collection of frame objects.

        Each object in the returned collection holds the atom positions (in Blender world
        units) for a single frame, in the same vertex order as this molecule's mesh,
        along with the attributes written on frame changes (see `frame_attributes`):
        those read from the trajectory, selections and calculations. The collection
        can then be read inside geometry nodes (e.g. with the *Animate Frames* node)
        to drive positions from baked data, rather than updating the ``Universe``
        every scene frame.

        As in playback, a trajectory attribute that a frame lacks keeps the values
        of the last frame that had it, and selections and calculations that don't
        update on frame changes keep the values they have now.

        Re-baking replaces any objects previously baked into the same collection.

        Parameters
        ----------
        start : int, default=0
            First trajectory frame to bake (inclusive).
        stop : int | None, optional
            One past the last frame to bake. Defaults to the number of frames in the
            trajectory (i.e. bake through the final frame).
        step : int, default=1
            Stride between baked frames.
        name : str | None, optional
            Name for the frames collection. Defaults to this molecule's name.

        Returns
        -------
        bpy.types.Collection
            The collection of baked frame objects.

        Raises
        ------
        ValueError
            If ``step`` is not a positive integer.
        """
        if step < 1:
            raise ValueError("step must be a positive integer")

        n_frames = self.universe.trajectory.n_frames
        if stop is None:
            stop = n_frames
        start = max(0, min(start, n_frames))
        stop = max(0, min(stop, n_frames))

        frames = coll.frames(name or self.name)
        # clear any previously-baked frames so re-baking replaces rather than appends
        for obj in list(frames.objects):
            bpy.data.objects.remove(obj, do_unlink=True)

        attributes = self.object.data.attributes
        static = {
            name: db.named_attribute(self.object, name)
            for name in self.calculations
            if not self.calculations.updating(name) and name in attributes
        }
        values: dict[str, np.ndarray] = {}
        sources: dict[str, str] = {}

        # bake each requested frame, restoring the current frame afterwards
        original_uframe = self.uframe
        try:
            for i, frame in enumerate(range(start, stop, step)):
                self.uframe = frame
                obj = db.create_object(
                    vertices=self._scaled_position,
                    name=f"{self.name}_frame_{i:04d}",
                    collection=frames,
                )
                # trajectory values a frame lacks carry over from earlier frames
                trajectory = self.frame_manager._read_frame(frame)
                trajectory.pop("position")
                selections = {
                    name: _ag_to_bool(ag)
                    for name, ag in self.selections.atomgroups.items()
                }
                calculations = {**self.calculations._frame_values(), **static}
                for source, frame_values in (
                    ("TRAJECTORY", trajectory),
                    ("SELECTION", selections),
                    ("CALCULATION", calculations),
                ):
                    values.update(frame_values)
                    sources.update(dict.fromkeys(frame_values, source))
                self._store_baked_attributes(obj, values, sources)
        finally:
            self.uframe = original_uframe

        return frames

    def _store_baked_attributes(
        self,
        obj: bpy.types.Object,
        values: dict[str, np.ndarray],
        sources: dict[str, str],
    ) -> None:
        """Store the per-frame attributes of a baked frame on its object.

        Each attribute gets the type it has on this molecule's mesh.

        Parameters
        ----------
        obj : bpy.types.Object
            The baked frame's object
        values : dict[str, np.ndarray]
            Values by attribute name
        sources : dict[str, str]
            Source of each attribute, by name: `"TRAJECTORY"`, `"SELECTION"` or
            `"CALCULATION"`
        """
        attributes = self.object.data.attributes
        for name, value in values.items():
            attribute = attributes.get(name)
            if attribute is not None:
                atype = attribute.data_type
            elif sources[name] == "SELECTION":
                atype = db.AttributeTypes.BOOLEAN
            else:
                atype = None
            try:
                db.store_named_attribute(obj, data=value, name=name, atype=atype)
            except Exception as e:
                # as in playback, trajectory values Blender can't store are expected
                expected = sources[name] == "TRAJECTORY" and isinstance(e, ValueError)
                level = logging.DEBUG if expected else logging.ERROR
                logger.log(level, f"Skipping baked attribute '{name}': {e}")

    def _update_trajectory_positions(self, frame: int) -> dict[str, np.ndarray]:
        """Update trajectory positions for the given frame.

        This method can be overridden by subclasses to implement custom
        position update logic (e.g., for streaming trajectories).

        Parameters
        ----------
        frame : int
            Scene frame number

        Returns
        -------
        dict[str, np.ndarray]
            Values of the frame's other trajectory attributes, to be stored
        """
        values = self.frame_manager.get_values_at_frame(frame)
        self.position = values.pop("position")
        return values

    def _update_box(self) -> None:
        """Update any Periodic Box nodes in the geometry node tree."""
        dimensions: tuple[float] | None = self.universe.trajectory.ts.dimensions
        if dimensions is None:
            return
        names = ["a", "b", "c", "alpha", "beta", "gamma"]
        nodes_to_update = ["Periodic Box", "Periodic Array"]
        for node in self.modifier_node_tree.nodes:
            if (
                not isinstance(node, bpy.types.GeometryNodeGroup)
                or node.node_tree is None
                or node.node_tree.name not in nodes_to_update
                or not node.inputs["Update"].default_value
            ):
                continue

            for name, value in zip(names, dimensions):
                node.inputs[name].default_value = value

    def __repr__(self) -> str:
        return f"<Molecule, `universe`: {self.universe}, `object`: {self.object}"

    def select(
        self, selection: str | AtomGroup | np.ndarray | None = None
    ) -> np.typing.NDArray[np.bool_]:
        """A per-atom boolean mask from a selection.

        Use it to index any attribute, combining masks with ``&``, ``|`` and ``~``.
        Indexing an attribute writes the changed values back to Blender::

            mol["b_factor"][mol.select("resid 19:26")] = 1.0
            mask = mol.select("protein") & ~mol.select("name H*")

        A new attribute has to be stored for every atom, such as with `np.where`::

            mol["highlight"] = np.where(mol.select("resid 19:26"), 1.0, 0.0)
            mol["is_site"] = mol.select("resid 19:26")

        For the ``Color`` attribute, `set_color` also takes colour names.

        Parameters
        ----------
        selection : str | AtomGroup | np.ndarray | None
            A string names an existing boolean attribute, or is otherwise an
            MDAnalysis selection phrase. An array is a boolean mask or atom indices.
            None selects every atom.

        Returns
        -------
        np.ndarray
            Boolean mask, True for the selected atoms

        Raises
        ------
        TypeError
            If the string names an attribute that isn't boolean
        """
        if isinstance(selection, AtomGroup):
            return _ag_to_bool(selection)
        if isinstance(selection, str):
            # an existing attribute is used directly, as in `add_style`
            if selection in self.list_attributes(drop_hidden=False):
                values = self.named_attribute(selection)
                if values.dtype != bool:
                    raise TypeError(
                        f"Attribute '{selection}' is {values.dtype}, not a boolean "
                        "attribute, so it can't be used as a selection"
                    )
                return values
            return _ag_to_bool(self.universe.select_atoms(selection))
        mask = np.zeros(self.universe.atoms.n_atoms, dtype=bool)
        mask[slice(None) if selection is None else selection] = True
        return mask

    def set_color(
        self,
        color: str | Sequence[float],
        selection: str | AtomGroup | np.ndarray | None = None,
    ) -> None:
        """Set the ``Color`` attribute of the selected atoms.

        Styles added without a ``color`` show this attribute::

            mol.set_color((0.1, 0.7, 0.7), selection="chainID A and resid 19:26")
            mol.set_color("#00FFCC", "protein")
            mol.set_color("pymol:slate", "chainID B")

        Parameters
        ----------
        color : str | Sequence[float]
            A hex string such as ``"#00FFCC"`` or a colour name such as
            ``"salmon"`` or ``"pymol:slate"``, which are sRGB and converted to
            linear (see `color.parse_color`), or linear RGB or RGBA from 0 to 1.
            Colours without an alpha keep the atoms' current alpha.
        selection : str | AtomGroup | np.ndarray | None
            Atoms to color, as for `select`. None colors every atom.
        """
        from ... import color as mn_color

        if isinstance(color, str):
            color = mn_color.parse_color(color)
        if len(color) not in (3, 4):
            raise ValueError(f"Color must be RGB or RGBA, not {color}")
        self["Color"][self.select(selection), : len(color)] = color

    def _resolve_style_selection(self, selection: str | AtomGroup | None) -> str | None:
        """Resolve an ``add_style`` selection to a boolean-attribute name (or None).

        A string is treated as an existing attribute name first; if no such attribute
        exists it is interpreted as an MDAnalysis selection phrase and stored as a new
        managed selection. A string that is neither raises a ``UserWarning`` and results
        in no selection. An ``AtomGroup`` always becomes a new managed selection.
        """
        if selection is None:
            return None
        if callable(selection):
            # evaluated later, inside the tree context, by `add_style`
            return None
        if isinstance(selection, AtomGroup):
            return self.selections.from_atomgroup(selection).name
        if isinstance(selection, str):
            # an existing boolean attribute is used directly, without creating a
            # managed selection
            if selection in self.list_attributes(drop_hidden=False):
                return selection
            # otherwise interpret it as an MDAnalysis selection phrase, validating it
            # first so an invalid phrase warns rather than storing a broken selection
            try:
                self.universe.select_atoms(selection)
            except Exception:
                warnings.warn(
                    f"Selection '{selection}' is neither an existing named attribute "
                    "nor a valid MDAnalysis selection. The style will be added but "
                    "nothing will be displayed unless that attribute is created.",
                    category=UserWarning,
                )
                return None
            return self.selections.from_string(selection).name
        return None

    @staticmethod
    def _style_color_groups(color) -> tuple[type, ...]:
        "Asset groups that `_style_color_input` builds for a color keyword."
        from ...nodes import geometry as g

        if not isinstance(color, str):
            return ()
        if color.lower() in ("common", "default"):
            return (g.ColorElement, g.RandomColor, g.ChainID)
        if color.lower() == "plddt":
            return (g.ColorPLDDT,)
        return ()

    def _style_color_input(self, color: str | Sequence[float] | Callable):
        """Resolve `add_style`'s color argument into the `Set Color` node input.

        Returns ``None`` if no ``Set Color`` node should be added, which happens when
        the argument names neither a known keyword, an existing attribute nor a colour.
        """
        from ... import color as mn_color
        from ...nodes import geometry as g

        # a callable is evaluated here, inside the tree context, so that it can build
        # nodes from `molecularnodes.nodes.geometry` (which require that context)
        if callable(color):
            return color()
        if not isinstance(color, str):
            return tuple(color)
        if color.lower() in ("common", "default"):
            # standard element colors, with carbons colored randomly per chain
            return g.ColorElement(c=g.RandomColor(g.ChainID(), 3))
        if color.lower() == "plddt":
            return g.ColorPLDDT()
        # otherwise it must name an existing *color* attribute on the geometry. Reading
        # a non-color attribute (a float, say) as a color silently renders black, so
        # check the data type too and not just that the name exists
        color_attributes = [
            attribute.name
            for attribute in self.object.data.attributes
            if attribute.data_type in ("FLOAT_COLOR", "BYTE_COLOR")
        ]
        if color in color_attributes:
            return NamedAttribute.color(color)
        # otherwise a hex colour or colour name, such as "#00FFCC" or "pymol:slate"
        try:
            rgba = mn_color.parse_color(color)
        except ValueError:
            pass
        else:
            alpha = rgba[3] if len(rgba) == 4 else 1.0
            return (*rgba[:3], alpha)
        warnings.warn(
            f"Color '{color}' is neither a known color keyword {COLOR_KEYWORDS}, a "
            f"color attribute on this molecule (available: {color_attributes}) nor a "
            "hex color or color name (see `mn.color.parse_color`). No color will be "
            "applied. For coloring beyond the keywords, pass a callable "
            "that builds the color nodes:\n"
            "    from molecularnodes.nodes import geometry as mg\n"
            "    mol.add_style('cartoon', color=lambda: mg.ColorRainbow())",
            category=UserWarning,
        )
        return None

    def add_style(
        self,
        style: STYLE_LITERALS | Callable = "spheres",
        selection: str | AtomGroup | Callable | None = None,
        material: bpy.types.Material
        | PresetMaterial
        | MaterialBuilder
        | str
        | None = _UNSET,
        color: "str | Sequence[float] | Callable | Density | None" = None,
        assembly: bool = False,
        name: str | None = None,
        **kwargs,
    ) -> "Molecule":
        """
        Add a visual style to the trajectory.

        Provides a simple interface for adding visual styles to the molecule. For more complex
        styling, use the manual node tree creation via the `with mol.tree:` context manager.

        Parameters
        ----------
        style : str | Callable, optional
            The style to apply. Either a string naming a predefined style
            ("spheres", "cartoon", "ribbon", "surface", "sticks", "ball_and_stick"),
            or a zero-argument callable returning a style node, evaluated inside the
            node tree context::

                mol.add_style(lambda: mg.StyleCartoon(quality=5, loop_radius=0.6))

            A callable defines the style node itself, so ``selection``, ``material``
            and extra keyword arguments cannot be combined with it - set them inside
            the callable. Default is "spheres".

        selection : str | AtomGroup | Callable | None, optional
            Apply the style only to atoms matching this selection. Can be:
            - A string naming an existing boolean attribute on the molecule (used directly)
            - A string MDAnalysis selection phrase (evaluated and stored as a new
              managed selection attribute)
            - A AtomGroup object defining a selection criteria
            - A zero-argument callable returning a boolean socket, evaluated inside
              the node tree context, e.g.
              ``lambda: mg.IsPeptide() & mg.IsSideChain()``
            - None to apply to all atoms (default)

            A string is treated as an existing attribute name first; only if no such
            attribute exists is it interpreted as an MDAnalysis selection phrase.

        material : bpy.types.Material | PresetMaterial | MaterialBuilder | str | None, optional
            The material to apply to the styled atoms. Can be one of the preset
            materials from `mn.material` (e.g. ``mn.material.AmbientOcclusion()``),
            a Blender Material object, a nodebpy MaterialBuilder, a string with a
            material name to append from the asset file, or None for no material.
            Default is "Default".

        color : str | Sequence[float] | Callable | None, optional
            Coloring to apply upstream of the style via a ``Set Color`` node. Can be:
            - ``"common"`` / ``"default"`` for standard element colors with carbons
              colored randomly per chain
            - ``"plddt"`` to color by pLDDT (B-factor) confidence
            - the name of an existing *color* attribute on the geometry
            - a hex color or color name for a single uniform color, such as
              ``"#00FFCC"``, ``"salmon"`` (CSS) or ``"pymol:slate"`` (see
              `color.parse_color`)
            - an RGBA sequence of floats for a single uniform color
            - a zero-argument callable returning a color socket, evaluated inside the
              node tree context, e.g. ``lambda: mg.ColorRainbow()``. This is the way
              to reach the full set of ``Color*`` nodes without writing out a whole
              node tree.
            - a density entity (e.g. an APBS potential loaded with
              ``mn.entities.density.Grids.load``) to color the styled geometry by
              sampling its grid, through a ``Color Sample Grid`` node. Unlike the
              other options this is applied *after* the style, at every point of
              its output (each surface vertex for ``"surface"``), so instanced
              geometry such as instanced spheres is not colored. The range is the
              node's default of -5 to 5 (red through light grey to blue); change
              ``Min`` and ``Max`` on the node for other grids. A density keeps
              its grid in the bundle of its styled output, which the node samples,
              so the density can keep its own style.
            - None (default) to add no color node, leaving the baked ``Color``
              attribute in use.

            A string that is neither a keyword nor an existing color attribute raises
            a ``UserWarning`` and no color is applied.

        assembly : bool, optional
            Instance the style over the biological assembly transforms parsed from
            the source file, via an ``Assembly Instance`` node. Default is False.

        name : str | None, optional
            Optional label for the added style node, shown in the node editor and
            style lists. Default is None.

        **kwargs : optional
            Additional keyword arguments passed to the style node, matching that
            node's inputs (e.g. ``quality``, ``scale``, ``sphere`` for
            ``spheres``). Unknown names raise a ``TypeError``. Cannot be combined
            with a callable ``style``.

        Returns
        -------
        Molecule
            Returns self for method chaining.

        Raises
        ------
        ValueError
            If an unsupported style string is passed
        TypeError
            If a callable ``style`` is combined with ``selection``, ``material`` or
            extra keyword arguments, or if a keyword argument does not match an input
            on the style node

        Notes
        -----
        If a selection is provided, it will be evaluated and stored as a new
        named attribute on the trajectory with an automatically generated name (sel_N).
        """

        from ...nodes import geometry as g
        from ..density import Density
        from . import OXDNA

        color_by_density = isinstance(color, Density)
        style_is_callable = callable(style)
        if not style_is_callable and style not in STYLE_NODE_MAPPING:
            raise ValueError(
                f"Invalid style '{style}'. Supported styles are "
                f"{sorted(STYLE_NODE_MAPPING)}"
            )
        if style_is_callable:
            # a callable builds the style node itself, so arguments belonging to that
            # node would be silently discarded - reject them instead of ignoring them
            owned = [
                argname
                for argname, was_given in (
                    ("selection", selection is not None),
                    ("material", material is not _UNSET),
                    *((key, True) for key in kwargs),
                )
                if was_given
            ]
            if owned:
                raise TypeError(
                    "When `style` is a callable it defines the style node itself, so "
                    f"{owned} cannot also be passed to `add_style()`. Set them inside "
                    "the callable instead:\n"
                    "    mol.add_style(lambda: mg.StyleCartoon(quality=5))"
                )

        material = "Default" if material is _UNSET else material
        # a callable selection is evaluated inside the tree context, below
        selection_is_callable = callable(selection)
        attribute_name = self._resolve_style_selection(selection)

        # oxDNA nucleotides get their own ribbon; override locally so the
        # module-level mapping used by every other entity is left untouched
        style_mapping = STYLE_NODE_MAPPING
        if isinstance(self, OXDNA):
            style_mapping = {**STYLE_NODE_MAPPING, "ribbon": g.OxDNAStyleClassic}

        if not style_is_callable and "sphere" not in kwargs:
            # spheres default to point clouds, which only Cycles can draw
            geometry = _sphere_for_engine(
                style_mapping[style], bpy.context.scene.render.engine
            )
            if geometry is not None:
                kwargs["sphere"] = geometry

        if isinstance(material, (PresetMaterial, MaterialBuilder)):
            material = material.material

        material = append_material(material) if isinstance(material, str) else material

        # link the asset groups this style is about to use with one library read
        link_asset_groups(
            [
                None if style_is_callable else style_mapping[style],
                g.SetColor if color is not None else None,
                *self._style_color_groups(color),
                g.ColorSampleGrid if color_by_density else None,
                g.AssemblyInstance if assembly else None,
            ]
        )

        with self.tree as tree:
            if style_is_callable:
                style_node = style()
            else:
                if selection_is_callable:
                    selection_input = selection()
                elif attribute_name:
                    selection_input = NamedAttribute.boolean(attribute_name)
                else:
                    selection_input = None
                style_node = style_mapping[style](
                    selection=selection_input,
                    material=material,
                    **kwargs,
                )
            if name:
                style_node.node.label = name
            color_input = (
                self._style_color_input(color)
                if color is not None and not color_by_density
                else None
            )
            (
                tree.atoms
                >> (g.SetColor(color=color_input) if color_input is not None else None)
                >> style_node
                >> (
                    g.SetColor(color=g.ColorSampleGrid(object=color.object).o.color)
                    if color_by_density
                    else None
                )
                >> (
                    g.AssemblyInstance(data_object=self.create_data_object())
                    if assembly
                    else None
                )
                >> tree.join
            )

        return self

    def __getstate__(self):
        """Custom serialization to handle MDAnalysis Universe objects."""
        state = super().__getstate__()

        # Store universe file paths for restoration
        if hasattr(self, "universe") and self.universe is not None:
            try:
                topology_filename = getattr(self.universe, "filename", None)
                trajectory_filename = getattr(
                    self.universe.trajectory, "filename", None
                )

                if isinstance(topology_filename, (str, Path)):
                    # MD universe backed by on-disk topology (+ trajectory) files
                    state["_universe_topology"] = str(topology_filename)
                    # multi-file universes (ChainReader) list every coordinate
                    # file in `filenames`, while `filename` only holds the first
                    trajectory_filenames = getattr(
                        self.universe.trajectory, "filenames", None
                    )
                    if trajectory_filenames is not None:
                        state["_universe_trajectory"] = [
                            str(f) for f in trajectory_filenames
                        ]
                    elif isinstance(trajectory_filename, (str, Path)):
                        state["_universe_trajectory"] = str(trajectory_filename)
                else:
                    # universe was converted from a structure file (in-memory, so the
                    # filename is a BiotiteWrapper or None); record the source file/code
                    # so it can be rebuilt via the biotite reader + converter on restore.
                    if self.props.filepath:
                        state["_structure_filepath"] = self.props.filepath
                        state["_structure_data_block"] = self.props.data_block
                    if self.props.code:
                        state["_structure_code"] = self.props.code
                        state["_structure_database"] = self.props.database

                state["_universe_frame"] = self.universe.trajectory.frame

            except AttributeError as e:
                logger.warning(
                    f"Could not extract file paths from universe during serialization: {e}"
                )

            del state["universe"]

        # Remove objects with circular references or PyCapsules
        if "frame_manager" in state:
            del state["frame_manager"]
        if "selections" in state:
            del state["selections"]
        if "annotations" in state:
            del state["annotations"]
        if "calculations" in state:
            # Preserve picklable calculations
            preserved_calculations = {}
            for name, calc_func in state["calculations"].items():
                # built-in calculations are rebuilt from their saved settings
                if name in ("sec_struct", "hbonds"):
                    continue
                try:
                    import pickle

                    pickle.dumps(calc_func)
                    preserved_calculations[name] = calc_func
                except (TypeError, AttributeError):
                    logger.debug(f"Skipping unpicklable calculation function: {name}")

            if preserved_calculations:
                state["_preserved_calculations"] = preserved_calculations
            del state["calculations"]

        return state

    def __setstate__(self, state):
        """Custom deserialization to recreate MDAnalysis Universe objects."""
        # Restore universe from saved file paths
        if "_universe_topology" in state:
            topology = state.pop("_universe_topology")
            trajectory = state.pop("_universe_trajectory", None)
            frame = state.pop("_universe_frame", None)
            if topology and trajectory:
                try:
                    self.universe = mda.Universe(
                        topology, trajectory, **self._universe_kwargs
                    )
                    if frame is not None:
                        self.universe.trajectory[frame]
                except Exception as e:
                    raise RuntimeError(
                        f"Failed to restore Molecule from saved session. "
                        f"Could not recreate MDAnalysis Universe from topology '{topology}' "
                        f"and trajectory '{trajectory}'. "
                        f"The files may have been moved, deleted, or corrupted. "
                        f"Original error: {e}"
                    ) from e
        elif "_structure_filepath" in state or "_structure_code" in state:
            # universe was converted from a structure file; rebuild it the same way
            from ...converters import universe_from_atoms
            from ...download import StructureDownloader
            from .reader import read_structure

            frame = state.pop("_universe_frame", None)
            filepath = state.pop("_structure_filepath", None)
            code = state.pop("_structure_code", None)
            database = state.pop("_structure_database", "rcsb")
            data_block = state.pop("_structure_data_block", None) or None
            try:
                if filepath:
                    source = path_resolve(filepath)
                else:
                    source = StructureDownloader().download(
                        code=code, format="bcif", database=database or "rcsb"
                    )
                reader = read_structure(source, data_block=data_block)
                self.universe = universe_from_atoms(
                    reader.array, frame_data=reader.model_data()
                )
                if frame is not None:
                    self.universe.trajectory[frame]
            except Exception as e:
                raise RuntimeError(
                    f"Failed to restore Molecule from saved session. Could not rebuild "
                    f"the structure from source '{filepath or code}'. The file may have "
                    f"been moved, deleted, or corrupted. Original error: {e}"
                ) from e

        self.__dict__.update(state)

        # Recreate objects with circular references
        if not hasattr(self, "frame_manager"):
            self.frame_manager = FrameManager(self)
        if not hasattr(self, "selections"):
            self.selections = SelectionManager(self)
        if not hasattr(self, "annotations"):
            self.annotations = MoleculeAnnotationManager(self)
        if not hasattr(self, "calculations"):
            self.calculations = CalculationManager(self)
            preserved = self.__dict__.pop("_preserved_calculations", {})
            self.calculations._funcs.update(preserved)

    def _world_positions(self, selection: mda.AtomGroup) -> np.ndarray:
        """
        Atom positions in Blender world space.

        The universe holds positions in angstroms while the object is built at
        ``world_scale``, so a selection's positions have to be scaled and put
        through the object's transform to end up in the same space as the
        geometry drawn from them.
        """
        positions = np.asarray(selection.positions, dtype=np.float64)
        positions = positions * self.world_scale
        matrix = np.array(self.object.matrix_world)
        return positions @ matrix[:3, :3].T + matrix[:3, 3]

    def _view_points(self, selection: mda.AtomGroup | None) -> list[tuple]:
        """
        The world-space positions that make up a view.

        With no selection this is the geometry the molecule actually renders;
        with one it is the positions of the selected atoms.
        """
        if selection is None:
            points = blender_utils.evaluated_points(self.object)
        else:
            points = self._world_positions(selection)
        return [tuple(point) for point in points]

    def get_view(
        self, selection: str | AtomGroup | None = None, frame: int | None = None
    ) -> list[tuple]:
        """
        The world-space positions of a view onto this molecule.

        Without a selection these are the positions of the geometry the molecule
        renders; with one they are the positions of the selected atoms. Pass the
        result to [](`~mn.Canvas.look_at`) to frame it, or combine views from
        several selections with ``+`` to frame all of them together.

        Parameters
        ----------

        selection : str | AtomGroup, optional
            A selection phrase or AtomGroup
            When not specified, the whole entity is considered

        frame: int, optional
            Frame number of trajectory to use for calculating bounds.
            When not specified, current trajectory frame is used

        """
        if frame is not None:
            if frame < 0 or frame >= self.universe.trajectory.n_frames:
                raise ValueError(
                    f"{frame} is not within range [0, {self.universe.trajectory.n_frames - 1}]"
                )
        else:
            frame = self.uframe
        # temporarily set trajectory frame to specified value
        with temp_override_property(self, "uframe", frame):
            if selection is None:
                # return bbox of object when no selection specified
                return self._view_points(None)
            if isinstance(selection, AtomGroup):
                atom_group = selection
            elif isinstance(selection, str):
                # allow multiple comma separated selection phrases as well
                selection_array = selection.split(",")
                try:
                    atom_group = self.universe.select_atoms(*selection_array)
                except Exception:
                    raise ValueError(f"Invalid {selection} phrase")
            else:
                raise ValueError(f"{selection} is neither a str or AtomGroup")

            if atom_group.n_atoms == 0:
                # return bbox of object when selection is empty
                return self._view_points(None)

            # return the 3D bounding box vertices of the selected AtomGroup
            return self._view_points(atom_group)


def _on_add_assembly_instance(
    node: bpy.types.GeometryNodeGroup, obj: bpy.types.Object | None
) -> None:
    """
    Setup for an ``Assembly Instance`` node added from the asset library.

    Builds (or reuses) the assembly-transforms data object for the entity whose
    modifier tree received the node and assigns it to the node's Data Object
    input, so dropping the node onto a structure with assembly data just works.
    A no-op when a data object is already assigned or the entity has no
    parseable assembly data.
    """
    from ...session import get_session

    socket = node.inputs.get("Data Object")
    if socket is None or socket.default_value is not None or obj is None:
        return
    entity = get_session().match(obj)
    if not isinstance(entity, Molecule):
        return
    if entity.assemblies(as_array=True) is None:
        return
    socket.default_value = entity.create_data_object()


node_handlers.register_on_add(AssemblyInstance._name, _on_add_assembly_instance)
