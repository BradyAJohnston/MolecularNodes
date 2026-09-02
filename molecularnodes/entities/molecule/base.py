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
from typing import Callable, Dict, Sequence
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
from ...nodes.material import PresetMaterial, append_material
from ...nodes.nodes import STYLE_LITERALS, STYLE_NODE_MAPPING
from ...utils import _UNSET, count_value_changes, temp_override_property
from ..base import EntityType, MolecularEntity
from ..utilities import (
    BoolObjectMNProperty,
    IntObjectMNProperty,
    StringObjectMNProperty,
    _validate_non_negative,
)
from .annotations import MoleculeAnnotationManager
from .dssp import DSSPManager
from .helpers import FrameManager, _ag_to_bool
from .selections import SelectionManager

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
    calculations : dict
        Custom per-frame calculations
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
    dssp : DSSPManager
        A DSSP Manager to compute and show secondary structures

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
    _entity_type: EntityType = EntityType.MD

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
        self.calculations: Dict[str, Callable] = {}
        self.world_scale = world_scale
        self._updating_in_progress = False
        self.annotations = MoleculeAnnotationManager(self)
        self.frame_manager = FrameManager(self)
        self.dssp = DSSPManager(self)
        if create_object:
            self.create_object(name=name)

    @property
    def _is_orthorhombic(self) -> bool:
        """Check if simulation box is orthorhombic.

        Orthorhombic boxes (all angles = 90°) enable optimized periodic boundary handling.

        Returns
        -------
        bool
            True if box angles are all 90 degrees
        """
        dim = self.universe.dimensions
        if dim is None:
            return False
        return np.allclose(dim[3:], 90.0)

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

    def _compute_elements(self) -> np.ndarray:
        """Return cached elements (for backwards compatibility)"""
        return self._elements

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
            return np.array([x.mass for x in self.atoms])
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
        # corresponds to Charges topology attr
        if hasattr(self.atoms, "charges"):
            return self.atoms.charges
        else:
            return np.zeros(len(self))

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
        return np.isin(self.atoms.resnames, data.RESNAMES_LIPID)

    def _compute_is_solvent(self) -> np.ndarray:
        resname_is_solvent = np.isin(self.atoms.resnames, data.RESNAMES_SOLVENT)
        name_is_solvent = np.isin(self.atoms.names, data.NAMES_SOLVENT)
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

    def _sel_bool(self, selection: str) -> np.ndarray:
        """Evaluate an MDAnalysis selection string to a per-atom boolean mask."""
        return _ag_to_bool(self.universe.select_atoms(selection))

    def _compute_is_backbone(self) -> np.ndarray:
        is_backbone_atom = np.isin(self.atoms.names, self._BACKBONE_ATOM_NAMES)
        return np.logical_and(is_backbone_atom, ~self._compute_is_solvent())

    def _compute_is_side_chain(self) -> np.ndarray:
        # side chain = polymer atoms that are not backbone, but the alpha carbon
        # (or CG backbone bead) is counted as side chain. Mirrors the biotite reader.
        backbone = self._compute_is_backbone()
        is_alpha_carbon = np.isin(self.atoms.names, ("CA", "BB"))
        is_polymer = np.logical_or(
            self._sel_bool("protein or (name BB SC*)"), self._sel_bool("nucleic")
        )
        return np.logical_and(np.logical_or(~backbone, is_alpha_carbon), is_polymer)

    def _compute_lipophobicity(self) -> np.ndarray:
        return np.array(
            [
                data.lipophobicity.get(res, {}).get(atom, 0)
                for res, atom in zip(self.atoms.resnames, self.atoms.names)
            ],
            dtype=float,
        )

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
        # carried across from the structure file by the converter
        return self.atoms.sec_structs

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
            "lipophobicity": self._compute_lipophobicity,
            "Color": self._compute_color,
            "is_alpha_carbon": "name CA or name BB",
            "is_backbone": self._compute_is_backbone,
            "is_side_chain": self._compute_is_side_chain,
            "is_solvent": self._compute_is_solvent,
            "is_nucleic": "nucleic",
            "is_lipid": self._compute_is_lipid,
            "is_peptide": "protein or (name BB SC*)",
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
                    data = _ag_to_bool(self.universe.select_atoms(item))
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
        pass

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
        **kwargs,
    ) -> "Molecule":
        """Load a single structure file, or an MD topology + trajectory.

        With only ``topology`` given, it is treated as a single structure file
        (``.pdb``/``.cif``/``.bcif``/``.sdf``/``.mol``) and routed through the biotite
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
        kwargs : dict, optional
            Additional keyword arguments to pass to the `MDAnalysis.Universe()` constructor.

        Returns
        -------
        Molecule
            The loaded entity.
        """
        if coordinates is None:
            entity = cls.from_file(topology, name=name)
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
    ) -> "Molecule":
        """Load a single structure file into a Universe-backed entity.

        The file (``.pdb``/``.cif``/``.bcif``/``.sdf``/``.mol``) is parsed by the biotite
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

        reader = read_structure(file_path)
        universe = universe_from_atoms(reader.array)
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
            The database accession code (e.g. a PDB id).
        format : str, optional
            File format to download, by default ``".bcif"``.
        cache : Path | str | None, optional
            Directory to cache downloads in.
        database : str, optional
            The database to fetch from, by default ``"rcsb"``.

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
        # a structure loaded from a single file is a "molecule" rather than an MD
        # trajectory; whether playback UI is shown is decided by the frame count.
        self.props.entity_type = EntityType.MOLECULE.value
        self.props.entity_ids = reader.entity_ids()
        self.props.chain_ids = reader.chain_ids()
        self.props.biological_assemblies = reader.assemblies(as_json_string=True)
        # record the source so the entity can be reloaded into a fresh session
        if not isinstance(file_path, io.BytesIO):
            self.props.filepath = str(file_path)

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
            except (ValueError, TypeError):
                return None
        return assemblies_info

    def create_data_object(self) -> bpy.types.Object:
        """Create the data object holding the biological assembly transforms."""
        from ...blender import mesh

        data_obj_name = f".data_{self.name}_assemblies"
        data_obj = bpy.data.objects.get(data_obj_name)
        if not data_obj:
            transforms = self.assemblies(as_array=True)
            data_obj = mesh.create_data_object(array=transforms, name=data_obj_name)

        return data_obj

    def _update_calculations(self) -> None:
        """Update all registered calculations for the current frame"""
        for name, func in self.calculations.items():
            try:
                self.store_named_attribute(data=func(self.universe), name=name)
            except Exception as e:
                logger.error(
                    f"Failed to update calculation '{name}': {e}", exc_info=True
                )

    def _update_selections(self) -> None:
        """Update all selections for the current frame."""
        self.selections.update_attributes()

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
            self._update_positions(frame)
            self._update_selections()
            self._update_calculations()
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
        units) for a single frame, in the same vertex order as this molecule's mesh. The
        collection can then be read inside geometry nodes (e.g. with the *Animate Frames*
        node) to drive positions from baked data, rather than updating the ``Universe``
        every scene frame.

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

        # bake each requested frame, restoring the current frame afterwards
        original_uframe = self.uframe
        try:
            for i, frame in enumerate(range(start, stop, step)):
                self.uframe = frame
                db.create_object(
                    vertices=self._scaled_position,
                    name=f"{self.name}_frame_{i:04d}",
                    collection=frames,
                )
        finally:
            self.uframe = original_uframe

        return frames

    def _update_positions(self, frame: int) -> None:
        """
        Internal method to update atom positions.

        Delegates to FrameManager which handles caching, interpolation,
        averaging, and periodic corrections.

        Args:
            frame: Scene frame number
        """
        self._update_trajectory_positions(frame)

    def _update_trajectory_positions(self, frame: int) -> None:
        """Update trajectory positions for the given frame.

        This method can be overridden by subclasses to implement custom
        position update logic (e.g., for streaming trajectories).

        Parameters
        ----------
        frame : int
            Scene frame number
        """
        self.position = self.frame_manager.get_positions_at_frame(frame)

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

    def _style_color_input(self, color: str | Sequence[float] | Callable):
        """Resolve `add_style`'s color argument into the `Set Color` node input.

        Returns ``None`` if no ``Set Color`` node should be added, which happens when
        the argument names neither a known keyword nor an existing attribute.
        """
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
        warnings.warn(
            f"Color '{color}' is neither a known color keyword {COLOR_KEYWORDS} nor a "
            f"color attribute on this molecule (available: {color_attributes}). No "
            "color will be applied. For coloring beyond the keywords, pass a callable "
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
        color: str | Sequence[float] | Callable | None = None,
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
            Default is "MN Default".

        color : str | Sequence[float] | Callable | None, optional
            Coloring to apply upstream of the style via a ``Set Color`` node. Can be:
            - ``"common"`` / ``"default"`` for standard element colors with carbons
              colored randomly per chain
            - ``"plddt"`` to color by pLDDT (B-factor) confidence
            - the name of an existing *color* attribute on the geometry
            - an RGBA sequence of floats for a single uniform color
            - a zero-argument callable returning a color socket, evaluated inside the
              node tree context, e.g. ``lambda: mg.ColorRainbow()``. This is the way
              to reach the full set of ``Color*`` nodes without writing out a whole
              node tree.
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
        from . import OXDNA

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

        material = "MN Default" if material is _UNSET else material
        # a callable selection is evaluated inside the tree context, below
        selection_is_callable = callable(selection)
        attribute_name = self._resolve_style_selection(selection)

        if isinstance(self, OXDNA):
            STYLE_NODE_MAPPING["ribbon"] = g.OxDNAStyleRibbon  # ty: ignore[invalid-assignment]

        if not style_is_callable and "sphere" not in kwargs:
            # spheres default to point clouds, which only Cycles can draw
            geometry = _sphere_for_engine(
                STYLE_NODE_MAPPING[style], bpy.context.scene.render.engine
            )
            if geometry is not None:
                kwargs["sphere"] = geometry

        if isinstance(material, (PresetMaterial, MaterialBuilder)):
            material = material.material

        material = append_material(material) if isinstance(material, str) else material

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
                style_node = STYLE_NODE_MAPPING[style](
                    selection=selection_input,
                    material=material,
                    **kwargs,
                )
            if name:
                style_node.node.label = name
            color_input = self._style_color_input(color) if color is not None else None
            (
                tree.atoms
                >> (g.SetColor(color=color_input) if color_input is not None else None)
                >> style_node
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
        if "dssp" in state:
            del state["dssp"]
        if "calculations" in state:
            # Preserve picklable calculations
            preserved_calculations = {}
            for name, calc_func in state["calculations"].items():
                if name == "sec_struct":
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
                    self.universe = mda.Universe(topology, trajectory)
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
            try:
                if filepath:
                    source = path_resolve(filepath)
                else:
                    source = StructureDownloader().download(
                        code=code, format="bcif", database=database or "rcsb"
                    )
                self.universe = universe_from_atoms(read_structure(source).array)
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
        if not hasattr(self, "dssp"):
            self.dssp = DSSPManager(self)
        if not hasattr(self, "calculations"):
            self.calculations = state.pop("_preserved_calculations", {})

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
