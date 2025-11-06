"""Trajectory entity for molecular dynamics simulations.

Provides Trajectory class for loading and visualizing MD trajectories in Blender
using MDAnalysis.
"""

import functools
import logging
from pathlib import Path
from typing import Callable, Dict
import bpy
import databpy as db
import MDAnalysis as mda
import numpy as np
from MDAnalysis.core.groups import AtomGroup
from ...assets import data
from ...blender import coll, path_resolve, set_obj_active
from ...blender import utils as blender_utils
from ...nodes.geometry import (
    add_style_branch,
)
from ...nodes.nodes import styles_mapping
from ...nodes.styles import (
    StyleBase,
)
from ...utils import (
    temp_override_property,
)
from ..base import EntityType, MolecularEntity
from ..utilities import (
    BoolObjectMNProperty,
    IntObjectMNProperty,
    StringObjectMNProperty,
    _unique_aname,
    _validate_non_negative,
)
from .annotations import TrajectoryAnnotationManager
from .helpers import FrameManager, _ag_to_bool
from .selections import SelectionManager

logger = logging.getLogger(__name__)


class Trajectory(MolecularEntity):
    """MD trajectory entity for Blender visualization.

    Complete interface for loading, visualizing, and manipulating MD trajectories
    in Blender using MDAnalysis.

    Features: trajectory loading, attribute computation, selection management,
    visual styling, frame interpolation/averaging, periodic boundary handling,
    Blender animation integration.

    Attributes
    ----------
    universe : mda.Universe
        MDAnalysis Universe with topology and trajectory
    frame_manager : FrameManager
        Position caching and frame updates
    selections : SelectionManager
        Dynamic atom selections
    calculations : dict
        Custom per-frame calculations
    annotations : TrajectoryAnnotationManager
        Trajectory annotations
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
    traj = mn.entities.Trajectory(u)
    traj.add_style(mn.StyleSpheres(geometry="Mesh"), selection="resname LYS")
    canvas.frame_view(traj)
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
    _mn_styles_active_index = IntObjectMNProperty(
        "styles_active_index", _validate_non_negative
    )
    _mn_entity_type = StringObjectMNProperty("entity_type")
    _mn_filepath_topology = StringObjectMNProperty("filepath_topology")
    _mn_filepath_trajectory = StringObjectMNProperty("filepath_trajectory")
    _mn_n_frames = IntObjectMNProperty("n_frames", _validate_non_negative)

    def __init__(
        self,
        universe: mda.Universe,
        name: str = "NewUniverseObject",
        world_scale: float = 0.01,
        create_object: bool = True,
    ):
        """Initialize Trajectory from MDAnalysis Universe.

        Parameters
        ----------
        universe : mda.Universe
            MDAnalysis Universe with topology and trajectory
        name : str, default="NewUniverseObject"
            Name for the Blender object
        world_scale : float, default=0.01
            Scale factor from Angstroms to Blender units
        create_object : bool, default=True
            Whether to immediately create the Blender object

        Notes
        -----
        Default world_scale of 0.01 converts Angstroms to Blender units.
        """
        super().__init__()
        self.universe: mda.Universe = universe
        self.selections: SelectionManager = SelectionManager(self)
        self.calculations: Dict[str, Callable] = {}
        self.world_scale = world_scale
        self._updating_in_progress = False
        self.annotations = TrajectoryAnnotationManager(self)
        self.frame_manager = FrameManager(self)
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
                for x in self.atoms.names  # type: ignore
            ]
            return np.array(guessed_elements)
        except Exception as e:
            logger.warning(f"Failed to compute elements, using placeholder 'X': {e}")
            return np.repeat("X", len(self))

    def _compute_elements(self) -> np.ndarray:
        """Return cached elements (for backwards compatibility)"""
        return self._elements

    def _compute_atomic_number(self) -> np.ndarray:
        return np.array(
            [
                data.elements.get(element, data.elements["X"]).get("atomic_number")
                for element in self._elements
            ]
        )

    def _compute_vdw_radii(self) -> np.ndarray:
        return (
            np.array(
                [
                    data.elements.get(element, {}).get("vdw_radii", 100)
                    for element in self._elements
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

    def _compute_atom_id(self) -> np.ndarray:
        return self.atoms.ids

    def _compute_segindices(self) -> np.ndarray:
        segs = []
        for seg in self.atoms.segments:
            segs.append(seg.atoms[0].segid)

        else:
            try:
                self.object["segments"] = segs
            except db.LinkedObjectError:
                logger.warning("Failed to store segments metadata on object")

        return self.atoms.segindices

    def _compute_chain_id_int(self) -> np.ndarray:
        chain_ids, chain_id_index = np.unique(self.atoms.chainIDs, return_inverse=True)

        try:
            self.object["chain_ids"] = chain_ids.astype(str).tolist()
        except db.LinkedObjectError:
            logger.warning("Failed to store chain_ids metadata on object")

        return chain_id_index

    def _compute_atom_type_int(self) -> np.ndarray:
        atom_type_unique, atom_type_index = np.unique(
            self.atoms.types, return_inverse=True
        )

        try:
            self.object["atom_type_unique"] = atom_type_unique
        except db.LinkedObjectError:
            logger.warning("Failed to store atom_type_unique metadata on object")

        return atom_type_index

    def _compute_atom_name_int(self) -> np.ndarray:
        if hasattr(self.atoms, "names"):
            return np.array(
                [data.atom_names.get(x, -1) for x in self.atoms.names],  # type: ignore
                dtype=int,
            )
        else:
            return np.repeat(int(-1), len(self))

    def _compute_is_lipid(self) -> np.ndarray:
        return np.isin(self.atoms.resnames, data.lipid_names)

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
            "segid": self._compute_segindices,
            "res_name": self._compute_res_name_int,
            "atom_id": self._compute_atom_id,
            "b_factor": self._compute_b_factor,
            "occupancy": self._compute_occupancy,
            "charge": self._compute_charge,
            "chain_id": self._compute_chain_id_int,
            "atom_types": self._compute_atom_type_int,
            "atom_name": self._compute_atom_name_int,
            "is_backbone": "backbone or nucleicbackbone or name BB",
            "is_alpha_carbon": "name CA or name BB",
            "is_solvent": "name OW or name HW1 or name HW2 or resname W or resname PW",
            "is_nucleic": "nucleic",
            "is_lipid": self._compute_is_lipid,
            "is_peptide": "protein or (name BB SC*)",
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
                self.store_named_attribute(
                    data=data,
                    name=name,
                )
            except (mda.NoDataError, AttributeError) as e:
                logger.debug(f"Skipping attribute '{name}': {e}")
            except Exception as e:
                logger.warning(f"Failed to compute attribute '{name}': {e}")

    def _store_extra_attributes(self) -> None:
        # TODO: enable adding of arbitrary mda.Universe attirbutes not currently applied
        pass

    def _create_object(
        self,
        name: str = "NewUniverseObject",
    ) -> None:
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
        self._mn_entity_type = EntityType.MD.value

        self._store_default_attributes()
        self._store_extra_attributes()
        self._setup_modifiers()

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

        self._mn_n_frames = self.universe.trajectory.n_frames
        self._save_filepaths_on_object()
        set_obj_active(self.object)

        return self.object

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
        for item in self.selections.items:
            try:
                # Skip non-updating selections - they use static masks
                if not item.updating:
                    continue

                # Lazy initialization will occur if needed
                selection = self.selections.get(item.name)
                if selection is None:
                    raise KeyError(f"Selection '{item.name}' not found")

                # Don't recreate atomgroup for immutable selections (created from atomgroups)
                # These selections maintain their own atomgroup reference
                if not item.immutable:
                    selection.set_atom_group(item.string)
                selection.set_selection()
            except KeyError as e:
                logger.warning(
                    f"Failed to update selection '{item.name}': {e}. Skipping this selection."
                )
            except mda.SelectionError as e:
                logger.error(
                    f"Invalid selection syntax for '{item.name}': {e}. Skipping this selection."
                )
            except Exception as e:
                logger.error(
                    f"Error updating selection '{item.name}': {e}. Skipping this selection.",
                    exc_info=True,
                )

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
        finally:
            self._updating_in_progress = False

    def _update_positions(self, frame: int) -> None:
        """
        Internal method to update atom positions.

        Delegates to FrameManager which handles caching, interpolation,
        averaging, and periodic corrections.

        Args:
            frame: Scene frame number
        """
        self.position = self.frame_manager.get_positions_at_frame(frame)

    def __repr__(self) -> str:
        return f"<Trajectory, `universe`: {self.universe}, `object`: {self.object}"

    def add_style(
        self,
        style: StyleBase | str = "spheres",
        color: str | None = "common",
        selection: str | AtomGroup | None = None,
        material: bpy.types.Material | str | None = None,
        name: str | None = None,
    ) -> "Trajectory":
        """
        Add a visual style to the trajectory.

        Parameters
        ----------
        style : bpy.types.GeometryNodeTree | str, optional
            The style to apply to the trajectory. Can be a GeometryNodeTree or a string
            identifying a predefined style (e.g., "spheres", "sticks", "ball_stick").
            Default is "spheres".

        color : str | None, optional
            The coloring scheme to apply. Can be "common" (element-based coloring),
            "chain", "residue", or other supported schemes. If None, no coloring
            is applied. Default is "common".

        selection : str | AtomGroup | None, optional
            Apply the style only to atoms matching this selection. Can be:
            - A string referring to an existing boolean attribute on the trajectory
            - A AtomGroup object defining a selection criteria
            - None to apply to all atoms (default)

        material : bpy.types.Material | str | None, optional
            The material to apply to the styled atoms. Can be a Blender Material object,
            a string with a material name, or None to use default materials. Default is None.

        name: str, optional
            The label for this style

        Returns
        -------
        Trajectory
            Returns self for method chaining.

        Raises
        ------
        ValueError
            If an unsupported style string is passed

        Notes
        -----
        If a selection is provided, it will be evaluated and stored as a new
        named attribute on the trajectory with an automatically generated name (sel_N).
        """
        if style is None:
            return self

        if isinstance(style, str) and style not in styles_mapping:
            raise ValueError(
                f"Invalid style '{style}'. Supported styles are {[key for key in styles_mapping.keys()]}"
            )
        if selection is None:
            sel_name = None
        else:
            att_name = _unique_aname(self.object, "sel")
            if isinstance(selection, str):
                # TODO: There are currently no validations for the selection phrase
                sel_name = self.selections.add(selection, att_name).name
            elif isinstance(selection, AtomGroup):
                sel_name = self.selections.from_atomgroup(selection, name=att_name).name
            # TODO: Delete these named attributes when style is deleted
            # Currently, styles are removed using GeometryNodeInterFace.remove(),

        node_style = add_style_branch(
            tree=self.tree,
            style=style,
            color=color,
            selection=sel_name,
            material=material,
            name=name,
        )

        # set the active index for UI to the newly added style
        self._mn_styles_active_index = self.tree.nodes.find(node_style.name)

        return self

    def _get_3d_bbox(self, selection: mda.AtomGroup | None) -> list[tuple]:
        """Get the 3D bounding box vertices of atoms in an AtomGroup"""
        if selection is None:
            return blender_utils.get_bounding_box(self.object)
        v0, v1 = selection.bbox() * self.world_scale
        bb_verts_3d = [
            (v0[0], v0[1], v0[2]),
            (v0[0], v0[1], v1[2]),
            (v0[0], v1[1], v0[2]),
            (v0[0], v1[1], v1[2]),
            (v1[0], v0[1], v0[2]),
            (v1[0], v0[1], v1[2]),
            (v1[0], v1[1], v0[2]),
            (v1[0], v1[1], v1[2]),
        ]
        return bb_verts_3d

    def get_view(
        self, selection: str | AtomGroup | None = None, frame: int | None = None
    ) -> list[tuple]:
        """
        Get the 3D bounding box of a selection within the trajectory

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
                return self._get_3d_bbox(None)
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
                return self._get_3d_bbox(None)

            # return the 3D bounding box vertices of the selected AtomGroup
            return self._get_3d_bbox(atom_group)
