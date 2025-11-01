"""
Trajectory entity for molecular dynamics simulations.

This module provides the Trajectory class for loading and visualizing molecular
dynamics trajectories in Blender using MDAnalysis as the backend.
"""

import functools
import inspect
import logging
from collections import OrderedDict
from pathlib import Path
from typing import Callable, Dict
import bpy
import databpy as db
import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from MDAnalysis.core.groups import AtomGroup
from ...assets import data
from ...blender import coll, path_resolve
from ...blender import utils as blender_utils
from ...nodes.geometry import (
    add_style_branch,
)
from ...nodes.nodes import styles_mapping
from ...nodes.styles import (
    StyleBase,
)
from ...utils import (
    frame_mapper,
    temp_override_property,
)
from ..base import EntityType, MolecularEntity
from .annotations import TrajectoryAnnotationManager
from .helpers import (
    AttributeMetadata,
    BlenderProperty,
    FrameManager,
    _validate_frame,
    _validate_non_negative,
)
from .selections import SelectionManager

logger = logging.getLogger(__name__)
# ============================================================================
# Utility Functions
# ============================================================================


def _unique_aname(obj: bpy.types.Object, prefix: str = "sel") -> str:
    attributes = db.list_attributes(obj)
    counter = 0
    aname = "{}_{}".format(prefix, counter)
    while aname in attributes:
        counter += 1
        aname = "{}_{}".format(prefix, counter)

    return aname


def _ag_to_bool(ag: mda.AtomGroup) -> np.ndarray:
    """Convert AtomGroup to boolean mask for the entire universe."""
    return np.isin(ag.universe.atoms.ix, ag.ix).astype(bool)


class Trajectory(MolecularEntity):
    """
    Molecular dynamics trajectory entity for Blender.

    This class provides a complete interface for loading, visualizing, and
    manipulating molecular dynamics trajectories in Blender using MDAnalysis
    as the backend. It handles:

    - Trajectory loading and frame management
    - Attribute computation (positions, elements, residues, etc.)
    - Selection management and updates
    - Visual styling and rendering
    - Frame interpolation and averaging
    - Periodic boundary condition handling
    - Integration with Blender's animation system

    Attributes:
        universe: MDAnalysis Universe containing the trajectory
        frame_manager: Handles position caching and frame updates
        selections: Manager for dynamic atom selections
        calculations: Dictionary of custom per-frame calculations
        annotations: Manager for trajectory annotations
        world_scale: Scale factor from Angstroms to Blender units
        frame_mapping: Optional custom frame number mapping

    Properties (synced with Blender):
        frame: Current animation frame
        subframes: Number of interpolation steps between frames
        offset: Frame offset for playback
        average: Number of frames to average (smoothing)
        correct_periodic: Apply periodic boundary corrections
        interpolate: Enable position interpolation between frames

    Example:
        >>> import MDAnalysis as mda
        >>> import molecularnodes as mn
        >>> u = mda.Universe("topology.pdb", "trajectory.dcd")
        >>> traj = mn.entities.Trajectory(u)
        >>> traj.add_style("ribbon", color="chain")
        >>> traj.selections.add("protein and name CA", name="alpha_carbons")
    """

    # Blender property descriptors with validation
    frame = BlenderProperty("frame", validate_fn=_validate_frame)
    subframes = BlenderProperty("subframes", validate_fn=_validate_non_negative)
    offset = BlenderProperty("offset")
    average = BlenderProperty("average", validate_fn=_validate_non_negative)
    correct_periodic = BlenderProperty("correct_periodic")
    interpolate = BlenderProperty("interpolate")
    _frame = BlenderProperty("frame_hidden")

    def __init__(
        self,
        universe: mda.Universe,
        name: str = "NewUniverseObject",
        world_scale: float = 0.01,
        create_object: bool = True,
    ):
        """
        Initialize a Trajectory entity from an MDAnalysis Universe.

        Args:
            universe: MDAnalysis Universe containing topology and trajectory
            name: Name for the Blender object
            world_scale: Scale factor from Angstroms to Blender units (default: 0.01)
            create_object: If True, immediately create the Blender object

        Note:
            The default world_scale of 0.01 converts Angstroms to Blender units,
            which works well for typical molecular sizes.
        """
        super().__init__()
        self.universe: mda.Universe = universe
        self.selections: SelectionManager = SelectionManager(self)
        self.calculations: Dict[str, Callable] = {}
        self.world_scale = world_scale
        self.frame_mapping: npt.NDArray[np.int64] | None = None
        self._entity_type = EntityType.MD
        self._updating_in_progress = False
        self.annotations = TrajectoryAnnotationManager(self)
        self.frame_manager = FrameManager(self)
        if create_object:
            self.create_object(name=name)

    @property
    def is_orthorhombic(self) -> bool:
        """
        Check if the simulation box has orthorhombic geometry.

        Orthorhombic boxes have all angles equal to 90 degrees, which
        enables optimized periodic boundary condition handling.

        Returns:
            True if box angles are all 90 degrees, False otherwise
        """
        dim = self.universe.dimensions
        if dim is None:
            return False
        return np.allclose(dim[3:], 90.0)

    @property
    def atoms(self) -> mda.AtomGroup:
        """All atoms in the trajectory as an MDAnalysis AtomGroup."""
        return self.universe.atoms

    @property
    def univ_positions(self) -> np.ndarray:
        """
        Current atom positions scaled to Blender world units.

        Returns:
            Array of shape (n_atoms, 3) with positions in Blender units
        """
        return self.atoms.positions * self.world_scale

    @property
    def bonds(self) -> np.ndarray | None:
        """
        Bond connectivity as pairs of atom indices.

        Returns:
            Array of shape (n_bonds, 2) with atom index pairs, or None if no bonds
        """
        if hasattr(self.atoms, "bonds"):
            return self.atoms.bonds.indices
        else:
            return None

    @property
    def n_frames(self) -> int:
        """Total number of frames in the trajectory."""
        return self.universe.trajectory.n_frames

    @property
    def cache(self) -> OrderedDict[int, np.ndarray]:
        """
        Access to the position cache for advanced use.

        The cache stores precomputed positions for recently accessed frames.
        This is primarily for internal use and testing.

        Returns:
            OrderedDict mapping frame numbers to position arrays
        """
        return self.frame_manager.cache._cache

    @property
    def uframe(self) -> int:
        """
        Current frame number in the MDAnalysis Universe.

        This differs from the `frame` property, which tracks the Blender
        scene frame. Use this to query the actual trajectory frame being displayed.

        Returns:
            Current Universe frame number (0-indexed)
        """
        return self.universe.trajectory.frame

    @uframe.setter
    def uframe(self, value: int) -> None:
        """
        Set the current frame in the MDAnalysis Universe.

        Only updates if the frame has changed to avoid redundant operations.

        Args:
            value: Target frame number (0-indexed, will be clamped to valid range)
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
            default_guesser = mda.guesser.default_guesser.DefaultGuesser(None)
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

    def _compute_atomic_number(self) -> np.ndarray:
        return np.array(
            [
                data.elements.get(element, data.elements.get("X")).get("atomic_number")
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
                data.residues.get(name, data.residues.get("UNK")).get("res_name_num")
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

    def _compute_segindices(
        self, metadata: AttributeMetadata | None = None
    ) -> np.ndarray:
        segs = []
        for seg in self.atoms.segments:
            segs.append(seg.atoms[0].segid)

        if metadata is not None:
            metadata.add("segments", segs)
        else:
            try:
                self.object["segments"] = segs
            except db.LinkedObjectError:
                logger.warning("Failed to store segments metadata on object")

        return self.atoms.segindices

    def _compute_chain_id_int(
        self, metadata: AttributeMetadata | None = None
    ) -> np.ndarray:
        chain_ids, chain_id_index = np.unique(self.atoms.chainIDs, return_inverse=True)

        if metadata is not None:
            metadata.add("chain_ids", chain_ids.astype(str).tolist())
        else:
            try:
                self.object["chain_ids"] = chain_ids.astype(str).tolist()
            except db.LinkedObjectError:
                logger.warning("Failed to store chain_ids metadata on object")

        return chain_id_index

    def _compute_atom_type_int(
        self, metadata: AttributeMetadata | None = None
    ) -> np.ndarray:
        atom_type_unique, atom_type_index = np.unique(
            self.atoms.types, return_inverse=True
        )

        if metadata is not None:
            metadata.add("atom_type_unique", atom_type_unique)
        else:
            try:
                self.object["atom_type_unique"] = atom_type_unique
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
        return np.isin(self.atoms.resnames, data.lipid_names)

    @property
    def _blender_attributes(self) -> Dict[str, Callable | str]:
        """
        Registry of default attributes to store on the Blender object.

        This property defines the standard molecular attributes that are
        computed and stored when a trajectory is created. Each entry maps
        an attribute name to either:
        - A callable that computes the attribute values
        - A selection string for MDAnalysis boolean selections

        Returns:
            Dictionary mapping attribute names to compute functions or selection strings
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

    def save_filepaths_on_object(self) -> None:
        """Save file paths to the Blender object for reference"""
        obj = self.object
        if isinstance(self.universe.filename, (str, Path)):
            obj.mn.filepath_topology = str(path_resolve(self.universe.filename))
        if isinstance(self.universe.trajectory.filename, (str, Path)):
            obj.mn.filepath_trajectory = str(
                path_resolve(self.universe.trajectory.filename)
            )

    def reset_playback(self) -> None:
        """Set the playback settings to their default values"""
        self.subframes = 0
        self.offset = 0
        self.average = 0
        self.correct_periodic = False
        self.interpolate = False

    def _store_default_attributes(self) -> None:
        """Store default attributes with batch metadata collection"""
        metadata = AttributeMetadata()

        for name, value in self._blender_attributes.items():
            try:
                if isinstance(value, str):
                    data = _ag_to_bool(self.universe.select_atoms(value))
                elif callable(value):
                    # Check if function accepts metadata parameter
                    sig = inspect.signature(value)
                    if "metadata" in sig.parameters:
                        data = value(metadata)
                    else:
                        data = value()
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

        # Batch apply all collected metadata to the object
        metadata.apply_to_object(self.object)

    def _store_extra_attributes(self) -> None:
        # TODO: enable adding of arbitrary mda.Universe attirbutes not currently applied
        pass

    def _create_object(
        self,
        name: str = "NewUniverseObject",
    ) -> None:
        """
        Internal method to create the Blender mesh object.

        Creates the mesh with initial positions and bonds, then stores
        all computed attributes and sets up modifiers.

        Args:
            name: Name for the Blender object
        """
        self.object = db.create_object(
            name=name,
            collection=coll.mn(),
            vertices=self.univ_positions,
            edges=self.bonds,
        )

        self._store_default_attributes()
        self._store_extra_attributes()
        self._setup_modifiers()

    def create_object(self, name: str = "NewUniverseObject") -> bpy.types.Object:
        """
        Create and initialize the Blender object for this trajectory.

        This method:
        1. Creates the mesh object with positions and bonds
        2. Computes and stores all molecular attributes
        3. Sets up geometry node modifiers
        4. Registers the object with MolecularNodes
        5. Sets it as the active object

        Args:
            name: Name for the Blender object

        Returns:
            The created Blender object
        """
        self._create_object(name=name)

        self.object.mn.entity_type = self._entity_type.value
        self.object.mn.n_frames = self.n_frames
        self.save_filepaths_on_object()
        bpy.context.view_layer.objects.active = self.object

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

    def _frame_range(self, frame: int) -> npt.NDArray[np.int64]:
        """
        Get frame numbers to include when averaging.

        Args:
            frame: Center frame number

        Returns:
            Array of frame numbers to average over
        """
        return self.frame_manager._frame_range(frame)

    def _cache_ordered(self) -> np.ndarray:
        """
        Get cached frames as a chronologically ordered 3D array.

        Returns:
            Array of shape (n_cached_frames, n_atoms, 3)
        """
        return self.frame_manager.cache.get_ordered_array()

    def adjust_periodic_positions(
        self, pos1: np.ndarray, pos2: np.ndarray
    ) -> np.ndarray:
        """
        Apply periodic boundary corrections to positions.

        Corrects for atoms that crossed periodic boundaries between pos1 and pos2,
        ensuring smooth visualization across boundaries.

        Args:
            pos1: Reference positions
            pos2: Positions to potentially correct

        Returns:
            Corrected positions (or unchanged if correction not needed/applicable)
        """
        return self.frame_manager.adjust_periodic_positions(pos1, pos2)

    def position_cache_mean(self, frame: int) -> np.ndarray:
        """
        Get mean position from the averaging window.

        Computes the average over frames specified by the `average` property,
        applying periodic corrections if enabled.

        Args:
            frame: Center frame for averaging window

        Returns:
            Averaged atom positions
        """
        return self.frame_manager.position_cache_mean(frame)

    def set_frame(self, frame: int) -> None:
        """
        Update trajectory state for a given scene frame.

        This is the main entry point called by Blender's animation system.
        It updates positions, selections, and custom calculations while
        preventing recursive updates.

        Args:
            frame: Scene frame number (not Universe frame - mapping is applied)

        Note:
            This method is typically called automatically by frame change handlers,
            not directly by user code.
        """
        if self._updating_in_progress:
            logger.debug("Update already in progress, skipping nested update")
            return

        try:
            self._updating_in_progress = True
            self._update_positions(frame)
            self._update_selections()
            self._update_calculations()
        finally:
            self._updating_in_progress = False

    def _position_at_frame(self, frame: int) -> np.ndarray:
        """
        Get raw atom positions at a specific Universe frame.

        Args:
            frame: Universe frame number

        Returns:
            Scaled atom positions
        """
        return self.frame_manager._position_at_frame(frame)

    def update_position_cache(self, frame: int, cache_ahead: bool = True) -> None:
        """
        Update the position cache for the current frame.

        Intelligently caches positions based on averaging and interpolation settings.

        Args:
            frame: Current frame number
            cache_ahead: If True, prefetch next frame for interpolation
        """
        self.frame_manager.update_position_cache(frame, cache_ahead)

    def frame_mapper(self, frame: int) -> int:
        """
        Map scene frame to Universe frame number.

        Applies offset, subframe subdivision, and custom frame mapping.

        Args:
            frame: Scene frame number

        Returns:
            Corresponding Universe frame number
        """
        return frame_mapper(
            frame=frame,
            subframes=self.subframes,
            offset=self.offset,
            mapping=self.frame_mapping,
        )

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
        self.object.mn.styles_active_index = self.tree.nodes.find(node_style.name)

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
            if frame < 0 or frame >= self.n_frames:
                raise ValueError(
                    f"{frame} is not within range [0, {self.n_frames - 1}]"
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
