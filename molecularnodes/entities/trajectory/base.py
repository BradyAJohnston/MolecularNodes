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
    correct_periodic_positions,
    fraction,
    frame_mapper,
    frames_to_average,
    temp_override_property,
)
from ..base import EntityType, MolecularEntity
from .annotations import TrajectoryAnnotationManager
from .selections import SelectionManager


def _unique_aname(obj: bpy.types.Object, prefix: str = "sel") -> str:
    attributes = db.list_attributes(obj)
    counter = 0
    aname = "{}_{}".format(prefix, counter)
    while aname in attributes:
        counter += 1
        aname = "{}_{}".format(prefix, counter)

    return aname


class Trajectory(MolecularEntity):
    def __init__(
        self,
        universe: mda.Universe,
        name: str = "NewUniverseObject",
        world_scale: float = 0.01,
        create_object: bool = True,
    ):
        super().__init__()
        self.universe: mda.Universe = universe
        self.selections: SelectionManager = SelectionManager(self)
        self.calculations: Dict[str, Callable] = {}
        self.world_scale = world_scale
        self.frame_mapping: npt.NDArray[np.int64] | None = None
        self.cache: dict = {}
        self._entity_type = EntityType.MD
        self._updating_in_progress = False
        self.annotations = TrajectoryAnnotationManager(self)
        if create_object:
            self.create_object(name=name)

    @property
    def is_orthorhombic(self):
        dim = self.universe.dimensions
        if dim is None:
            return False

        return np.allclose(dim[3:], 90.0)

    @property
    def atoms(self) -> mda.AtomGroup:
        return self.universe.atoms

    @staticmethod
    def bool_selection(ag, selection, **kwargs) -> np.ndarray:
        return np.isin(ag.ix, ag.select_atoms(selection, **kwargs).ix).astype(bool)

    @property
    def univ_positions(self) -> np.ndarray:
        return self.atoms.positions * self.world_scale

    @property
    def bonds(self) -> np.ndarray | None:
        if hasattr(self.atoms, "bonds"):
            return self.atoms.bonds.indices
        else:
            return None

    def _compute_elements(self) -> np.ndarray:
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

        except Exception:
            return np.repeat("X", len(self))

    def _compute_atomic_number(self) -> np.ndarray:
        return np.array(
            [
                data.elements.get(element, data.elements.get("X")).get("atomic_number")
                for element in self._compute_elements()
            ]
        )

    def _compute_vdw_radii(self) -> np.ndarray:
        return (
            np.array(
                [
                    data.elements.get(element, {}).get("vdw_radii", 100)
                    for element in self._compute_elements()
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
                for element in self._compute_elements()
            ]
            return np.array(masses)

    @property
    def n_frames(self) -> int:
        return self.universe.trajectory.n_frames

    @property
    def uframe(self) -> int:
        """
        Get the current frame number of the linked `Universe.trajectory`.

        Returns:
            int: Current frame number in the trajectory.
        """
        return self.universe.trajectory.frame

    @uframe.setter
    def uframe(self, value) -> None:
        """
        Set the current frame number of the linked `Universe.trajectory`.

        The frame number is clamped between 0 and n_frames-1 to prevent
        out-of-bounds access.

        Args:
            value (int): Target frame number to set.

        Returns:
            None
        """
        if self.universe.trajectory.frame != value:
            self.universe.trajectory[value]

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

    def _compute_segindices(self) -> np.ndarray:
        segs = []
        for seg in self.atoms.segments:
            segs.append(seg.atoms[0].segid)

        try:
            self.object["segments"] = segs
        except db.LinkedObjectError:
            pass

        return self.atoms.segindices

    def _compute_chain_id_int(self) -> np.ndarray:
        chain_ids, chain_id_index = np.unique(self.atoms.chainIDs, return_inverse=True)

        self.object["chain_ids"] = chain_ids.astype(str).tolist()

        return chain_id_index

    def _compute_atom_type_int(self) -> np.ndarray:
        atom_type_unique, atom_type_index = np.unique(
            self.atoms.types, return_inverse=True
        )
        try:
            self.object["atom_type_unique"] = atom_type_unique
        except db.LinkedObjectError:
            pass

        return atom_type_index

    def _compute_atom_name_int(self) -> np.ndarray:
        if hasattr(self.atoms, "names"):
            return np.array(
                [data.atom_names.get(x, -1) for x in self.atoms.names],
                dtype=int,
            )
        else:
            return np.repeat(int(-1), len(self))

    def _compute_is_nucleic(self) -> np.ndarray:
        return self.bool_selection(self.atoms, "nucleic")

    def _compute_is_peptide(self) -> np.ndarray:
        return self.bool_selection(self.atoms, "protein or (name BB SC*)")

    def _compute_is_lipid(self) -> np.ndarray:
        return np.isin(self.atoms.resnames, data.lipid_names)

    def _compute_is_backbone(self) -> np.ndarray:
        return self.bool_selection(self.atoms, "backbone or nucleicbackbone or name BB")

    def _compute_is_alpha_carbon(self) -> np.ndarray:
        return self.bool_selection(self.atoms, "name CA or name BB")

    def _compute_is_solvent(self) -> np.ndarray:
        return self.bool_selection(
            self.atoms, "name OW or name HW1 or name HW2 or resname W or resname PW"
        )

    @property
    def _blender_attributes(self):
        """
        Attribute names and the methods to generate their values for the Blender object
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
            "is_backbone": self._compute_is_backbone,
            "is_alpha_carbon": self._compute_is_alpha_carbon,
            "is_solvent": self._compute_is_solvent,
            "is_nucleic": self._compute_is_nucleic,
            "is_lipid": self._compute_is_lipid,
            "is_peptide": self._compute_is_peptide,
        }

    def save_filepaths_on_object(self) -> None:
        obj = self.object
        if isinstance(self.universe.filename, (str, Path)):
            obj.mn.filepath_topology = str(path_resolve(self.universe.filename))
        if isinstance(self.universe.trajectory.filename, (str, Path)):
            obj.mn.filepath_trajectory = str(
                path_resolve(self.universe.trajectory.filename)
            )

    def reset_playback(self) -> None:
        "Set the playback settings to their default values"
        self.subframes = 0
        self.offset = 0
        self.average = 0
        self.correct_periodic = False
        self.interpolate = False

    def _store_default_attributes(self) -> None:
        for name, func in self._blender_attributes.items():
            try:
                self.store_named_attribute(
                    data=func(),
                    name=name,
                )
            except (mda.NoDataError, AttributeError):
                pass

    def _store_extra_attributes(self) -> None:
        # TODO: enable adding of arbitrary mda.Universe attirbutes not currently applied
        pass

    def _create_object(
        self,
        name: str = "NewUniverseObject",
    ) -> None:
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
        self._create_object(name=name)

        self.object.mn.entity_type = self._entity_type.value
        self.object.mn.n_frames = self.n_frames
        self.save_filepaths_on_object()
        bpy.context.view_layer.objects.active = self.object

        return self.object

    def _update_calculations(self):
        for name, func in self.calculations.items():
            try:
                self.store_named_attribute(data=func(self.universe), name=name)
            except Exception as e:
                print(e)

    def _update_selections(self):
        """Update all selections for the current frame."""
        for item in self.selections.items:
            try:
                # Skip non-updating selections - they use static masks
                if not item.updating:
                    continue

                # Lazy initialization will occur if needed
                selection = self.selections.get(item.name)
                if selection is None:
                    raise KeyError
                # Don't recreate atomgroup for immutable selections (created from atomgroups)
                # These selections maintain their own atomgroup reference
                if not item.immutable:
                    selection.set_atom_group(item.string)
                selection.set_selection()
            except KeyError as e:
                print(
                    f"Warning: Failed to update selection '{item.name}': {e}. "
                    "Skipping this selection."
                )
            except Exception as e:
                print(
                    f"Warning: Error updating selection '{item.name}': {e}. "
                    "Skipping this selection."
                )

    @property
    def _frame(self) -> int:
        return self.object.mn.frame_hidden

    @_frame.setter
    def _frame(self, value: int) -> None:
        self.object.mn.frame_hidden = value

    @property
    def frame(self) -> int:
        return self.object.mn.frame

    @frame.setter
    def frame(self, value: int) -> None:
        self.object.mn.frame = value

    @property
    def subframes(self) -> int:
        return self.object.mn.subframes

    @subframes.setter
    def subframes(self, value: int) -> None:
        self.object.mn.subframes = value

    @property
    def offset(self) -> int:
        return self.object.mn.offset

    @offset.setter
    def offset(self, value: int) -> None:
        self.object.mn.offset = value

    @property
    def average(self) -> int:
        return self.object.mn.average

    @average.setter
    def average(self, value: int) -> None:
        self.object.mn.average = value

    @property
    def correct_periodic(self) -> bool:
        return self.object.mn.correct_periodic

    @correct_periodic.setter
    def correct_periodic(self, value: bool) -> None:
        self.object.mn.correct_periodic = value

    @property
    def interpolate(self) -> bool:
        return self.object.mn.interpolate

    @interpolate.setter
    def interpolate(self, value: bool) -> None:
        self.object.mn.interpolate = value

    def _frame_range(self, frame: int):
        "Get the trajectory frame numbers over which we will average values"
        return frames_to_average(frame, self.n_frames, average=self.average)

    def _cache_ordered(self) -> np.ndarray:
        "Return the cached frames as a 3D array, in chronological order"
        keys = list(self.cache.keys())
        keys.sort()
        return np.array([self.cache[k] for k in keys])

    def adjust_periodic_positions(
        self, pos1: np.ndarray, pos2: np.ndarray
    ) -> np.ndarray:
        "Returns the input pos2 with a periodic correction potentially applied"
        if self.correct_periodic and self.is_orthorhombic:
            return correct_periodic_positions(pos1, pos2, self.universe.dimensions[:3])
        else:
            return pos2

    def position_cache_mean(self, frame: int) -> np.ndarray:
        "Return the mean position from the currently cached positions"
        self.update_position_cache(frame)

        if self.average == 0:
            return self.cache[frame]

        array = self._cache_ordered()
        if self.correct_periodic and self.is_orthorhombic:
            # we want to correct the periodic boundary crossing in refernce to the fist
            # frame we are averaging
            for i, pos in enumerate(array):
                if i == 0:
                    continue
                array[i] = self.adjust_periodic_positions(array[0], pos)

        return np.mean(array, axis=0)

    def set_frame(self, frame: int) -> None:
        """
        Update the positions, selections and calculations for this trajectory, based on
        frame number of the current scene, not the frame number of the Universe
        """
        self._update_positions(frame)
        self._update_selections()
        self._update_calculations()

    def _position_at_frame(self, frame: int) -> np.ndarray:
        "Return the atom positions at the given universe frame number"
        self.uframe = frame
        return self.univ_positions

    def update_position_cache(self, frame: int, cache_ahead: bool = True) -> None:
        "Update the currently cached positions, based on the new frame"
        # get the individual frame numbers that we will be caching
        frames_to_cache = self._frame_range(frame)

        # if we should be looking ahead by 1 for interpolating, ensure we are caching 1
        # frame ahead so when the frame changes we already have it stored and aren't
        # double dipping
        if (
            len(frames_to_cache) == 1
            and frames_to_cache[0] != (self.n_frames - 1)
            and cache_ahead
        ):
            frames_to_cache = np.array(
                (frames_to_cache[0], frames_to_cache[0] + 1), dtype=int
            )

        # only cleanup the cache if we have more than 2 frame stored, helps when moving
        # forward or back a single frame
        if len(self.cache) > 2:
            # remove any frames that no longer need to be cached
            to_remove = [f for f in self.cache if f not in frames_to_cache]
            for f in to_remove:
                del self.cache[f]

        # update the cache with any frames that are not yet cached
        for f in frames_to_cache:
            if f not in self.cache:
                self.cache[f] = self._position_at_frame(f)

    def frame_mapper(self, frame: int):
        return frame_mapper(
            frame=frame,
            subframes=self.subframes,
            offset=self.offset,
            mapping=self.frame_mapping,
        )

    def _update_positions(self, frame):
        """
        The function that will be called when the frame changes.
        It will update the positions and selections of the atoms in the scene.
        """
        if not self.update_with_scene:
            # get the current positions for the original frame and set
            # those on the object
            self.position = self._position_at_frame(frame)
            return

        # the rest of the code below is when update_with_scene is enabled
        # get the two frames of the trajectory to potentially access data from
        # uframe_current, uframe_next = [self.frame_mapper(x) for x in (frame, frame + 1)]
        uframe_current = self.frame_mapper(frame)
        uframe_next = uframe_current + 1
        last_frame = self.n_frames - 1
        if uframe_current >= last_frame:
            uframe_current = last_frame
            uframe_next = uframe_current
        # update the frame_hidden property for the UI
        # TODO: changing a bpy.prop here could lead to a
        #   `AttributeError: Writing to ID classes in this context is not allowed:` error
        self._frame = uframe_current

        if self.subframes > 0 and self.interpolate:
            # if we are adding subframes and interpolating, then we get the positions
            # at the two universe frames, then interpolate between them, potentially
            # correcting for any periodic boundary crossing
            pos_current = self.position_cache_mean(
                uframe_current,
            )
            pos_next = self.position_cache_mean(uframe_next)

            # if we are averaging, then we have already applied periodic correction
            # and we can skip this step
            if self.correct_periodic and self.is_orthorhombic and self.average == 0:
                pos_next = correct_periodic_positions(
                    pos_current,
                    pos_next,
                    dimensions=self.universe.dimensions[:3] * self.world_scale,
                )

            # interpolate between the two sets of positions
            self.position = db.lerp(
                pos_current, pos_next, t=fraction(frame, self.subframes + 1)
            )
        elif self.average > 0:
            # if we have subframes then we get the potential mean positions for the cached
            # frames that we are looking at
            self.position = self.position_cache_mean(uframe_current)
        else:
            # otherwise just get the current positions for the relevant frame and set
            # those on the object
            self.position = self._position_at_frame(uframe_current)

    def __repr__(self):
        return f"<Trajectory, `universe`: {self.universe}, `object`: {self.object}"

    def add_style(
        self,
        style: StyleBase | str = "spheres",
        color: str | None = "common",
        selection: str | AtomGroup | None = None,
        material: bpy.types.Material | str | None = None,
        name: str | None = None,
    ):
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

    def _get_3d_bbox(self, selection: mda.AtomGroup) -> list[tuple]:
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

    def get_view(self, selection: str | AtomGroup = None, frame: int = None) -> None:
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
