from typing import Dict, Callable

import bpy
import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from math import floor, remainder

from enum import Enum

from ... import data
from ..entity import MolecularEntity, EntityType
from ...blender import coll, nodes, path_resolve
import databpy
from ...utils import (
    correct_periodic_positions,
    frame_mapper,
    frames_to_average,
    fraction,
)
from .selections import Selection


class Trajectory(MolecularEntity):
    def __init__(self, universe: mda.Universe, world_scale: float = 0.01):
        super().__init__()
        self.universe: mda.Universe = universe
        self.selections: Dict[str, Selection] = {}
        self.calculations: Dict[str, Callable] = {}
        self.world_scale = world_scale
        self.frame_mapping: npt.NDArray[np.int64] | None = None
        self.cache: dict = {}
        self._entity_type = EntityType.MD

    def selection_from_ui(self, ui_item) -> Selection:
        self.selections[ui_item.name] = Selection(
            universe=self.universe,
            selection_str=ui_item.selection_str,
            name=ui_item.name,
            updating=ui_item.updating,
            periodic=ui_item.periodic,
        )
        self.apply_selection(self.selections[ui_item.name])

        return self.selections[ui_item.name]

    def add_selection(
        self,
        selection_str: str,
        name: str,
        updating: bool = True,
        periodic: bool = True,
    ):
        "Adds a new selection with the given name, selection string and selection parameters."
        obj = self.object

        obj.mn_trajectory_selections.add()
        sel = obj.mn_trajectory_selections[-1]
        sel.name = name
        sel.selection_str = selection_str
        sel.updating = updating
        sel.periodic = periodic

        self.selection_from_ui(sel)

        return sel

    def add_selection_from_atomgroup(self, atomgroup: mda.AtomGroup, name: str = ""):
        "Create a Selection object from an AtomGroup"
        selection = Selection.from_atomgroup(atomgroup, name=name)

        obj = self.object
        obj.mn_trajectory_selections.add()
        sel = obj.mn_trajectory_selections[-1]

        if not atomgroup.__class__.__name__ == "UpdatingAtomGroup":
            sel.immutable = True
        sel.name = selection.name
        sel.selection_str = selection.selection_str
        sel.updating = selection.updating
        sel.periodic = selection.periodic

        self.selections[selection.name] = selection
        self.apply_selection(selection)
        return sel

    def apply_selection(self, selection: Selection):
        "Set the boolean attribute for this selection on the mesh of the object"
        self.set_boolean(selection.to_mask(), name=selection.name)

    @property
    def is_orthorhombic(self):
        dim = self.universe.dimensions
        if dim is None:
            return False

        return np.allclose(dim[3:], 90.0)

    @property
    def atoms(self) -> mda.AtomGroup:
        return self.universe.atoms

    @property
    def n_atoms(self) -> int:
        return self.atoms.n_atoms

    @staticmethod
    def bool_selection(ag, selection, **kwargs) -> np.ndarray:
        return np.isin(ag.ix, ag.select_atoms(selection, **kwargs).ix).astype(bool)

    @property
    def univ_positions(self) -> np.ndarray:
        return self.atoms.positions * self.world_scale

    @property
    def bonds(self) -> np.ndarray:
        # the code to remap indices for a selection was removed as we don't subset the trajectory anymore
        # when importing it, everything is imported and the selections just update
        if hasattr(self.atoms, "bonds"):
            return self.atoms.bonds.indices
        else:
            return None

    @property
    def elements(self) -> np.ndarray:
        if hasattr(self.atoms, "elements"):
            return self.atoms.elements

        try:
            guessed_elements = [
                x
                if x in data.elements.keys()
                else mda.topology.guessers.guess_atom_element(x)
                for x in self.atoms.names
            ]
            return np.array(guessed_elements)

        except Exception:
            return np.repeat("X", self.n_atoms)

    @property
    def atomic_number(self) -> np.ndarray:
        return np.array(
            [
                data.elements.get(element, data.elements.get("X")).get("atomic_number")
                for element in self.elements
            ]
        )

    @property
    def vdw_radii(self) -> np.ndarray:
        return (
            np.array(
                [
                    data.elements.get(element, {}).get("vdw_radii", 100)
                    for element in self.elements
                ]
            )
            * 0.01  # pm to Angstrom
            * self.world_scale  # Angstrom to world scale
        )

    @property
    def mass(self) -> np.ndarray:
        # units: daltons
        if hasattr(self.atoms, "masses"):
            return np.array([x.mass for x in self.atoms])
        else:
            masses = [
                data.elements.get(element, {"standard_mass": 0}).get("standard_mass")
                for element in self.elements
            ]
            return np.array(masses)

    @property
    def n_frames(self) -> int:
        return self.universe.trajectory.n_frames

    @property
    def res_id(self) -> np.ndarray:
        return self.atoms.resids

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

    @property
    def res_name(self) -> np.ndarray:
        return np.array(list(map(lambda x: x[0:3], self.atoms.resnames)))

    @property
    def atom_id(self) -> np.ndarray:
        return self.universe.atoms.atom_id

    @property
    def res_num(self) -> np.ndarray:
        return np.array(
            [
                data.residues.get(res_name, data.residues.get("UNK")).get(
                    "res_name_num"
                )
                for res_name in self.res_name
            ]
        )

    @property
    def b_factor(self) -> np.ndarray:
        if hasattr(self.atoms, "tempfactors"):
            return self.atoms.tempfactors
        else:
            return np.zeros(self.n_atoms)

    @property
    def segindices(self) -> np.ndarray:
        if hasattr(self.atoms, "segindices"):
            return self.atoms.segindices

    @property
    def chain_id(self) -> np.ndarray:
        if hasattr(self.atoms, "chainIDs"):
            return self.atoms.chainIDs
        else:
            return np.zeros(self.n_atoms)

    @property
    def chain_ids(self) -> np.ndarray:
        return np.unique(self.chain_id)

    @property
    def chain_id_num(self) -> np.ndarray:
        chain_ids, chain_id_index = np.unique(self.chain_id, return_inverse=True)
        return chain_id_index

    @property
    def atom_type(self) -> np.ndarray:
        return self.atoms.types

    @property
    def atom_type_unique(self) -> np.ndarray:
        return np.unique(self.atom_type)

    @property
    def atom_type_num(self) -> np.ndarray:
        try:
            atom_type_unique, atom_type_index = np.unique(
                self.atom_type, return_inverse=True
            )
            return atom_type_index
        except AttributeError:
            return None

    @property
    def atom_name(self) -> np.ndarray:
        if hasattr(self.atoms, "names"):
            return self.atoms.names
        else:
            return np.zeros(self.n_atoms)

    @property
    def atom_name_num(self) -> np.ndarray:
        if hasattr(self.atoms, "names"):
            return np.array(
                list(map(lambda x: data.atom_names.get(x, -1), self.atom_name))
            )
        else:
            return np.repeat(-1, self.n_atoms)

    @property
    def is_nucleic(self) -> np.ndarray:
        return self.bool_selection(self.atoms, "nucleic")

    @property
    def is_peptide(self) -> np.ndarray:
        return self.bool_selection(self.atoms, "protein or (name BB SC*)")

    @property
    def is_lipid(self) -> np.ndarray:
        return np.isin(self.atoms.resnames, data.lipid_names)

    @property
    def is_backbone(self) -> np.ndarray:
        return self.bool_selection(self.atoms, "backbone or nucleicbackbone or name BB")

    @property
    def is_alpha_carbon(self) -> np.ndarray:
        return self.bool_selection(self.atoms, "name CA or name BB")

    @property
    def is_solvent(self) -> np.ndarray:
        return self.bool_selection(
            self.atoms, "name OW or name HW1 or name HW2 or resname W or resname PW"
        )

    @property
    def _attributes_2_blender(self):
        """
        The attributes that will be added to the Blender object.
        """
        return {
            "atomic_number": {
                "value": self.atomic_number,
                "type": "INT",
                "domain": "POINT",
            },
            "vdw_radii": {
                "value": self.vdw_radii,
                "type": "FLOAT",
                "domain": "POINT",
            },
            "mass": {
                "value": self.mass,
                "type": "FLOAT",
                "domain": "POINT",
            },
            "res_id": {
                "value": self.res_id,
                "type": "INT",
                "domain": "POINT",
            },
            "segid": {
                "value": self.segindices,
                "type": "INT",
                "domain": "POINT",
            },
            "res_name": {
                "value": self.res_num,
                "type": "INT",
                "domain": "POINT",
            },
            "b_factor": {
                "value": self.b_factor,
                "type": "FLOAT",
                "domain": "POINT",
            },
            "chain_id": {
                "value": self.chain_id_num,
                "type": "INT",
                "domain": "POINT",
            },
            "atom_types": {
                "value": self.atom_type_num,
                "type": "INT",
                "domain": "POINT",
            },
            "atom_name": {
                "value": self.atom_name_num,
                "type": "INT",
                "domain": "POINT",
            },
            "is_backbone": {
                "value": self.is_backbone,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
            "is_alpha_carbon": {
                "value": self.is_alpha_carbon,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
            "is_solvent": {
                "value": self.is_solvent,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
            "is_nucleic": {
                "value": self.is_nucleic,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
            "is_lipid": {
                "value": self.is_lipid,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
            "is_peptide": {
                "value": self.is_peptide,
                "type": "BOOLEAN",
                "domain": "POINT",
            },
        }

    def save_filepaths_on_object(self) -> None:
        obj = self.object
        obj.mn.filepath_topology = str(path_resolve(self.universe.filename))
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

    def _create_object(
        self, style: str = "vdw", name: str = "NewUniverseObject"
    ) -> bpy.types.Object:
        self.object = databpy.create_object(
            name=name,
            collection=coll.mn(),
            vertices=self.univ_positions,
            edges=self.bonds,
        )

        for att_name, att in self._attributes_2_blender.items():
            try:
                self.store_named_attribute(
                    data=att["value"],
                    name=att_name,
                    atype=att["type"],
                    domain=att["domain"],
                )
            except Exception as e:
                print(e)

        if hasattr(self.atoms, "segindices"):
            segs = []
            for seg in self.atoms.segments:
                segs.append(seg.atoms[0].segid)

            self.object["segments"] = segs

        if style is not None:
            nodes.create_starting_node_tree(
                self.object, style=style, name=f"MN_{self.name}"
            )

    def create_object(
        self,
        name: str = "NewUniverseObject",
        style: str = "vdw",
    ):
        self._create_object(style=style, name=name)

        self.object["chain_ids"] = self.chain_ids

        if hasattr(self, "atom_type_unique"):
            self.object["atom_type_unique"] = self.atom_type_unique

        self.object.mn.entity_type = self._entity_type.value
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
        objs_to_update = [obj for obj in bpy.data.objects if obj.uuid == self.uuid]

        # mark all selections for cleanup if they are no longer relevant
        for selection in self.selections.values():
            selection.cleanup = True

        for obj in objs_to_update:
            for sel in obj.mn_trajectory_selections:
                # try and get a corresponding selection for this named selection
                # if the selection can't be found we create one
                selection = self.selections.get(sel.name)
                if selection is None:
                    selection = self.selection_from_ui(sel)
                elif sel.updating:
                    # if the selection string has, or some of the parameters about it have
                    # changed then we have to change the selection's AtomGroup before we
                    # apply the selection to the mesh
                    if (
                        selection.selection_str != sel.selection_str
                        or selection.periodic != sel.periodic
                    ):
                        selection.change_selection(
                            selection_str=sel.selection_str,
                            name=sel.name,
                            updating=sel.updating,
                            periodic=sel.periodic,
                        )

                    # Apply the selection to the actual mesh in the form of a boolean
                    # named attribute
                    self.apply_selection(selection)

                # mark the selection to not be cleaned up, and add any message from the
                # selection to the UI selection item for display in the UI
                selection.cleanup = False
                sel.message = selection.message

        # remove all of the attributes and selections that are marked for cleanup because
        # they are no longer being used by any objects for selections
        for name in list(self.selections):
            if self.selections[name].cleanup:
                try:
                    self.object.data.attributes.remove(
                        self.object.data.attributes[name]
                    )
                    self.selections.pop(name)
                except Exception as e:
                    print(e)

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
    def update_with_scene(self) -> bool:
        return self.object.mn.update_with_scene

    @update_with_scene.setter
    def update_with_scene(self, value: bool) -> None:
        self.object.mn.update_with_scene = value

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
        return frames_to_average(frame, self.average)

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
        self._frame = self.frame_mapper(frame)
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
        if len(frames_to_cache) == 1 and cache_ahead:
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
        # get the two frames of the trajectory to potentially access data from
        # uframe_current, uframe_next = [self.frame_mapper(x) for x in (frame, frame + 1)]
        uframe_current = self.frame_mapper(frame)
        uframe_next = uframe_current + 1

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
            self.position = databpy.lerp(
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
