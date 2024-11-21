from typing import Dict, Callable

import bpy
import MDAnalysis as mda
import numpy as np
import numpy.typing as npt

from ... import data
from ..entity import MolecularEntity
from ...blender import coll, nodes, path_resolve
from ... import bpyd
from ...utils import correct_periodic_positions
from .selections import Selection, TrajectorySelectionItem


class Trajectory(MolecularEntity):
    def __init__(self, universe: mda.Universe, world_scale: float = 0.01):
        super().__init__()
        self.universe: mda.Universe = universe
        self.selections: Dict[str, Selection] = {}
        self.calculations: Dict[str, Callable] = {}
        self.world_scale = world_scale
        self.frame_mapping: npt.NDArray[np.in64] | None = None
        bpy.context.scene.MNSession.entities[self.uuid] = self

    def selection_from_ui(self, ui_item: TrajectorySelectionItem) -> Selection:
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
    ) -> TrajectorySelectionItem:
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
    def subframes(self):
        obj = self.object
        if obj is None:
            return None
        return obj.mn.subframes

    @subframes.setter
    def subframes(self, value: int):
        obj = self.object
        if obj is None:
            return None
        obj.mn.subframes = value

    @property
    def offset(self) -> int:
        try:
            return self.object.mn.offset
        except AttributeError:
            return None

    @offset.setter
    def offset(self, value: int):
        self.object.mn.offset = value

    @property
    def interpolate(self) -> bool:
        obj = self.object
        if obj is None:
            return None
        return obj.mn.interpolate

    @interpolate.setter
    def interpolate(self, value: bool):
        obj = self.object
        if obj is None:
            return None
        obj.mn.interpolate = value

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
    def res_id(self) -> np.ndarray:
        return self.atoms.resnums

    @property
    def frame(self) -> int:
        return self.universe.trajectory.frame

    @frame.setter
    def frame(self, value) -> None:
        if self.universe.trajectory.frame != value:
            self.universe.trajectory[value]

    @property
    def res_name(self) -> np.ndarray:
        return np.array(list(map(lambda x: x[0:3], self.atoms.resnames)))

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
        atom_type_unique, atom_type_index = np.unique(
            self.atom_type, return_inverse=True
        )
        return atom_type_index

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

    def create_object(
        self,
        style: str = "vdw",
        name: str = "NewUniverseObject",
        subframes: int = 0,
        # in_memory: bool = False,
    ):
        obj = bpyd.create_object(
            name=name,
            collection=coll.mn(),
            vertices=self.univ_positions,
            edges=self.bonds,
        )
        self.object = obj

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

            obj["segments"] = segs

        obj["chain_ids"] = self.chain_ids
        obj["atom_type_unique"] = self.atom_type_unique
        self.subframes = subframes
        obj.mn.molecule_type = "md"
        self.save_filepaths_on_object()

        if style is not None:
            nodes.create_starting_node_tree(obj, style=style, name=f"MN_{obj.name}")

        bpy.context.view_layer.objects.active = obj
        obj.mn.uuid = self.uuid

        return obj

    def _update_calculations(self):
        for name, func in self.calculations.items():
            try:
                self.store_named_attribute(data=func(self.universe), name=name)
            except Exception as e:
                print(e)

    def _update_selections(self):
        objs_to_update = [obj for obj in bpy.data.objects if obj.mn.uuid == self.uuid]

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
                else:
                    pass

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

    def _update_positions(self, frame):
        """
        The function that will be called when the frame changes.
        It will update the positions and selections of the atoms in the scene.
        """
        universe = self.universe
        frame_mapping = self.frame_mapping
        obj = self.object

        subframes: int = obj.mn.subframes
        interpolate: bool = obj.mn.interpolate
        offset: int = obj.mn.offset

        # we subtraect the offset, a negative offset value ensures that the trajectory starts
        # playback that many frames before 0 and a positive value ensures we start the
        # playback after 0
        frame -= offset
        # for actually getting frames from the trajectory we need to clamp it to a lower
        # bound of 0 which will be the start frame for the trajectory
        frame = max(frame, 0)

        if frame_mapping:
            # add the subframes to the frame mapping
            frame_map = np.repeat(frame_mapping, subframes + 1)
            # get the current and next frames
            frame_a = frame_map[frame]
            frame_b = frame_map[frame + 1]

        else:
            # get the initial frame
            if subframes == 0:
                frame_a = frame
            else:
                frame_a = int(frame / (subframes + 1))

            # get the next frame
            frame_b = frame_a + 1

        if frame_a >= universe.trajectory.n_frames:
            return None

        # set the trajectory at frame_a
        self.frame = frame_a

        if subframes > 0 and interpolate:
            fraction = frame % (subframes + 1) / (subframes + 1)

            # get the positions for the next frame
            positions_a = self.univ_positions

            if frame_b < universe.trajectory.n_frames:
                self.frame = frame_b
            positions_b = self.univ_positions

            if obj.mn.correct_periodic and self.is_orthorhombic:
                positions_b = correct_periodic_positions(
                    positions_a,
                    positions_b,
                    dimensions=universe.dimensions[:3] * self.world_scale,
                )

            # interpolate between the two sets of positions
            self.position = bpyd.lerp(positions_a, positions_b, t=fraction)
        else:
            self.position = self.univ_positions

    def __repr__(self):
        return f"<Trajectory, `universe`: {self.universe}, `object`: {self.object}"
