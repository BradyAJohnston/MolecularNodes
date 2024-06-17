from typing import Dict, List, Union

import bpy
import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from bpy.app.handlers import persistent

from ... import data
from ...blender import coll, nodes, obj
from ...logger import start_logging
from ...utils import lerp


@persistent
def update_universes(scene):
    for universe in scene.MNSession.universes:
        universe._update_trajectory(scene.frame_current)
        universe._update_selections()


class MNUniverseSelection:
    def __init__(
        self, universe: mda.Universe, selection_str, name, updating=True, periodic=True
    ):
        self.selection_str: str = selection_str
        self.selection_previous: str = selection_str
        self.name: str = name
        self.periodic: bool = periodic
        self.updating: bool = updating
        self.previous_name: str = name
        self.ag = universe.select_atoms(
            self.selection_str, updating=self.updating, periodic=self.periodic
        )
        self.mask_array = self._ag_to_mask()

    def _ag_to_mask(self) -> npt.NDArray[np.bool_]:
        "Uses the selection atom group to provide a boolean mask for the universe atoms."
        return np.isin(self.ag.universe.atoms.ix, self.ag.ix).astype(bool)

    def change_selection(self, selection_str, updating=True, periodic=True):
        self.periodic = periodic
        self.updating = updating
        self.selection_str = selection_str
        self.ag = self.ag.universe.select_atoms(
            selection_str, updating=updating, periodic=periodic
        )

    def selection_mask(self) -> npt.NDArray[np.bool_]:
        "Returns the selection as a 1D numpy boolean mask. If updating=True, recomputes selection."
        if self.updating:
            self.mask_array = self._ag_to_mask()
        return self.mask_array


class MNUniverse:
    def __init__(self, universe: mda.Universe, world_scale=0.01):
        self.universe: mda.Universe = universe
        self.bob: bpy.types.Object | None
        self.selections: Dict[str, MNUniverseSelection] = {}
        self.world_scale = world_scale
        self.object: bpy.types.Object | None = None
        self.name: str | None
        self.frame_mapping: npt.NDArray[np.in64] | None = None

    def add_selection(
        self,
        selection_str: str,
        name: str,
        updating: bool = True,
        periodic: bool = True,
    ) -> None:
        "Adds a new selection with the given name, selection string and selection parameters."
        self.selections[name] = MNUniverseSelection(
            universe=self.universe,
            selection_str=selection_str,
            name=name,
            updating=updating,
            periodic=periodic,
        )
        self.set_selection(self.selections[name])

    def set_selection(self, selection: MNUniverseSelection):
        obj.set_attribute(
            bob=self.object,
            name=selection.name,
            data=selection.selection_mask(),
            type="BOOLEAN",
        )

    def get_selection(self, name: str) -> npt.NDArray[np.bool_]:
        "Get a given selection as a boolean mask for the universe."
        try:
            for selection in self.selections:
                if selection.name == name:
                    return selection.selection_mask()
        except Exception as e:
            print(f"No matching selection. Error: {e}")

    def _update_selections(self):
        ui_selections = self.object.mn_universe_selections

        for sel in ui_selections:
            if not sel.updating:
                continue
            try:
                selection = self.selections[sel.name]
                # if the same selection string, test to see if updating. If not updating
                # then move on - if updating then call the set selection method to update
                # the selection on the mesh
                if selection.selection_str != sel.selection:
                    # if the selection string has changed, we need to delete the old selection
                    # and thus the old atom group, and create a new atom selection and
                    # AtomGroup
                    selection.change_selection(
                        selection_str=sel.selection,
                        updating=sel.updating,
                        periodic=sel.periodic,
                    )
                self.set_selection(selection)

            except KeyError as e:
                print(e)
                # if the selection doesn't exist for the universe, create a new one
                self.add_selection(
                    name=sel.name,
                    selection_str=sel.selection,
                    updating=sel.updating,
                    periodic=sel.periodic,
                )

        # cleanup the remaining selections - if they no longer exist in the selections
        # list then they should be removed
        for name, selection in self.selections.items():
            if name not in [sel.name for sel in ui_selections]:
                try:
                    self.object.data.attributes.remove(name)
                    self.selections.pop(name)
                except Exception as e:
                    print(e)
                    pass

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
    def positions(self) -> np.ndarray:
        return self.atoms.positions * self.world_scale

    @property
    def bonds(self) -> List[List[int]]:
        if hasattr(self.atoms, "bonds"):
            bond_indices = self.atoms.bonds.indices
            atm_indices = self.atoms.indices
            bond_filtering = np.all(np.isin(bond_indices, atm_indices), axis=1)
            bond_indices = bond_indices[bond_filtering]
            index_map = {
                index: i for i, index in enumerate(self.universe.atoms.indices)
            }

            bonds = [[index_map[bond[0]], index_map[bond[1]]] for bond in bond_indices]
        else:
            bonds = []
        return bonds

    @property
    def elements(self) -> List[str]:
        try:
            elements = self.atoms.elements.tolist()
        except Exception:
            try:
                elements = [
                    x
                    if x in data.elements.keys()
                    else mda.topology.guessers.guess_atom_element(x)
                    for x in self.atoms.names
                ]

            except Exception:
                elements = ["X"] * self.atoms.n_atoms
        return elements

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
        # pm to Angstrom
        return (
            np.array(
                [
                    data.elements.get(element, {}).get("vdw_radii", 100)
                    for element in self.elements
                ]
            )
            * 0.01
            * self.world_scale
        )

    @property
    def mass(self) -> np.ndarray:
        # units: daltons
        try:
            masses = np.array([x.mass for x in self.atoms])
        except mda.exceptions.NoDataError:
            masses = np.array(
                [
                    data.elements.get(element, {"standard_mass": 0}).get(
                        "standard_mass"
                    )
                    for element in self.elements
                ]
            )
        return masses

    @property
    def res_id(self) -> np.ndarray:
        return self.atoms.resnums

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

    def create_model(
        self,
        style: str = "vdw",
        name: str = "NewUniverseObject",
        subframes: int = 0,
        # in_memory: bool = False,
    ):
        bob = obj.create_object(
            name=name, collection=coll.mn(), vertices=self.positions, edges=self.bonds
        )
        self.object = bob
        self.name = bob.name

        for att_name, att in self._attributes_2_blender.items():
            obj.set_attribute(bob, att_name, att["value"], att["type"], att["domain"])

        bob["chain_ids"] = self.chain_ids
        bob["atom_type_unique"] = self.atom_type_unique
        bob.mn["subframes"] = subframes
        bob.mn["molecule_type"] = "md"

        if style is not None:
            nodes.create_starting_node_tree(bob, style=style, name=f"MN_{bob.name}")

        bpy.context.view_layer.objects.active = bob
        return bob

    @persistent
    def _update_trajectory(self, frame):
        """
        The function that will be called when the frame changes.
        It will update the positions and selections of the atoms in the scene.
        """
        if self.object is None:
            self.object = bpy.data.objects[self.name]
        universe = self.universe
        frame_mapping = self.frame_mapping
        subframes = self.object.mn["subframes"]

        if frame < 0:
            return None

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

        bob = self.object

        # set the trajectory at frame_a
        universe.trajectory[frame_a]

        if subframes > 0:
            fraction = frame % (subframes + 1) / (subframes + 1)

            # get the positions for the next frame
            locations_a = self.positions

            if frame_b < universe.trajectory.n_frames:
                universe.trajectory[frame_b]
            locations_b = self.positions

            # interpolate between the two sets of positions
            locations = lerp(locations_a, locations_b, t=fraction)
        else:
            locations = self.positions

        # update the positions of the underlying vertices
        obj.set_attribute(bob, "position", locations)
