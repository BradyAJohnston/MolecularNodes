import bpy
from bpy.app.handlers import persistent

try:
    import MDAnalysis as mda
except ImportError:
    HAS_mda = False
    import types

    class MockAtomGroup:
        pass

    class MockUniverse:
        pass

    mda = types.ModuleType("MDAnalysis")
    mda.Universe = MockUniverse
    mda.AtomGroup = MockAtomGroup
    mda.core = types.ModuleType("core")
    mda.topology = types.ModuleType("topology")

else:
    HAS_mda = True
import numpy as np
import warnings
import pickle
from typing import Union, List, Dict

from ... import data
from ...pkg import start_logging
from ...blender import (
    coll, obj, nodes
)
from ...utils import lerp
from .molecule import Molecule, AtomAttribute

class MDAnalysisSession:
    """mock class in order that the module not breaks"""
    pass

class MDA(Molecule): # consistent namingy
    """
    "sibling" class MoleculeAtomArray
    """

    def __init__(self, atoms: mda.AtomGroup, name="MDAnalysisSession", world_scale: float = 0.01, style:str= "vdw", in_memory: bool = False, **kwargs):

        log = start_logging(logfile_name=name)
        if not HAS_mda:
            raise ImportError("MDAnalysis is not installed.")
            exit(1)

        super().__init__(name=name, atoms=atoms, **kwargs)

        self.world_scale = world_scale
        self.style = style
        self.in_memory = in_memory

    def __len__(self) -> int:
        return len(self._atoms)
    
    @staticmethod
    def bool_selection(_atoms, selection) -> np.ndarray:
        return np.isin(_atoms.ix, _atoms.select_atoms(selection).ix).astype(bool)

    @property
    def positions(self) -> np.ndarray:
        return self._atoms.positions * self.world_scale

    @property
    def bonds(self) -> List[List[int]]:
        if hasattr(self._atoms, "bonds"):
            bond_indices = self._atoms.bonds.indices
            atm_indices = self._atoms.indices
            bond_filtering = np.all(np.isin(bond_indices, atm_indices), axis=1)
            bond_indices = bond_indices[bond_filtering]

            index_map = {index: i for i, index in enumerate(self._atoms.indices)}

            bonds = [[index_map[bond[0]], index_map[bond[1]]]
                     for bond in bond_indices]
        else:
            bonds = []
        return bonds

    @property
    def elements(self) -> List[str]:
        try:
            elements = self._atoms.elements.tolist()
        except:
            try:
                elements = [
                    "BB" if x == "BB" else
                    "SC" if x.startswith("SC") else
                    "GL" if x.startswith("GL") else
                    "CD" if x.startswith("D") else
                    mda.topology.guessers.guess_atom_element(x) for x in self._atoms.atoms.names
                ]

            except:
                elements = ['X'] * self._atoms.n_atoms
        return elements

    @property
    def atomic_number(self) -> np.ndarray:
        return np.array(
            [data.elements.get(element,
                               data.elements.get('X'))
             .get('atomic_number') for element in self.elements]
        )

    @property
    def vdw_radii(self) -> np.ndarray:
        # pm to Angstrom
        return np.array(
            [data.elements.get(element,
                               data.elements.get('X'))
             .get('vdw_radii') for element in self.elements]) * 0.01 * self.world_scale

    @property
    def res_id(self) -> np.ndarray:
        return self.ag.resnums

    @property
    def res_name(self) -> np.ndarray:
        return np.array(list(map(lambda x: x[0:3], self._atoms.resnames)))

    @property
    def res_num(self) -> np.ndarray:
        return np.array(
            [data.residues.get(res_name,
                               data.residues.get('UNK'))
             .get('res_name_num') for res_name in self.res_name]
        )

    @property
    def b_factor(self) -> np.ndarray:
        if hasattr(self._atoms, "tempfactors"):
            return self._atoms.tempfactors
        else:
            return np.zeros(self._atoms.n_atoms)

    @property
    def chain_id(self) -> np.ndarray:
        if hasattr(self._atoms, "chainIDs"):
            return self._atoms.chainIDs
        else:
            return np.zeros(self._atoms.n_atoms)

    @property
    def chain_ids(self) -> np.ndarray:
        return np.unique(self.chain_id)

    @property
    def chain_id_num(self) -> np.ndarray:
        chain_ids, chain_id_index = np.unique(
            self.chain_id, return_inverse=True)
        return chain_id_index

    @property
    def atom_type(self) -> np.ndarray:
        return self._atoms.types

    @property
    def atom_type_unique(self) -> np.ndarray:
        return np.unique(self.atom_type)

    @property
    def atom_type_num(self) -> np.ndarray:
        atom_type_unique, atom_type_index = np.unique(
            self.atom_type, return_inverse=True)
        return atom_type_index

    @property
    def atom_name(self) -> np.ndarray:
        if hasattr(self._atoms, "names"):
            return self._atoms.names
        else:
            return np.zeros(self._atoms.n_atoms)

    @property
    def atom_name_num(self) -> np.ndarray:
        if hasattr(self._atoms, "names"):
            return np.array(list(map(lambda x: data.atom_names.get(x, -1), self.atom_name)))
        else:
            return np.repeat(-1, self._atoms.n_atoms)

    @property
    def is_nucleic(self) -> np.ndarray:
        return self.bool_selection(self._atoms, "nucleic")

    @property
    def is_peptide(self) -> np.ndarray:
        return self.bool_selection(self._atoms, "protein or (name BB SC*)")

    @property
    def is_lipid(self) -> np.ndarray:
        return np.isin(self._atoms.resnames, data.lipid_names)

    @property
    def is_backbone(self) -> np.ndarray:
        return self.bool_selection(self._atoms, "backbone or nucleicbackbone or name BB")

    @property
    def is_alpha_carbon(self) -> np.ndarray:
        return self.bool_selection(self._atoms, "name CA or name BB")

    @property
    def is_solvent(self) -> np.ndarray:
        return self.bool_selection(self._atoms, "name OW or name HW1 or name HW2 or resname W or resname PW")

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
                "domain": "POINT"
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


@persistent
def _rejuvenate_universe(scene):
    """
    Rejuvenate the session when the old Blend file is loaded.
    It will search through all the objects in the scene and
    find the ones that are molecules.
    It requires the pkl file to be in the default location and
    still exist.

    Warning:
    -------
    When a Blend file saved from another computer is loaded,
    the session will likely be lost.
    The Blend file also need to be opened from the same place
    (working directory) as when it is saved.
    """
    mol_objects = {}
    for object in bpy.data.objects:
        try:
            obj_type = object["type"]
            if obj_type == "molecule":
                mol_objects[object.name] = object
        except KeyError:
            pass

    if len(mol_objects) > 0:
        bpy.types.Scene.mda_session = MDAnalysisSession._rejuvenate(
            mol_objects)


@persistent
def _sync_universe(scene):
    """
    Sync the universe when the Blend file is saved.
    It will be saved as a .mda_session file in the
    same place as the Blend file).
    """
    if hasattr(bpy.types.Scene, "mda_session"):
        blender_save_loc = bpy.data.filepath
        if bpy.types.Scene.mda_session is not None:
            bpy.types.Scene.mda_session._dump(blender_save_loc)
