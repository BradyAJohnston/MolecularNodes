from abc import ABC, abstractmethod
from typing import Optional, Any, Union, List, Tuple, Dict
import warnings
import time
import numpy as np
import bpy
import biotite.structure as struc

from ... import blender as bl
from ... import utils, data, color

class AtomAttribute:
    """
    Provides setter and getter for the atom attributes [biotite.structure.AtomArray (np.ndarray)](https://github.com/biotite-dev/biotite/blob/master/src/biotite/structure/atoms.py). 
    The AtomArray is a np.ndarray with a structured dtype. 

    The setter:

    The getter:
    """

    def __set_name__(self, owner_class, prop_name):
        self.prop_name = prop_name
        
    def __set__(self, instance, value):
        setattr(instance._atoms, self.prop_name, value)
        
    def __get__(self, instance, owner_class):
            return getattr(instance._atoms, self.prop_name)

   

class Molecule(ABC):
    """
    Represents a molecule.

    Attributes:
    ----
    name (str): Name of the molecule.
    coord (AtomField): Coordinates of the atoms.
    chain_id (AtomField): Chain IDs of the atoms.
    res_name (AtomField): Residue names of the atoms.
    entity_id (AtomField): Entity IDs of the atoms.
    hetero (AtomField): Hetero flags of the atoms.
    atom_name (AtomField): Atom names of the atoms.
    n_atoms (int): Number of atoms in the molecule.

    Methods:
    ----
    __init__(self, name: str, atoms: struc.AtomArray): Initialize the molecule.
    _get_structure(self): Get the structure of the molecule.
    assemblies(self, as_array=False): Get the biological assemblies of the molecule.
    __repr__(self): Return the string representation of the molecule.
    __len__(self): Return the number of atoms in the molecule.

    Arithmetic Operators - delegate to numpy


    """

    coord = AtomAttribute()
    chain_id = AtomAttribute()
    res_name = AtomAttribute()
    hetero = AtomAttribute()
    atom_name = AtomAttribute()


    # TODO write custom method to apply filters like vstack and filter, sort, apply, ... to all attributes
    def __init__(self, name: str, atoms: struc.AtomArray):
        self._name = name
        self._n_atoms = len(atoms)
        self._atoms = atoms
        self.entity_id: Optional[np.ndarray] = None
    
    @property
    def name(self) -> str:
        return self._name
    
    @property
    def n_atoms(self) -> int:
        if not self._n_atoms:
            self._n_atoms = len(self)
        return self._n_atoms

    @classmethod
    @abstractmethod
    def _get_structure(cls, fpath: str):
        """
        Get the structural information of the molecule with biotite.structure.io.

        Returns
        -------

        molecule : Molecule

        """
        ...

    @abstractmethod
    def _assemblies(self, as_array=False):
        """
        Get the biological assemblies of the molecule.

        Parameters
        ----------
        as_array : bool, optional
            Whether to return the assemblies as an array of quaternions.
            Default is False.

        Returns
        -------
        dict or None
            The biological assemblies of the molecule, as a dictionary of
            transformation matrices, or None if no assemblies are available.
        """
        
        ...

    def assemblies(self, as_array=False):
        """
        Get the biological assemblies of the molecule.

        Parameters
        ----------
        as_array : bool, optional
            Whether to return the assemblies as an array of quaternions.
            Default is False.

        Returns
        -------
        dict or None
            The biological assemblies of the molecule, as a dictionary of
            transformation matrices, or None if no assemblies are available.
        """
        from biotite import InvalidFileError
        try:
            assemblies_info = self._assemblies()
        except InvalidFileError:
            return None

        if isinstance(assemblies_info, dict) and as_array:
            return utils.array_quaternions_from_dict(assemblies_info)

        return assemblies_info

    def __repr__(self) -> str:
        return f"<Molecule object: {self.name}>"
    
    def __len__(self) -> int:
        return self._atoms.shape[1]
    
    def __getitem__(self, index: Union[int, list, slice]) -> "Molecule":
        return self.__class__(name=self.name, atoms=self._atoms[:,index])
        
    
    def __eq__(self, __value: object):
        pass

    def __add__(self, other):
        self.coord = self.coord + other

    def __sub__(self, other):
        self.coor = self.coord - other

    def __mul__(self, other):
        pass

    def __truediv__(self, other):
        pass

    def __iadd__(self, other):
        pass

    def __isub__(self, other):
        pass

    def __imul__(self, other):
        pass    

