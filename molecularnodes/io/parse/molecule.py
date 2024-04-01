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
    n_atoms (int): Number of atoms in the molecule
    _atoms (biotite.structure.AtomArrayStack): Contains the atom information of the molecule.
    entity_id (np.ndarray): The entity ID of the atoms in the molecule.

    Descriptors:
    ----
    coord (np.ndarray): The coordinates of the atoms in the molecule.
    bonds (np.ndarray): The bonds between the atoms in the molecule.
    box (np.ndarray): The box of the molecule.

    every attribute in the AtomArray._annot is also a descriptor.

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
    bonds = AtomAttribute()
    box = AtomAttribute()

    def __init__(self, name: str, atoms: struc.AtomArray):
        self._name = name
        self._n_atoms = len(atoms)
        self._atoms = atoms

        # set every attribute in the AtomArray._annot as a descriptor
        for attr in self._atoms._annot.keys():
            setattr(self.__class__, attr, AtomAttribute())

        self.entity_id: Optional[np.ndarray] = None
    
    @property
    def name(self) -> str:
        return self._name
    
    @property
    def n_atoms(self) -> int:
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
        return self._atoms.array_length()

    def __getattr__(self, name: str):
        """
        Get the attribute from the AtomArray if the attribute is not found in the Molecule instance.

        Raises
        ------
        AttributeError
            If the attribute is not found in the AtomArray.
        """
        if name in self._atoms._annot.keys():
            return getattr(self._atoms, name)
        else:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")
    
    def __getitem__(self, index: Union[int, list, slice]) -> "Molecule":
        return self.__class__(name=self.name, atoms=self._atoms[:,index], **self.__dict__)

    def __add__(self, other):
        self.coord = self.coord + other

    def __sub__(self, other):
        self.coord = self.coord - other

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

