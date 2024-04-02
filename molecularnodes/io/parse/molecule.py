from abc import ABC, abstractmethod
from typing import Optional, Any, Union, List, Tuple, Dict, TypeVar, Generic, Sequence
import warnings
import time
import numpy as np
import bpy
import biotite.structure as struc

import MDAnalysis as mda

from ... import blender as bl
from ... import utils, data, color

class AtomAttribute:
    """
        Provides setter and getter for the atom attributes, in a fashion, that the attributes
        are stored in the data-container of the class, but can be accessed directly from the class,
        so that there is no dubplicated data and it simplifies the interaction with the user.

        Example usage:

        class MDA:
            
            postions = AtomGroupAttribute() # attribute of the AtomGroup
            bonds = AtomGroupAttribute()    # ...

            def __init__(self, ag: struc.AtomGroup):
                self._atoms = atoms

        """

    def __set_name__(self, owner_class, prop_name):
        self.prop_name = prop_name
        
    def __set__(self, instance, value):
        setattr(instance._atoms, self.prop_name, value)
        
    def __get__(self, instance, owner_class):
            return getattr(instance._atoms, self.prop_name)

T = TypeVar('T', struc.AtomArray, struc.AtomArrayStack, mda.Universe, mda.AtomGroup)

class Molecule(ABC):

    """
    Every Molecule has the atrributes name, n_atoms, _atoms and entity_id.

    Attributes:
    ----

    - `name (str):` Name of the molecule.
    - `n_atoms (int):` Number of atoms in the molecule
    - `_atoms (Sequence[T]):` Data container which contains the atom information of the molecule.
        The Attributes are dynamically set as AtomAttribute descriptors in the child class.
    - `entity_id (np.ndarray):` The entity ID of the atoms in the molecule.

    Methods:
    ----

    - `__init__(self, name: str, atoms: List[T]):` Initialize the molecule.
    - `__len__(self) -> int:` Return the number of atoms in the molecule. Has to be defined in the child class.
    - `__repr__(self) -> str:` Return the string representation of the molecule.
    - `__getattr__(self, name: str):` Get the attribute from the AtomArray if the attribute is not found in the Molecule instance.
    - `__getitem__(self, index: Union[int, list, slice]) -> "Molecule":` Return a new Molecule instance with the sliced atoms.

    """

    def __init__(self, name : str, atoms : List[T]):
        self._name : str = name
        self._atoms : List[T ]= atoms
        self._n_atoms : int | None = None
        self.entity_id: np.ndarray | None = None
    
    @property
    def name(self) -> str:
        return self._name
    
    @property
    def n_atoms(self) -> int:
        if not self._n_atoms:
            self._n_atoms = len(self)
        return self._n_atoms
    
    def __repr__(self) -> str:
        return f"<Molecule object: {self.name}>"
    
    def _getattr__(self, name: str):
        """
        Get the attribute from the AtomArray if the attribute is not found in the Molecule instance.

        Raises
        ------
        AttributeError
            If the attribute is not found in the AtomArray.
        """
        try:
            return getattr(self._atoms, name)
        except AttributeError:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")
    
    def __getitem__(self, index: Union[int, list, slice]) -> "Molecule":
        return self.__class__(name=self.name, atoms=self._atoms[:,index], **self.__dict__)
    
    #TODO: __setattr__

    @abstractmethod
    def __len__(self) -> int:
        ...

class MoleculeAtomArray(Molecule):
    """
    Represents a molecule with an AtomArrayStruct as data container.

    Descriptors:
    ----
    - `coord (np.ndarray):` The coordinates of the atoms in the molecule.
    - `bonds (np.ndarray):` The bonds between the atoms in the molecule.
    - `box (np.ndarray):` The box of the molecule.

    - `every attribute` in the AtomArray._annot is also a descriptor.

    Methods:
    ----

    - `_get_structure(self):` Get the structure of the molecule.
    - `_assemblies(self, as_array=False):` Get the biological assemblies of the molecule.
    - `assemblies(self, as_array=False):` Get the biological assemblies of the molecule.


    Example usage:

    ```python
    class PDB(MoleculeAtomArray):
        def __init__(self, atoms: struc.AtomArrayStack, fpath:str, **kwargs):

            super().__init__(name="from-local-pdb-file", atoms=atoms)
            self.fpath = fpath

        @classmethod
        def _get_structure(cls, fpath: str) -> "PDB":
        ...
    
        mol = PDB._get_structure("path/to/file.pdb")
    """

    coord = AtomAttribute() # attributes in the AtomArray._annot
    bonds = AtomAttribute() # ...
    box = AtomAttribute()   # ...
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # set every attribute in the AtomArray._annot as a descriptor
        # -> the data container has to have the name _atoms
        for attr in self._atoms._annot.keys():
            setattr(self.__class__, attr, AtomAttribute())

    @classmethod
    @abstractmethod
    def _get_structure(cls, fpath: str) -> "MoleculeAtomArray":
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
    
    def __len__(self) -> int:
        return self._atoms.array_length()