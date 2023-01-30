"""
A subpackage for reading rotation matrices and translation vectors
for biological assemblies from different file formats.

The central functions are `get_transformations_`
"""

from abc import ABCMeta, abstractmethod


class AssemblyParser(metaclass=ABCMeta):

    @abstractmethod
    def list_assemblies(self):
        """
        Return a ``list`` of ``str`` containing the available assembly
        IDs.
        """
    
    @abstractmethod
    def get_transformations(self, assembly_id):
        """
        Parse the necessary transformations for a given
        assembly ID.
        
        Return a ``list`` of transformations for a set of chains
        transformations:

        transformations on sets of chains for this assembly
        |          chain IDs affected by the transformation
        |          |        3x3 rotation matrix
        |          |        |        translation vector 
        |          |        |        |
        list[tuple[ndarray, ndarray, ndarray]]]
        """