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
        |          |        4x4 rotation, translation & scale matrix
        |          |        |
        list[tuple[ndarray, ndarray]]]
        """

    @abstractmethod
    def get_assemblies(self):
        """
        Parse all the transformations for each assembly, returning a dictionary of
        key:value pairs of assembly_id:transformations. The transformations list
        comes from the `get_transformations(assembly_id)` method.

        Dictionary of all assemblies
        |     Assembly ID
        |     |   List of transformations to create biological assembly.
        |     |   |
        dict{'1', list[transformations]}

        """
