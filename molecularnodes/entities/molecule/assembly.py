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

        Each transformation is a ``dict`` with the keys:

        - ``"chain_ids"``: ``list[str]`` of chain IDs affected by the transformation
        - ``"matrix"``: 4x4 rotation, translation & scale matrix as nested lists
        - ``"pdb_model_num"``: ``int`` index of the chain set the transformation
          belongs to
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
