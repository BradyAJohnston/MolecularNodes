# Node group ".Sample Atomic Attributes" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    GeometrySocket,
    IntegerSocket,
    SocketAccessor,
)
from nodebpy.types import InputGeometry, InputInteger


class SampleAtomicAttributes(CustomGeometryGroup):
    """
    .Sample Atomic Attributes

    Parameters
    ----------
    atoms : InputGeometry
        Atoms
    sample_atoms : InputGeometry
        Sample Atoms
    index : InputInteger
        Index

    Inputs
    ------
    i.atoms : GeometrySocket
        Atoms
    i.sample_atoms : GeometrySocket
        Sample Atoms
    i.index : IntegerSocket
        Index

    Outputs
    -------
    o.atoms : GeometrySocket
        Atoms
    """

    _name = ".Sample Atomic Attributes"
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atoms"""
        sample_atoms: GeometrySocket
        """Sample Atoms"""
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        atoms: GeometrySocket
        """Atoms"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        sample_atoms: InputGeometry = None,
        index: InputInteger = 0,
    ):
        super().__init__(
            **{"Atoms": atoms, "Sample Atoms": sample_atoms, "Index": index}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry("Atoms")
        sample_atoms = tree.inputs.geometry("Sample Atoms")
        index = tree.inputs.integer("Index", 0, default_input="INDEX")
        atoms_1 = tree.outputs.geometry("Atoms")

        (
            atoms
            >> g.TransferAttributes(
                source=sample_atoms,
                source_point_id=index,
                attribute_names=g.FieldToList(items={"String": "bond_type"}),
                exclude_names=True,
            )
            >> atoms_1
        )
