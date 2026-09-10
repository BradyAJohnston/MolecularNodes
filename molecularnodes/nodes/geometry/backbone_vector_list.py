# Node-group asset "Backbone Vector List" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputGeometry, InputInteger
from .backbone_positions import BackbonePositions
from .sample_position import SamplePosition


class BackboneVectorList(AssetGeometryGroup):
    """
    Backbone Vector List

    Parameters
    ----------
    ca_atoms : InputGeometry
        CA Atoms
    index : InputInteger
        The `Index` at which to sample the `Position` field from

    Inputs
    ------
    i.ca_atoms : GeometrySocket
        CA Atoms
    i.index : IntegerSocket
        The `Index` at which to sample the `Position` field from

    Outputs
    -------
    o.value : VectorSocket
        Output list with evaluated field values
    o.ca_atoms : GeometrySocket
        CA Atoms
    """

    _name = "Backbone Vector List"
    _asset_name = "Backbone Vector List"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")

    class _Inputs(SocketAccessor):
        ca_atoms: GeometrySocket
        """CA Atoms"""
        index: IntegerSocket
        """The `Index` at which to sample the `Position` field from"""

    class _Outputs(SocketAccessor):
        value: VectorSocket
        """Output list with evaluated field values"""
        ca_atoms: GeometrySocket
        """CA Atoms"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        ca_atoms: InputGeometry = None,
        index: InputInteger = 0,
    ):
        super().__init__(**{"CA Atoms": ca_atoms, "Index": index})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        ca_atoms = tree.inputs.geometry("CA Atoms")
        index = tree.inputs.integer(
            "Index",
            0,
            description="The `Index` at which to sample the `Position` field from",
            hide_value=True,
            default_input="INDEX",
        )
        value = tree.outputs.vector(
            "Value", description="Output list with evaluated field values"
        )
        ca_atoms_1 = tree.outputs.geometry("CA Atoms")

        group = BackbonePositions(method="Read")
        index_switch = g.IndexSwitch.vector(
            g.Index(),
            (
                SamplePosition(geometry=ca_atoms, position=group.o.c, index=index),
                SamplePosition(geometry=ca_atoms, position=group.o.ca, index=index),
                SamplePosition(geometry=ca_atoms, position=group.o.n, index=index),
            ),
        )
        g.FieldToList(count=3, items={"Value": index_switch}) >> value

        ca_atoms >> ca_atoms_1


ASSET = BackboneVectorList

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
