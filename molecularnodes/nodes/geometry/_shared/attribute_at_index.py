# Node group 'Attribute at Index' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    IntegerSocket,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputInteger, InputString


class AttributeAtIndex(CustomGeometryGroup):
    """
    Attribute at Index

    Parameters
    ----------
    index : InputInteger
        Index
    name : InputString
        Name

    Inputs
    ------
    i.index : IntegerSocket
        Index
    i.name : StringSocket
        Name

    Outputs
    -------
    o.value : IntegerSocket
        Value
    """

    _name = "Attribute at Index"
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""
        name: StringSocket
        """Name"""

    class _Outputs(SocketAccessor):
        value: IntegerSocket
        """Value"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        index: InputInteger = 0,
        name: InputString = "res_id",
    ):
        super().__init__(**{"Index": index, "Name": name})

    def _build_group(self, tree):
        index = tree.inputs.integer("Index", 0, min_value=0, default_input="INDEX")
        name = tree.inputs.string("Name", "res_id")
        value = tree.outputs.integer("Value")

        g.NamedAttribute.integer(name).o.attribute.point.at(index) >> value
