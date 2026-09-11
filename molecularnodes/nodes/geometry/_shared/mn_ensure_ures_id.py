# Node group ".MN_ensure_ures_id" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import CustomGeometryGroup, GeometrySocket, SocketAccessor
from nodebpy.types import InputGeometry
from ..set_ures_id import SetUResID


class MN_ensure_ures_id(CustomGeometryGroup):
    """
    If the `ures_id` attribute doesn't exist, create it

    Parameters
    ----------
    input : InputGeometry
        Input

    Inputs
    ------
    i.input : GeometrySocket
        Input

    Outputs
    -------
    o.output : GeometrySocket
        Output
    """

    _name = ".MN_ensure_ures_id"
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "If the `ures_id` attribute doesn't exist, create it"
    }

    class _Inputs(SocketAccessor):
        input: GeometrySocket
        """Input"""

    class _Outputs(SocketAccessor):
        output: GeometrySocket
        """Output"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        input: InputGeometry = None,
    ):
        super().__init__(**{"Input": input})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        input = tree.inputs.geometry("Input")
        output = tree.outputs.geometry("Output")

        with g.Frame("Check if `ures_id` exists"):
            get_attribute_names = g.GetAttributeNames(
                geometry=input,
                filter_data_type=True,
                data_type="Integer",
                filter_domain=True,
            )
            filter_list = get_attribute_names.o.names.filter(
                g.Compare.string.equal(get_attribute_names, "ures_id")
            )
            list_length = filter_list.list_length()
        switch = g.Switch.geometry(list_length, SetUResID(geometry=input), input)

        switch >> output
