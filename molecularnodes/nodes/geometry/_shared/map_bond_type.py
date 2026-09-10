# Node group '.map_bond_type' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import CustomGeometryGroup, IntegerSocket, SocketAccessor


class Map_bond_type(CustomGeometryGroup):
    """
    .map_bond_type

    Outputs
    -------
    o.bond_count : IntegerSocket
        Bond Count
    """

    _name = ".map_bond_type"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        bond_count: IntegerSocket
        """Bond Count"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        bond_count = tree.outputs.integer("Bond Count")

        (
            g.IndexSwitch.integer(
                g.NamedAttribute.integer("bond_type").o.attribute, (1, 1, 2, 3, 1, 1, 2)
            )
            >> bond_count
        )
