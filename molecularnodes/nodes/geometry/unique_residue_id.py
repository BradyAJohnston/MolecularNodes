# Node-group asset 'Unique Residue ID' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from .atom_name import AtomName
from .integer_run import IntegerRun
from .residue_id import ResidueID


class UniqueResidueID(AssetGeometryGroup):
    """
    Unique Residue ID

    Outputs
    -------
    o.group_id : IntegerSocket
        A unique Group ID for eash residue
    """

    _name = "Unique Residue ID"
    _asset_name = "Unique Residue ID"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.unique_residue_id"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        group_id: IntegerSocket
        """A unique Group ID for eash residue"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        group_id = tree.outputs.integer(
            "Group ID", description="A unique Group ID for eash residue"
        )

        group = AtomName()
        boolean_math = (
            g.Compare.integer.equal(group, 1).o.result
            | g.Compare.integer.equal(group, 50)
            | IntegerRun(value=ResidueID()).o.is_different
        )
        accumulate_field = g.AccumulateField.point.integer(boolean_math)
        with g.Frame():
            accumulate_field.o.leading - 1 >> group_id


ASSET = UniqueResidueID

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
