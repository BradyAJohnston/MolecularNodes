# Node-group asset "Import Volume" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputBoolean, InputInteger, InputString


class ImportVolume(AssetGeometryGroup):
    """
    Import Volume

    Parameters
    ----------
    path : InputString
        Path of the .vdb file to import. `{index}` is replaced by the index, to pick one of several numbered files
    index : InputInteger
        Index of the file to import, when Follow Frame is off
    follow_frame : InputBoolean
        Take the index from the scene frame, counted from Frame Start, to play the files as an animation
    frame_start : InputInteger
        Scene frame that shows the first file, when Follow Frame is on
    name : InputString
        Name to store the `density` grid under in the bundle of the output, as `MN Grid <Name>`. Color Sample Grid finds it by this name after a style has turned the volume into a mesh
    count : InputInteger
        Number of numbered files. The index is clamped to them, so frames past the last file hold it

    Inputs
    ------
    i.path : StringSocket
        Path of the .vdb file to import. `{index}` is replaced by the index, to pick one of several numbered files
    i.index : IntegerSocket
        Index of the file to import, when Follow Frame is off
    i.follow_frame : BooleanSocket
        Take the index from the scene frame, counted from Frame Start, to play the files as an animation
    i.frame_start : IntegerSocket
        Scene frame that shows the first file, when Follow Frame is on
    i.name : StringSocket
        Name to store the `density` grid under in the bundle of the output, as `MN Grid <Name>`. Color Sample Grid finds it by this name after a style has turned the volume into a mesh
    i.count : IntegerSocket
        Number of numbered files. The index is clamped to them, so frames past the last file hold it

    Outputs
    -------
    o.volume : GeometrySocket
        The imported volume
    """

    _name = "Import Volume"
    _asset_name = "Import Volume"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        path: StringSocket
        """Path of the .vdb file to import. `{index}` is replaced by the index, to pick one of several numbered files"""
        index: IntegerSocket
        """Index of the file to import, when Follow Frame is off"""
        follow_frame: BooleanSocket
        """Take the index from the scene frame, counted from Frame Start, to play the files as an animation"""
        frame_start: IntegerSocket
        """Scene frame that shows the first file, when Follow Frame is on"""
        name: StringSocket
        """Name to store the `density` grid under in the bundle of the output, as `MN Grid <Name>`. Color Sample Grid finds it by this name after a style has turned the volume into a mesh"""
        count: IntegerSocket
        """Number of numbered files. The index is clamped to them, so frames past the last file hold it"""

    class _Outputs(SocketAccessor):
        volume: GeometrySocket
        """The imported volume"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        path: InputString = "",
        index: InputInteger = 0,
        follow_frame: InputBoolean = False,
        frame_start: InputInteger = 1,
        name: InputString = "density",
        count: InputInteger = 1,
    ):
        super().__init__(
            **{
                "Path": path,
                "Index": index,
                "Follow Frame": follow_frame,
                "Frame Start": frame_start,
                "Name": name,
                "Count": count,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        path = tree.inputs.string(
            "Path",
            "",
            description="Path of the .vdb file to import. `{index}` is replaced by the index, to pick one of several numbered files",
        )
        index = tree.inputs.integer(
            "Index",
            0,
            description="Index of the file to import, when Follow Frame is off",
            min_value=0,
        )
        follow_frame = tree.inputs.boolean(
            "Follow Frame",
            False,
            description="Take the index from the scene frame, counted from Frame Start, to play the files as an animation",
        )
        frame_start = tree.inputs.integer(
            "Frame Start",
            1,
            description="Scene frame that shows the first file, when Follow Frame is on",
        )
        name = tree.inputs.string(
            "Name",
            "density",
            description="Name to store the `density` grid under in the bundle of the output, as `MN Grid <Name>`. Color Sample Grid finds it by this name after a style has turned the volume into a mesh",
        )
        count = tree.inputs.integer(
            "Count",
            1,
            description="Number of numbered files. The index is clamped to them, so frames past the last file hold it",
            min_value=1,
        )
        volume = tree.outputs.geometry("Volume", description="The imported volume")

        switch = follow_frame.switch.integer(
            index, g.FloatToInteger(g.SceneTime().o.frame).o.integer - frame_start
        )
        import_vdb = g.ImportVDB(path.format({"index": switch.max(0).min(count - 1)}))
        with g.Frame("Keep the grid in the bundle, where density styles pass it on"):
            store_bundle_item = g.StoreBundleItem.float(
                import_vdb.o.volume.bundle(),
                g.FormatString("MN Grid {name}", items={"name": name}),
                g.GetNamedGrid.float(import_vdb, "density").o.grid,
                structure_type="GRID",
            )
        import_vdb.o.volume.set_bundle(store_bundle_item) >> volume


ASSET = ImportVolume

ASSET_METADATA = {
    "description": "Import a volume from a .vdb file, or one of several numbered files picked by index or by the scene frame",
    "catalog_id": "35b7cc56-45fd-4113-8e88-3c7387edb2d8",
}
