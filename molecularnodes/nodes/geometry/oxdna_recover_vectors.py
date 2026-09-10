# Node-group asset "oxDNA Recover Vectors" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputGeometry


class OxDNARecoverVectors(AssetGeometryGroup):
    """
    oxDNA Recover Vectors

    Parameters
    ----------
    geometry : InputGeometry
        Geometry

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry

    Outputs
    -------
    o.atoms : GeometrySocket
        Atoms
    """

    _name = "oxDNA Recover Vectors"
    _asset_name = "oxDNA Recover Vectors"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"is_modifier": True}

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

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
        geometry: InputGeometry = None,
    ):
        super().__init__(**{"Geometry": geometry})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        atoms = tree.outputs.geometry("Atoms")

        with g.Frame("Isolate verts pointing toward normals, get positions"):
            delete_geometry = g.DeleteGeometry.point(
                geometry,
                g.Compare.integer.not_equal(
                    g.NamedAttribute.integer("attribute_ID").o.attribute, 2
                ),
            )
            sample_index = g.SampleIndex(
                geometry=delete_geometry,
                value=g.Position(),
                index=g.Index(),
                data_type="FLOAT_VECTOR",
            )
        with g.Frame("Isolate verts pointing toward base vecs, get positions"):
            delete_geometry_1 = g.DeleteGeometry.point(
                geometry,
                g.Compare.integer.not_equal(
                    g.NamedAttribute.integer("attribute_ID").o.attribute, 1
                ),
            )
            sample_index_1 = g.SampleIndex(
                geometry=delete_geometry_1,
                value=g.Position(),
                index=g.Index(),
                data_type="FLOAT_VECTOR",
            )
        with g.Frame("Isolate verts containing bbone positions, get positions & geom"):
            delete_geometry_2 = g.DeleteGeometry.point(
                geometry,
                g.Compare.integer.not_equal(
                    g.NamedAttribute.integer("attribute_ID").o.attribute, 0
                ),
            )
            sample_index_2 = g.SampleIndex(
                geometry=delete_geometry_2,
                value=g.Position(),
                index=g.Index(),
                data_type="FLOAT_VECTOR",
            )
        with g.Frame("Recover base & normal vectors from vert positions, store"):
            (
                g.StoreNamedAttribute.point.vector(
                    delete_geometry_2,
                    name="base_vector",
                    value=(sample_index_1.o.value - sample_index_2).normalize(),
                )
                >> g.StoreNamedAttribute.point.vector(
                    name="base_normal",
                    value=(sample_index.o.value - sample_index_2).normalize(),
                )
                >> g.RemoveNamedAttribute(name="attribute_ID")
                >> atoms
            )


ASSET = OxDNARecoverVectors

ASSET_METADATA = {
    "catalog_id": "0094c3e0-7885-427b-81b4-187a84dcff18",
}
