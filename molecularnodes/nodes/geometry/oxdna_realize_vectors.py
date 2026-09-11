# Node-group asset "oxDNA Realize Vectors" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .oxdna_vectors import OxDNAVectors


class OxDNARealizeVectors(AssetGeometryGroup):
    """
    oxDNA Realize Vectors

    Parameters
    ----------
    atoms : InputGeometry
        Geometry to duplicate elements of

    Inputs
    ------
    i.atoms : GeometrySocket
        Geometry to duplicate elements of

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "oxDNA Realize Vectors"
    _asset_name = "oxDNA Realize Vectors"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Geometry to duplicate elements of"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
    ):
        super().__init__(**{"Atoms": atoms})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Geometry to duplicate elements of"
        )
        geometry = tree.outputs.geometry("Geometry")

        value = g.Value(0.15)
        group = OxDNAVectors()
        for_each = g.ForEachGeometryElementZone(geometry=atoms)
        position = for_each.inputs.vector("Position", g.Position())
        base_vector = for_each.inputs.vector("base_vector", group.o.base_vector)
        base_normal = for_each.inputs.vector("base_normal", group.o.base_normal)
        with g.Frame("Create edges pointing along base and normal vectors"):
            mesh_line = g.MeshLine.end_points(
                2, position.output, position.output + base_vector.output * value
            )
            mesh_line_1 = g.MeshLine.end_points(
                2, position.output, position.output + base_normal.output * value
            )
        with g.Frame("Store attribute_ID on edge verts: 0=bbone, 1=base vec, 2=normal"):
            string = g.String(string="attribute_ID")
            string_1 = g.String(string="attribute_ID")
            sample_nearest = g.SampleNearest.point(mesh_line, position.output)
            store_named_attribute = g.StoreNamedAttribute.point.integer(
                mesh_line, g.Compare.integer.equal(g.Index(), sample_nearest), string
            )
            store_named_attribute_1 = g.StoreNamedAttribute.point.integer(
                store_named_attribute,
                g.Compare.integer.not_equal(g.Index(), sample_nearest),
                string_1,
                1,
            )
            sample_nearest_1 = g.SampleNearest.point(mesh_line_1, position.output)
            store_named_attribute_2 = g.StoreNamedAttribute.point.integer(
                mesh_line_1,
                g.Compare.integer.equal(g.Index(), sample_nearest_1),
                string,
            )
            store_named_attribute_3 = g.StoreNamedAttribute.point.integer(
                store_named_attribute_2,
                g.Compare.integer.not_equal(g.Index(), sample_nearest_1),
                string_1,
                2,
            )
        join_geometry = g.JoinGeometry(
            geometry=(store_named_attribute_3, store_named_attribute_1)
        )
        with g.Frame("Transfer all attrib to new verts to retain thru Merge by Dist"):
            integer = g.Integer(integer=0)
            transfer_attributes = g.TransferAttributes(
                target=join_geometry,
                target_point_id=integer,
                source=for_each.input.o.element,
                source_point_id=integer,
                attribute_names=g.String(string="position,attribute_ID").o.string.split(
                    ","
                ),
                pattern_mode="Exact",
                exclude_names=True,
            )
        transfer_attributes >> for_each.generation.input
        store_named_attribute_4 = (
            for_each.output
            >> g.StoreNamedAttribute.point.integer(name=g.String(string="attribute_ID"))
        )
        (
            g.JoinGeometry(
                geometry=(for_each.generation.output, store_named_attribute_4)
            )
            >> g.MergeByDistance(distance=0.001)
            >> geometry
        )


ASSET = OxDNARealizeVectors

ASSET_METADATA = {
    "catalog_id": "0094c3e0-7885-427b-81b4-187a84dcff18",
}
