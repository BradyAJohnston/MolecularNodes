import os
import tempfile

os.environ.setdefault("BLENDER_USER_EXTENSIONS", tempfile.mkdtemp())
import bpy
from nodebpy import TreeBuilder
from nodebpy import geometry as g

bpy.ops.wm.read_homefile(use_empty=True)
for layers in (True, False):
    for wrap in (True, False):
        with TreeBuilder.geometry(f"p{layers}{wrap}", is_modifier=True) as tree:
            o = tree.outputs.geometry("Geometry")
            c = g.JoinGeometry((g.CurveCircle(radius=1.0), g.CurveLine()))
            if wrap:
                c = g.GeometryToInstance(c)
            g.CurvesToGreasePencil(curves=c, instances_as_layers=layers) >> o
        ob = bpy.data.objects.new("GP", bpy.data.grease_pencils.new("gpd"))
        bpy.context.scene.collection.objects.link(ob)
        ob.modifiers.new("GN", "NODES").node_group = bpy.data.node_groups[
            f"p{layers}{wrap}"
        ]
        print(
            "layers",
            layers,
            "wrap",
            wrap,
            ob.evaluated_get(
                bpy.context.evaluated_depsgraph_get()
            ).evaluated_geometry(),
        )
