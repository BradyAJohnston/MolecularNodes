import os
import tempfile

os.environ.setdefault("BLENDER_USER_EXTENSIONS", tempfile.mkdtemp())
import bpy
from nodebpy import TreeBuilder
from nodebpy import geometry as g

bpy.ops.wm.read_homefile(use_empty=True)
scene = bpy.context.scene
mat = bpy.data.materials.new("gp")
bpy.data.materials.create_gpencil_data(mat)
mat.grease_pencil.show_fill = True
mat.grease_pencil.fill_color = (1, 0.5, 0, 1)
print(
    "gp mat props",
    [
        p.identifier
        for p in mat.grease_pencil.bl_rna.properties
        if "fill" in p.identifier or "show" in p.identifier or "node" in p.identifier
    ],
)
print("mat use_nodes", mat.use_nodes, mat.node_tree)
with TreeBuilder.geometry("probe", is_modifier=True) as tree:
    o = tree.outputs.geometry("Geometry")
    c = g.CurveCircle(radius=1.0) >> g.SetCurveRadius(radius=0.05)
    (
        g.CurvesToGreasePencil(curves=c, instances_as_layers=False)
        >> g.SetMaterial(material=mat)
        >> o
    )
ob = bpy.data.objects.new("GP", bpy.data.grease_pencils.new("gpd"))
scene.collection.objects.link(ob)
ob.data.materials.append(mat)
ob.modifiers.new("GN", "NODES").node_group = bpy.data.node_groups["probe"]
ev = ob.evaluated_get(bpy.context.evaluated_depsgraph_get())
print("eval mats", [m and m.name for m in ev.data.materials])
d = ev.data.layers[0].current_frame().drawing
print("mat idx", [s.material_index for s in d.strokes], [s.cyclic for s in d.strokes])
cam = bpy.data.objects.new("cam", bpy.data.cameras.new("cam"))
scene.collection.objects.link(cam)
cam.location = (0, 0, 10)
scene.camera = cam
scene.render.engine = "BLENDER_EEVEE"
scene.render.resolution_x = scene.render.resolution_y = 200
scene.render.filepath = os.path.join(os.path.dirname(__file__), "probe4.png")
bpy.ops.render.render(write_still=True)
