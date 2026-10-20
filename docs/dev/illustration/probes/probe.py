import os
import sys
import tempfile

os.environ.setdefault("BLENDER_USER_EXTENSIONS", tempfile.mkdtemp())
import bpy
from nodebpy import TreeBuilder
from nodebpy import geometry as g

out = sys.argv[1]
bpy.ops.wm.read_homefile(use_empty=True)
scene = bpy.context.scene
mat = bpy.data.materials.new("gp")
bpy.data.materials.create_gpencil_data(mat)
mat.grease_pencil.show_fill = True
mat.grease_pencil.fill_color = (1, 0.5, 0, 1)

with TreeBuilder.geometry("probe", is_modifier=True) as tree:
    o = tree.outputs.geometry("Geometry")
    c = g.CurveCircle(radius=1.0) >> g.SetCurveRadius(radius=0.05)
    gp = g.CurvesToGreasePencil(curves=g.GeometryToInstance(c)) >> g.SetMaterial(
        material=mat
    )
    gp >> g.SetGreasePencilColor(color=(0, 0.2, 1, 1), mode="STROKE") >> o

for kind in ["MESH", "GREASEPENCIL"]:
    if kind == "MESH":
        data = bpy.data.meshes.new("m")
    else:
        data = bpy.data.grease_pencils.new("gpd")
    ob = bpy.data.objects.new(kind, data)
    scene.collection.objects.link(ob)
    ob.data.materials.append(mat)
    mod = ob.modifiers.new("GN", "NODES")
    mod.node_group = bpy.data.node_groups["probe"]
    dg = bpy.context.evaluated_depsgraph_get()
    ev = ob.evaluated_get(dg)
    print(kind, "eval type", ev.type, "data", type(ev.data).__name__)
    if kind == "MESH":
        ob.location.x = -3
    else:
        ob.location.x = 3
for o in bpy.context.view_layer.objects:
    print(o.name, o.type)
# export
for ob in scene.objects:
    ob.select_set(True)
bpy.context.view_layer.objects.active = scene.objects["GREASEPENCIL"]
cam = bpy.data.objects.new("cam", bpy.data.cameras.new("cam"))
scene.collection.objects.link(cam)
cam.location = (0, 0, 15)
scene.camera = cam
try:
    r = bpy.ops.wm.grease_pencil_export_svg(
        filepath=out + "/probe.svg", selected_object_type="VISIBLE", use_fill=True
    )
    print("svg", r)
except Exception as e:
    print("svg error", e)
scene.render.engine = "BLENDER_EEVEE"
scene.render.resolution_x = scene.render.resolution_y = 300
scene.render.filepath = out + "/probe.png"
bpy.ops.render.render(write_still=True)
