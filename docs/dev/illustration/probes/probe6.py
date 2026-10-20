import math
import os
import tempfile

os.environ.setdefault("BLENDER_USER_EXTENSIONS", tempfile.mkdtemp())
import bpy

bpy.ops.wm.read_homefile(use_empty=True)
scene = bpy.context.scene
w = bpy.data.worlds.new("w")
scene.world = w
for p in ("show_fill", "show_stroke", "is_fill_visible", "is_stroke_visible"):
    print(p, bpy.types.MaterialGPencilStyle.bl_rna.properties[p].description)
combos = [(sf, ss) for sf in (True, False) for ss in (True, False)]
for k, (sf, ss) in enumerate(combos):
    mat = bpy.data.materials.new(f"m{k}")
    bpy.data.materials.create_gpencil_data(mat)
    mat.grease_pencil.show_fill = sf
    mat.grease_pencil.show_stroke = ss
    mat.grease_pencil.fill_color = (1, 0.5, 0, 1)
    mat.grease_pencil.color = (0, 0, 1, 1)
    print(
        k,
        sf,
        ss,
        "->",
        mat.grease_pencil.show_fill,
        mat.grease_pencil.show_stroke,
        mat.grease_pencil.is_fill_visible,
        mat.grease_pencil.is_stroke_visible,
    )
    gpd = bpy.data.grease_pencils.new("gpd")
    ob = bpy.data.objects.new("GP", gpd)
    scene.collection.objects.link(ob)
    gpd.materials.append(mat)
    d = gpd.layers.new("L").frames.new(1).drawing
    d.add_strokes([32])
    s = d.strokes[0]
    s.cyclic = True
    for i, p in enumerate(s.points):
        a = 2 * math.pi * i / 32
        p.position = (-4.5 + 3 * k + math.cos(a), math.sin(a), 0)
        p.radius = 0.05
cam = bpy.data.objects.new("cam", bpy.data.cameras.new("cam"))
scene.collection.objects.link(cam)
cam.location = (0, 0, 14)
scene.camera = cam
scene.render.engine = "BLENDER_EEVEE"
scene.render.resolution_x, scene.render.resolution_y = 400, 120
scene.render.filepath = os.path.join(os.path.dirname(__file__), "probe6.png")
bpy.ops.render.render(write_still=True)
