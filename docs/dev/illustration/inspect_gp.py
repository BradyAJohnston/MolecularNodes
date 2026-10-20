import os
import sys
import tempfile

os.environ.setdefault("BLENDER_USER_EXTENSIONS", tempfile.mkdtemp())
import bpy
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

bpy.ops.wm.open_mainfile(filepath=sys.argv[1])
from illustrate import STYLES, illustration  # noqa: E402

mol_ob = bpy.data.objects[sys.argv[3]]
print(
    "mol obj",
    mol_ob,
    len(mol_ob.evaluated_get(bpy.context.evaluated_depsgraph_get()).data.vertices),
)
cam = bpy.context.scene.camera


class M:
    pass


m = M()
m.object = mol_ob
m.name = sys.argv[3]
ob = illustration(m, cam, STYLES[sys.argv[2]])
ev = ob.evaluated_get(bpy.context.evaluated_depsgraph_get())
for layer in ev.data.layers:
    d = layer.current_frame().drawing
    print("strokes", len(d.strokes))
    pos = np.array([p.position for s in d.strokes for p in s.points])
    rad = np.array([p.radius for s in d.strokes for p in s.points])
    print("pos range", pos.min(0), pos.max(0), "radius", rad.min(), rad.max())
    for s in list(d.strokes)[:3]:
        print(
            len(s.points),
            s.material_index,
            s.cyclic,
            [tuple(round(x, 3) for x in p.position) for p in s.points][:2],
            s.fill_color[:],
            s.points[0].vertex_color[:],
        )
    print({a.name: a.domain for a in d.attributes})
print(
    "slots",
    [m.name for m in ob.data.materials],
    [
        (
            m.name,
            m.grease_pencil.show_fill,
            m.grease_pencil.fill_style,
            m.grease_pencil.fill_color[:],
        )
        for m in ob.data.materials
    ],
)
