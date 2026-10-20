import os
import sys
import tempfile

os.environ.setdefault("BLENDER_USER_EXTENSIONS", tempfile.mkdtemp())
import bpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from illustrate import STYLES, illustration

bpy.ops.wm.open_mainfile(filepath=sys.argv[1])


class M:
    pass


m = M()
m.object = bpy.data.objects["CFF"]
m.name = "CFF"
ob = illustration(m, bpy.context.scene.camera, STYLES["classic"])
wrapper = ob.modifiers[0].node_group
tree = bpy.data.node_groups["Style Illustration"]
out = next(n for n in tree.nodes if n.bl_idname == "NodeGroupOutput")
cob = ob


for n in tree.nodes:
    if n.bl_idname == "NodeFrame":
        continue
    for s in n.outputs:
        if s.type == "GEOMETRY":
            tree.links.new(s, out.inputs[0])
            gs = cob.evaluated_get(
                bpy.context.evaluated_depsgraph_get()
            ).evaluated_geometry()
            print(f"{n.bl_idname:40s} {n.name:30s} {gs}")
for n in tree.nodes:
    for w in getattr(n, "warnings", []) or []:
        print("WARN", n.name, w)
n = next(n for n in tree.nodes if n.bl_idname == "GeometryNodeCurvesToGreasePencil")
for i in n.inputs:
    print(
        i.name,
        i.is_linked,
        getattr(i, "default_value", None),
        [link.from_node.name for link in i.links],
    )
