import bpy
from . import color
import tempfile
import os
from mathutils import Vector

def setup_scene():
    # load template scene with nice HDRI lighting
    bpy.ops.wm.read_homefile(app_template = "MolecularNodes")

    # change the background to a custom color
    try:
        world_nodes = bpy.data.worlds['World Shader'].node_tree.nodes
        world_nodes['MN_world_shader'].inputs['BG Color'].default_value = color.random_rgb()
    except KeyError:
        print("Oh no, didn't set up the base scene.")

def clear_scene():
    bpy.ops.object.select_all(action="DESELECT")
    bpy.ops.object.select_by_type(type="MESH")
    bpy.ops.object.delete()
    for node in bpy.data.node_groups:
        if node.type == "GEOMETRY":
            bpy.data.node_groups.remove(node)

def orient_camera(object, lens = 85, dof = True, f = 2, zoom = 0, focus = 0):
    object.select_set(True)
    camera = bpy.data.objects['Camera']
    distance = Vector(object.location) - Vector(camera.location)
    camera.data.lens = lens
    camera.data.dof.aperture_fstop = f
    camera.data.dof.focus_object = None
    bpy.ops.view3d.camera_to_view_selected()
    camera.data.lens = lens + zoom
    camera.data.dof.use_dof = dof
    camera.data.dof.focus_distance = distance.length + focus
    # camera.data.dof.focus_distance = 1.2

def render_image(engine = 'eevee', x = 1000, y = 500):
    # setup render engine
    if engine == "eevee":
        bpy.context.scene.render.engine = "BLENDER_EEVEE"
    elif engine == "cycles":
        
        bpy.context.scene.render.engine = "CYCLES"
        try:
            bpy.context.scene.cycles.device = "GPU"
        except:
            print("GPU Rendering not available")
    

    # Render
    with tempfile.TemporaryDirectory() as temp:

        path = os.path.join(temp, "test.png")
        bpy.context.scene.render.resolution_x = x
        bpy.context.scene.render.resolution_y = y
        bpy.context.scene.render.image_settings.file_format = "PNG"
        bpy.context.scene.render.filepath = path
        bpy.ops.render.render(write_still=True)
        display(Image(filename=path))