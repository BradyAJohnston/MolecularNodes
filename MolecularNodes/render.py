import bpy

def frame_object(object: bpy.types.Object):
    # old_active_object = bpy.data.
    
    # set as active object
    object.select_set(True)
    
    # frame the active object
    bpy.ops.view3d.camera_to_view_selected()

def render(file_path: str = None):
    if not file_path:
        file_path = "/tmp/"
    
    bpy.context.scene.render.filepath = file_path
    bpy.ops.render.render(write_still = True)

def render_object(object: bpy.types.Object, file_path: str = None):
    frame_object(object)
    render(file_path)