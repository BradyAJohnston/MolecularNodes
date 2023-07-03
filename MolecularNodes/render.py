import bpy
from . import md
from . import nodes
import tempfile

def render_universe(universe, starting_style = 2, engine = "eevee", fstop = 2):
    from IPython.display import Image
    obj, frames = md.load_universe(universe)
    nodes.create_starting_node_tree(obj, frames, starting_style)
    with tempfile.NamedTemporaryFile(suffix='.png', delete=True) as temp_file:
        filename = temp_file.name
        render_object(obj, filename, engine=engine, fstop = fstop)
        return Image(filename)

def camera_dof(obj = None, fstop = 2):
    cam = bpy.context.scene.camera
    cam.data.dof.use_dof = True
    
    if obj:
        diff = obj.location - cam.location
        cam.data.dof.focus_distance = diff.length / 2
        cam.data.dof.aperture_fstop = fstop

def default_settings_cycles(
    samples = 128,
    resolution_x = 1280, 
    resolution_y = 720, 
    device = "GPU",
    dof = True
    ):
    scene = bpy.context.scene
    scene.render.engine = "CYCLES"
    scene.cycles.device = device
    scene.cycles.samples = samples
    scene.render.resolution_x = resolution_x
    scene.render.resolution_y = resolution_y

def default_settings_eevee(
    ao: bool = True, 
    ao_dis: float = 2, 
    bloom: bool = True, 
    resolution_x = 1280, 
    resolution_y = 720
    ):
    """
    Set useful default render settings for Eevee.

    Parameters:
        ao (bool): Whether to enable ambient occlusion. Defaults to True.
        ao_dis (float): The distance for ambient occlusion. Defaults to 2.
        bloom (bool): Whether to enable bloom. Defaults to True.
    """
    
    eevee = bpy.context.scene.eevee
    scene = bpy.context.scene
    
    scene.render.resolution_x = resolution_x
    scene.render.resolution_y = resolution_y

    eevee.use_gtao = ao
    eevee.gtao_distance = ao_dis
    eevee.gtao_factor = 2
    eevee.use_bloom = bloom

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



def render_object(object: bpy.types.Object, file_path: str = None, engine = "eevee", fstop = 2):
    # remove default cube for rendering
    objects = bpy.data.objects
    if objects.get('Cube'):
        objects.remove(objects['Cube'])
    object.select_set(True)
    bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='MEDIAN')
    frame_object(object)
    camera_dof(object, fstop = fstop)
    if engine == "eevee":
        # frame the camera and render the image
        default_settings_eevee()
    elif engine == "cycles":
        default_settings_cycles()
    render(file_path)