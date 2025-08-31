import bpy
import numpy as np
from bpy.app.handlers import persistent
from PIL import Image
from .annotations.utils import render_annotations
from .scene.compositor import annotations_image, setup_compositor


# this update function requires a self and context input, as funcitons with these inputs
# have ot be passed to the `update` arguments of UI properties. When the UI is updated,
# the function is called with the UI element being `self` and the current context being
# passed into the function
def _update_entities(self, context: bpy.types.Context) -> None:
    """
    Function for being called at various points in the updating of the UI, to ensure
    positions and selections of the trajectories are udpated with the new inputs
    """
    update_entities(context.scene)


def _selection_update_trajectories(self, context: bpy.types.Context) -> None:
    """
    Function for selection changing. If the selection is immutable e.g.
    when it is generated from an AtomGroup,
    update universe will not be called to avoid invalid selection str.
    """
    if self.immutable:
        return
    else:
        update_entities(context.scene)


# this is the 'perisisent' function which can be appended onto the
# `bpy.app.handlers.frame_change_*` functions. Either before or after the frame changes
# this function will then be called - ensuring all of the trajectories are up to date. We
# use the `frame_change_pre` handler as we want the frame to change first, then we update
# the universe based on the current frame value
@persistent
def update_entities(scene):
    "Call the `set_frame()` method of all entities in the current session"
    session = scene.MNSession
    session.prune()

    for entity in session.entities.values():
        if entity.update_with_scene:
            frame_to_set = scene.frame_current
        else:
            frame_to_set = entity.frame

        # do the entity setting, if the method isn't implemented, just pass
        try:
            entity.set_frame(frame_to_set)
        except NotImplementedError:
            pass


@persistent
def render_pre_handler(scene: bpy.types.Scene) -> None:
    """
    Blender's on render (before) handler
    Any changes needed before the rendering of a frame need to go in here

    """
    if scene.mn.auto_setup_compositor:
        # Setup compositor if not already done
        setup_compositor(scene)
    # Render annotations to an image
    bpy.context.view_layer.update()
    render_scale = scene.render.resolution_percentage / 100
    width = int(scene.render.resolution_x * render_scale)
    height = int(scene.render.resolution_y * render_scale)
    image_scale = 2  # To workaround anti-aliasing issues with lines
    image = Image.new(
        "RGBA",
        (width * image_scale, height * image_scale),
        (0, 0, 0, 0),
    )
    # render annotations of all entities
    render_annotations(scene, render_scale, image, image_scale)
    if image_scale != 1:
        # scale down to actual render size
        image = image.resize((width, height))
    # create blender image
    if annotations_image not in bpy.data.images:
        bpy.data.images.new(annotations_image, width, height)
    bpy_image = bpy.data.images[annotations_image]
    bpy_image.scale(width, height)
    # update from PIL image
    bpy_image.pixels[:] = (np.flipud(np.array(image).astype(float)) / 255.0).ravel()
