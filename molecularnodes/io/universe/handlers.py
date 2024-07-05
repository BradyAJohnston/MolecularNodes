import bpy
from bpy.app.handlers import persistent


# this update function requires a self and context input, as funcitons with these inputs
# have ot be passed to the `update` arguments of UI properties. When the UI is updated,
# the function is called with the UI element being `self` and the current context being
# passed into the function
def _update_universes(self, context: bpy.types.Context) -> None:
    """
    Function for being called at various points in the updating of the UI, to ensure
    positions and selections of the universes are udpated with the new inputs
    """
    update_universes(context.scene)


def _selection_update_universes(self, context: bpy.types.Context) -> None:
    """
    Function for selection changing. If the selection is immutable e.g.
    when it is generated from an AtomGroup,
    update universe will not be called to avoid invalid selection str.
    """
    if self.immutable:
        return
    else:
        update_universes(context.scene)


# this is the 'perisisent' function which can be appended onto the
# `bpy.app.handlers.frame_change_*` functions. Either before or after the frame changes
# this function will then be called - ensuring all of the universes are up to date. We use
# the `frame_change_post` handler as we want the frame to change first, then we update the
# universe based on the current frame value
@persistent
def update_universes(scene):
    "Updatins all positions and selections for each universe."
    for universe in scene.MNSession.universes.values():
        try:
            universe._update_trajectory(scene.frame_current)
            universe._update_selections()
        except Exception as e:
            print(f"Error updating {universe}: {e}")
