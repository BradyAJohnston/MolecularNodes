import bpy
from bpy.app.handlers import persistent


# this update function requires a self and context input, as funcitons with these inputs
# have ot be passed to the `update` arguments of UI properties. When the UI is updated,
# the function is called with the UI element being `self` and the current context being
# passed into the function
def _udpate_entities(self, context: bpy.types.Context) -> None:
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
    for entity in session.entities:
        # use the updated method if it exists but otherwise fallback on the old method
        # of updating the trajectories
        
        if hasattr(entity, "update_with_scene"):
            if entity.update_with_scene:
                frame_to_set = scene.frame_current
            else:
                frame_to_set = entity.frame
            
            # do the entity setting
            try:
                entity.set_frame(scene.frame_current)
            except NotImplementedError:
                pass

        else:
            # this is the old method of updating the trajectories and is maintained for
            # backwards compatibility # TODO: takeout for later release
            entity._update_positions(scene.frame_current)
            entity._update_selections()
            entity._update_calculations()
