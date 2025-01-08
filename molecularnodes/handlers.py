import bpy
from bpy.app.handlers import persistent
# from .session import update_trajectories

# from .session import MNSession


# this update function requires a self and context input, as funcitons with these inputs
# have ot be passed to the `update` arguments of UI properties. When the UI is updated,
# the function is called with the UI element being `self` and the current context being
# passed into the function
def _update_trajectories(self, context: bpy.types.Context) -> None:
    """
    Function for being called at various points in the updating of the UI, to ensure
    positions and selections of the trajectories are udpated with the new inputs
    """
    update_trajectories(context.scene)


def _update_trajectories_on_frame_change(self, context: bpy.types.Context) -> None:
    """
    Function for being called at various points in the updating of the UI, to ensure
    positions and selections of the trajectories are udpated with the new inputs
    """
    update_trajectories(context.scene)


def _selection_update_trajectories(self, context: bpy.types.Context) -> None:
    """
    Function for selection changing. If the selection is immutable e.g.
    when it is generated from an AtomGroup,
    update universe will not be called to avoid invalid selection str.
    """
    if self.immutable:
        return
    else:
        update_trajectories(context.scene)


# this is the 'perisisent' function which can be appended onto the
# `bpy.app.handlers.frame_change_*` functions. Either before or after the frame changes
# this function will then be called - ensuring all of the trajectories are up to date. We
# use the `frame_change_pre` handler as we want the frame to change first, then we update
# the universe based on the current frame value
@persistent
def update_trajectories(scene):
    "Call the set_frame method of all trajectories in the current session"
    session = scene.MNSession
    for traj in session.trajectories.values():
        # use the updated method if it exists but otherwise fallback on the old method
        # of updating the trajectories
        if hasattr(traj, "update_with_scene"):
            if traj.update_with_scene:
                traj.set_frame(scene.frame_current)
            else:
                traj.set_frame(traj.frame)

        else:
            traj._update_positions(scene.frame_current)
            traj._update_selections()
            traj._update_calculations()

        # except NotImplementedError:
        #     pass
        # except Exception as e:
        #     # print(f"Error updating {traj}: {e}")
        #     raise e
