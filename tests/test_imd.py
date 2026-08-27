import MDAnalysis as mda
from molecularnodes.entities.molecule.imd import StreamingTrajectory
from .constants import data_dir


def _streaming_trajectory() -> StreamingTrajectory:
    # a regular multi-frame universe stands in for an IMD stream, as `.next()`
    # behaves the same and lets the test observe the reader's frame
    universe = mda.Universe(
        data_dir / "md_ppr/box.gro",
        data_dir / "md_ppr/first_5_frames.xtc",
    )
    return StreamingTrajectory(universe, name="stream")


def test_stream_advances_once_per_scene_frame():
    traj = _streaming_trajectory()
    assert traj.universe.trajectory.frame == 0

    traj.set_frame(1)
    assert traj.universe.trajectory.frame == 1

    # repeated updates for the same scene frame (render handlers, UI property
    # changes) must not advance the stream
    traj.set_frame(1)
    traj.set_frame(1)
    assert traj.universe.trajectory.frame == 1

    traj.set_frame(2)
    assert traj.universe.trajectory.frame == 2
