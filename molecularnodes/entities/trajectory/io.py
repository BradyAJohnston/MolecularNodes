"""
Importing molecular dynamics trajectories and associated files.
"""

__name__ = "MolecularNodes.trajectory"
__author__ = "Brady Johnston"

from pathlib import Path
import MDAnalysis as mda
from . import StreamingTrajectory, imd, oxdna
from .base import Trajectory


def load(
    top: str | Path,
    traj: str | Path,
    name: str = "NewTrajectory",
    style: str | None = "spheres",
    selection: str | None = None,
):
    """Load a trajectory with automatic format detection.

    Detects and loads trajectories in various formats including standard
    MD trajectories and real-time IMD streaming connections.

    Parameters
    ----------
    top : str | Path
        Path to topology file
    traj : str | Path
        Path to trajectory file or IMD URL (e.g., 'imd://localhost:3000')
    name : str, optional
        Name for the trajectory object, by default "NewTrajectory"
    style : str | None, optional
        Visual style to apply, by default "spheres"
    selection : str | None, optional
        Atom selection string, by default None (all atoms)

    Returns
    -------
    Trajectory
        Loaded trajectory object (may be StreamingTrajectory for IMD)
    """
    if imd.is_imd_url(traj):
        trajectory = StreamingTrajectory.load(top, traj, name=name)
    else:
        universe = mda.Universe(top, traj)
        trajectory = Trajectory(universe=universe, name=name)

    trajectory.add_style(style=style, selection=selection)

    return trajectory


def load_oxdna(top, traj, name="oxDNA", style="oxdna", world_scale=0.01):
    """
    Load an oxDNA trajectory.

    Parameters
    ----------
    top : str
        Path to topology file
    traj : str
        Path to trajectory file
    name : str, optional
        Name for the created object, by default "oxDNA"
    style : str, optional
        Style of representation, by default "oxdna"
    world_scale : float, optional
        Scaling factor for world coordinates, by default 0.01

    Returns
    -------
    OXDNA
        The created trajectory object
    """
    univ = mda.Universe(
        top, traj, topology_format=oxdna.OXDNAParser, format=oxdna.OXDNAReader
    )
    traj = oxdna.OXDNA(univ, name=name, world_scale=world_scale).add_style(style=style)
    return traj
