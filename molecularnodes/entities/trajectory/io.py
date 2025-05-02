"""
Importing molecular dynamics trajectories and associated files.
"""

__name__ = "MolecularNodes.trajectory"
__author__ = "Brady Johnston"

from pathlib import Path
import MDAnalysis as mda
from . import oxdna
from .base import Trajectory


def load(
    top: str | Path,
    traj: str | Path,
    name: str = "NewTrajectory",
    style: str | None = "spheres",
):
    universe = mda.Universe(top, traj)
    trajectory = Trajectory(universe=universe)
    trajectory.create_object(name=name, style=style)

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
    traj = oxdna.OXDNA(univ, world_scale=world_scale)
    traj.create_object(name=name, style=style)
    return traj
