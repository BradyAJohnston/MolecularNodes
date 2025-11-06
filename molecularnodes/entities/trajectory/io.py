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
    selection: str | None = None,
):
    universe = mda.Universe(top, traj)
    trajectory = Trajectory(universe=universe, name=name).add_style(
        style=style, selection=selection
    )

    return trajectory


def load_oxdna(top, traj, name="oxDNA", style="oxdna"):
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

    Returns
    -------
    OXDNA
        The created trajectory object

    Notes
    -----
    World scale is now read from the global scene property bpy.context.scene.mn.world_scale.
    """
    univ = mda.Universe(
        top, traj, topology_format=oxdna.OXDNAParser, format=oxdna.OXDNAReader
    )
    traj = oxdna.OXDNA(univ, name=name).add_style(style=style)
    return traj
