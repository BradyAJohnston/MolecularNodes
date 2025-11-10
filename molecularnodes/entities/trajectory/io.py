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
    if str(traj).startswith("imd:"):
        imd = str(traj).replace("imd:/", "imd://")
        universe = mda.Universe(top, imd, format="IMD")
        trajectory = Trajectory(universe, name=name)
        trajectory._mn_is_streaming = True
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
