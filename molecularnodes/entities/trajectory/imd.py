"""Interactive Molecular Dynamics (IMD) streaming trajectory support.

Provides support for real-time trajectory streaming from molecular dynamics
simulations using the IMD protocol via MDAnalysis and imdclient.
"""

import logging
from pathlib import Path
from typing import Optional
import databpy as db
import MDAnalysis as mda
from bpy.types import Object
from ...blender import set_obj_active
from ..base import EntityType
from .base import Trajectory

logger = logging.getLogger(__name__)


def normalize_imd_url(url: str) -> str:
    """Normalize IMD URL format.

    Converts 'imd:/host:port' to 'imd://host:port' to ensure
    proper URL formatting for MDAnalysis.

    Parameters
    ----------
    url : str
        IMD URL in any accepted format

    Returns
    -------
    str
        Normalized IMD URL with proper protocol formatting
    """
    url_str = str(url)
    if url_str.startswith("imd:/") and not url_str.startswith("imd://"):
        return url_str.replace("imd:/", "imd://", 1)
    return url_str


def is_imd_url(path: str | Path) -> bool:
    """Check if path is an IMD streaming URL.

    Parameters
    ----------
    path : str | Path
        Path or URL to check

    Returns
    -------
    bool
        True if path is an IMD streaming URL
    """
    return str(path).startswith("imd:")


class StreamingTrajectory(Trajectory):
    """Trajectory subclass for IMD streaming connections.

    Handles real-time molecular dynamics trajectories streamed via the IMD
    protocol. Unlike regular trajectories with random frame access, streaming
    trajectories read frames sequentially as they become available.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        MDAnalysis Universe with IMD reader
    name : str, optional
        Name for the trajectory object, by default "StreamingTrajectory"
    create_object: bool
        Create the mesh object to be linked to this trajectory
    """

    _entity_type = EntityType.MD_STREAMING

    def __init__(
        self,
        universe: mda.Universe,
        name: str = "StreamingTrajectory",
        create_object=True,
    ):
        super().__init__(universe, name=name, create_object=create_object)

    @property
    def n_frames(self) -> Optional[int]:
        """Return None for streaming trajectories (unknown frame count).

        Streaming trajectories don't have a predetermined number of frames
        since new frames arrive continuously from the simulation.

        Returns
        -------
        Optional[int]
            None, as frame count is not known in advance
        """
        return None

    def _update_trajectory_positions(self, frame: int) -> None:
        """Override to handle streaming frame updates.

        For streaming trajectories, we ignore the requested frame number
        and simply advance to the next available frame from the stream.

        Parameters
        ----------
        frame : int
            Ignored for streaming trajectories

        Raises
        ------
        StopIteration
            If the stream has ended
        Exception
            If there's an error reading from the IMD stream
        """
        try:
            self.universe.trajectory.next()
            self.position = self._scaled_position
        except StopIteration:
            logger.warning("Stream ended or connection lost")
            raise
        except Exception as e:
            logger.error(f"Error reading IMD stream: {e}")
            raise

    @classmethod
    def load(
        cls,
        topology: Path | str,
        coordinates: Path | str,
        name: str = "StreamingTrajectory",
        style: str | None = "spheres",
        selection: str | None = None,
        create_object: bool = True,
    ):
        url = normalize_imd_url(str(coordinates))
        logger.info(f"Attempting IMD connection to {url}")
        try:
            u = mda.Universe(topology, url)
            logger.info("IMD connection established successfully: {}".format(u))
        except Exception as e:
            logger.error(f"Failed to create IMD connection to {url}: {e}")
            raise ValueError(f"Could not connect to IMD server: {e}") from e
        traj = cls(u, name=name, create_object=create_object)
        if style and create_object:
            traj.add_style(style=style, selection=selection)

        return traj
