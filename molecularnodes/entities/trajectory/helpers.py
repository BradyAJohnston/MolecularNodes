"""Helper classes for trajectory data management and Blender integration.

Provides utility classes for position caching, frame management, Blender property
synchronization, and attribute metadata collection.
"""

import logging
from collections import OrderedDict
from typing import Callable
import databpy as db
import numpy as np
import numpy.typing as npt
from MDAnalysis import AtomGroup
from ...utils import (
    correct_periodic_positions,
    fraction,
    frame_mapper,
    frames_to_average,
)

logger = logging.getLogger(__name__)


def _ag_to_bool(ag: AtomGroup) -> np.ndarray:
    """Convert AtomGroup to boolean mask for the entire universe."""
    return np.isin(ag.universe.atoms.ix, ag.ix).astype(bool)


# ============================================================================
# Position Cache Management
# ============================================================================


class PositionCache:
    """Position cache with automatic size limits and cleanup.

    OrderedDict-based cache with automatic eviction for memory efficiency.

    Parameters
    ----------
    max_size : int, default=10
        Maximum number of frames to cache
    """

    def __init__(self, max_size: int = 10):
        """Initialize the position cache.

        Parameters
        ----------
        max_size : int
            Maximum number of frames to cache
        """
        self._cache: OrderedDict[int, np.ndarray] = OrderedDict()
        self._max_size = max_size

    def __contains__(self, frame: int) -> bool:
        """Check if a frame is cached."""
        return frame in self._cache

    def __getitem__(self, frame: int) -> np.ndarray:
        """Get cached positions for a frame."""
        return self._cache[frame]

    def __setitem__(self, frame: int, positions: np.ndarray) -> None:
        """Cache positions for a frame."""
        self._cache[frame] = positions
        self._enforce_size_limit()

    def __delitem__(self, frame: int) -> None:
        """Remove a frame from the cache."""
        del self._cache[frame]

    def keys(self):
        """Get all cached frame numbers."""
        return self._cache.keys()

    def _enforce_size_limit(self) -> None:
        """Remove oldest entries if cache exceeds max size."""
        while len(self._cache) > self._max_size:
            self._cache.popitem(last=False)

    def get_or_compute(
        self, frame: int, compute_fn: Callable[[int], np.ndarray]
    ) -> np.ndarray:
        """Get cached position or compute and cache it.

        Parameters
        ----------
        frame : int
            Frame number
        compute_fn : callable
            Function to compute positions if not cached

        Returns
        -------
        np.ndarray
            Atom positions for the frame
        """
        if frame not in self._cache:
            self._cache[frame] = compute_fn(frame)
            self._enforce_size_limit()
        return self._cache[frame]

    def clear(self) -> None:
        """Clear all cached positions."""
        self._cache.clear()

    def remove_frames_except(self, frames_to_keep: npt.NDArray[np.int64]) -> None:
        """Remove all frames except the specified ones.

        Parameters
        ----------
        frames_to_keep : np.ndarray
            Frame numbers to retain
        """
        to_remove = [f for f in self._cache if f not in frames_to_keep]
        for f in to_remove:
            del self._cache[f]

    def get_ordered_array(self) -> np.ndarray:
        """Return cached frames as a 3D array in chronological order.

        Returns
        -------
        np.ndarray
            Shape (n_frames, n_atoms, 3)
        """
        keys = list(self._cache.keys())
        keys.sort()
        return np.array([self._cache[k] for k in keys])


# ============================================================================
# Frame Management
# ============================================================================


class FrameManager:
    """Frame updates, position caching, interpolation, and periodic corrections.

    Encapsulates:
    - Position caching and retrieval
    - Frame-to-frame interpolation
    - Periodic boundary corrections
    - Frame averaging

    Maps between Blender scene frames and trajectory universe frames.

    Parameters
    ----------
    trajectory : Trajectory
        Parent Trajectory instance
    """

    def __init__(self, trajectory):
        """Initialize the frame manager.

        Parameters
        ----------
        trajectory : Trajectory
            Parent Trajectory instance
        """
        self.trajectory = trajectory
        self.cache = PositionCache()

    @property
    def n_frames(self) -> int:
        return self.trajectory.universe.trajectory.n_frames

    def _position_at_frame(self, frame: int) -> np.ndarray:
        """Get atom positions at a specific universe frame.

        Parameters
        ----------
        frame : int
            Universe frame number

        Returns
        -------
        np.ndarray
            Scaled atom positions
        """
        self.trajectory.uframe = frame
        return self.trajectory._scaled_position

    def adjust_periodic_positions(
        self, pos1: np.ndarray, pos2: np.ndarray
    ) -> np.ndarray:
        """Apply periodic boundary correction to positions.

        Corrects for atoms crossing periodic boundaries, ensuring smooth interpolation.

        Parameters
        ----------
        pos1 : np.ndarray
            Reference positions
        pos2 : np.ndarray
            Positions to correct

        Returns
        -------
        np.ndarray
            Corrected positions (unchanged if correction not needed)
        """
        if self.trajectory.correct_periodic and self.trajectory._is_orthorhombic:
            return correct_periodic_positions(
                pos1, pos2, self.trajectory.universe.dimensions[:3]
            )
        else:
            return pos2

    def _frame_range(self, frame: int) -> npt.NDArray[np.int64]:
        """Get frame numbers to average over.

        Parameters
        ----------
        frame : int
            Center frame number

        Returns
        -------
        np.ndarray
            Frame numbers to include in average
        """
        return frames_to_average(
            frame,
            self.n_frames,
            average=self.trajectory.average,
        )

    def update_position_cache(self, frame: int, cache_ahead: bool = True) -> None:
        """Update the position cache for the current frame.

        Intelligently caches based on averaging and interpolation settings,
        prefetching the next frame when needed.

        Parameters
        ----------
        frame : int
            Current frame number
        cache_ahead : bool, default=True
            Whether to cache next frame for interpolation
        """
        frames_to_cache = self._frame_range(frame)

        # If interpolating, ensure we cache 1 frame ahead
        if (
            len(frames_to_cache) == 1
            and frames_to_cache[0] != (self.n_frames - 1)
            and cache_ahead
        ):
            frames_to_cache = np.array(
                (frames_to_cache[0], frames_to_cache[0] + 1), dtype=int
            )

        # Only cleanup the cache if we have more than 2 frames stored
        if len(self.cache.keys()) > 2:
            self.cache.remove_frames_except(frames_to_cache)

        # Update the cache with any frames that are not yet cached
        for f in frames_to_cache:
            if f not in self.cache:
                self.cache[f] = self._position_at_frame(f)

    def position_cache_mean(self, frame: int) -> np.ndarray:
        """Get mean position from cached frames.

        Computes average over the averaging window with periodic corrections if needed.

        Parameters
        ----------
        frame : int
            Center frame number

        Returns
        -------
        np.ndarray
            Mean positions (or single frame if not averaging)
        """
        self.update_position_cache(frame)

        if self.trajectory.average == 0:
            return self.cache[frame]

        array = self.cache.get_ordered_array()
        if self.trajectory.correct_periodic and self.trajectory._is_orthorhombic:
            # Correct periodic boundary crossing relative to the first frame
            for i, pos in enumerate(array):
                if i == 0:
                    continue
                array[i] = self.adjust_periodic_positions(array[0], pos)

        return np.mean(array, axis=0)

    def get_positions_at_frame(self, frame: int) -> np.ndarray:
        """Get positions for a given frame with all processing applied.

        Main entry point for position retrieval. Handles frame mapping, interpolation,
        averaging, and periodic boundary corrections.

        Parameters
        ----------
        frame : int
            Scene frame number

        Returns
        -------
        np.ndarray
            Processed atom positions
        """
        if not self.trajectory.update_with_scene:
            # Just return positions at the frame without any special handling
            return self._position_at_frame(frame)

        # Map scene frame to universe frame
        uframe_current = frame_mapper(
            frame=frame,
            subframes=self.trajectory.subframes,
            offset=self.trajectory.offset,
        )
        uframe_next = uframe_current + 1
        last_frame = self.n_frames - 1

        if uframe_current >= last_frame:
            uframe_current = last_frame
            uframe_next = uframe_current

        # Update the frame_hidden property for the UI
        try:
            self.trajectory._frame = uframe_current
        except AttributeError:
            # Silently ignore if we can't write to Blender property
            pass

        if self.trajectory.subframes > 0 and self.trajectory.interpolate:
            # Interpolate between current and next frame
            pos_current = self.position_cache_mean(uframe_current)
            pos_next = self.position_cache_mean(uframe_next)

            # Apply periodic correction if needed (and not already applied via averaging)
            if (
                self.trajectory.correct_periodic
                and self.trajectory._is_orthorhombic
                and self.trajectory.average == 0
            ):
                pos_next = correct_periodic_positions(
                    pos_current,
                    pos_next,
                    dimensions=self.trajectory.universe.dimensions[:3]
                    * self.trajectory.world_scale,
                )

            # Interpolate between the two sets of positions
            return db.lerp(
                pos_current, pos_next, t=fraction(frame, self.trajectory.subframes + 1)
            )
        elif self.trajectory.average > 0:
            # Return mean positions for cached frames
            return self.position_cache_mean(uframe_current)
        else:
            # Just return current positions
            return self._position_at_frame(uframe_current)
