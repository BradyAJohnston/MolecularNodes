"""
Helper classes for trajectory data management and Blender integration.

This module contains utility classes that support the Trajectory class by handling:
- Position caching and management
- Frame calculations and interpolation
- Blender property synchronization
- Attribute metadata collection
"""

import logging
from collections import OrderedDict
from dataclasses import dataclass
from typing import Any, Callable, Dict, Protocol
import bpy
import databpy as db
import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from ...utils import (
    correct_periodic_positions,
    fraction,
    frames_to_average,
)

logger = logging.getLogger(__name__)


# ============================================================================
# Protocols and Type Definitions
# ============================================================================


class ComputeAttributeFunc(Protocol):
    """Protocol for attribute computation functions."""

    def __call__(self) -> np.ndarray: ...


# ============================================================================
# Attribute Management
# ============================================================================


@dataclass
class AttributeSpec:
    """
    Metadata specification for a Blender attribute.

    Attributes:
        name: The attribute name in Blender
        compute_fn: Optional function to compute the attribute values
        selection_str: Optional MDAnalysis selection string
        dtype: Python type for the attribute data
        domain: Blender attribute domain (POINT, EDGE, etc.)
    """

    name: str
    compute_fn: Callable[[], np.ndarray] | None = None
    selection_str: str | None = None
    dtype: type = float
    domain: str = "POINT"

    def compute(self, universe: mda.Universe | None = None) -> np.ndarray:
        """
        Compute attribute values.

        Args:
            universe: MDAnalysis Universe for selection-based computation

        Returns:
            Computed attribute values as numpy array

        Raises:
            ValueError: If neither compute_fn nor selection_str is provided
        """
        if self.compute_fn is not None:
            return self.compute_fn()
        elif self.selection_str is not None and universe is not None:
            from .base import _ag_to_bool

            return _ag_to_bool(universe.select_atoms(self.selection_str))
        else:
            raise ValueError(
                f"Cannot compute attribute {self.name}: "
                "no compute function or selection string"
            )


class AttributeMetadata:
    """
    Collects metadata during attribute computation for batch application.

    This class accumulates metadata entries during attribute computation
    and applies them all at once to minimize Blender API calls.
    """

    def __init__(self):
        self.metadata: Dict[str, Any] = {}

    def add(self, key: str, value: Any) -> None:
        """
        Add a metadata entry.

        Args:
            key: Metadata key
            value: Metadata value (must be JSON-serializable)
        """
        self.metadata[key] = value

    def apply_to_object(self, obj: bpy.types.Object) -> None:
        """
        Batch apply all collected metadata to a Blender object.

        Args:
            obj: Target Blender object
        """
        for key, value in self.metadata.items():
            try:
                obj[key] = value
            except (TypeError, db.LinkedObjectError) as e:
                logger.warning(f"Failed to set metadata '{key}' on object: {e}")

    def clear(self) -> None:
        """Clear all collected metadata."""
        self.metadata.clear()


# ============================================================================
# Blender Property Management
# ============================================================================


class BlenderProperty:
    """
    Descriptor for properties stored on Blender objects with optional validation.

    This descriptor provides a clean interface for accessing Blender custom
    properties while supporting validation and type checking.

    Example:
        >>> class MyClass:
        ...     frame = BlenderProperty("frame", validate_fn=lambda x: x >= 0)
    """

    def __init__(
        self, attr_name: str, validate_fn: Callable[[Any], None] | None = None
    ):
        """
        Initialize the descriptor.

        Args:
            attr_name: Name of the property on the Blender object
            validate_fn: Optional validation function that raises on invalid values
        """
        self.attr_name = attr_name
        self.validate_fn = validate_fn

    def __set_name__(self, owner, name):
        """Store the Python attribute name."""
        self.name = name

    def __get__(self, obj, objtype=None):
        """Get property value from Blender object."""
        if obj is None:
            return self
        return getattr(obj.object.mn, self.attr_name)

    def __set__(self, obj, value):
        """Set property value on Blender object with validation."""
        if self.validate_fn is not None:
            self.validate_fn(value)
        setattr(obj.object.mn, self.attr_name, value)


class BlenderPropertyBridge:
    """
    Handles bidirectional data exchange between trajectory state and Blender properties.

    This class centralizes the logic for reading and writing trajectory
    properties to/from Blender objects.
    """

    @staticmethod
    def sync_to_blender(obj: bpy.types.Object, properties: Dict[str, Any]) -> None:
        """
        Write trajectory properties to a Blender object.

        Args:
            obj: Target Blender object
            properties: Dictionary of property names and values
        """
        for key, value in properties.items():
            try:
                setattr(obj.mn, key, value)
            except AttributeError as e:
                logger.warning(f"Failed to set property '{key}' on Blender object: {e}")

    @staticmethod
    def sync_from_blender(
        obj: bpy.types.Object, property_names: list[str]
    ) -> Dict[str, Any]:
        """
        Read trajectory properties from a Blender object.

        Args:
            obj: Source Blender object
            property_names: List of property names to read

        Returns:
            Dictionary of property names and their values
        """
        properties = {}
        for name in property_names:
            try:
                properties[name] = getattr(obj.mn, name)
            except AttributeError:
                logger.warning(f"Property '{name}' not found on Blender object")
        return properties


# ============================================================================
# Position Cache Management
# ============================================================================


class PositionCache:
    """
    Manages position caching with automatic size limits and cleanup.

    This class provides an OrderedDict-based cache for trajectory positions
    with automatic eviction of old entries to maintain memory efficiency.
    """

    def __init__(self, max_size: int = 10):
        """
        Initialize the position cache.

        Args:
            max_size: Maximum number of frames to cache
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
        """
        Get cached position or compute and cache it.

        Args:
            frame: Frame number
            compute_fn: Function to compute positions if not cached

        Returns:
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
        """
        Remove all frames except the specified ones.

        Args:
            frames_to_keep: Array of frame numbers to retain
        """
        to_remove = [f for f in self._cache if f not in frames_to_keep]
        for f in to_remove:
            del self._cache[f]

    def get_ordered_array(self) -> np.ndarray:
        """
        Return cached frames as a 3D array in chronological order.

        Returns:
            3D numpy array of shape (n_frames, n_atoms, 3)
        """
        keys = list(self._cache.keys())
        keys.sort()
        return np.array([self._cache[k] for k in keys])


# ============================================================================
# Frame Management
# ============================================================================


class FrameManager:
    """
    Manages frame updates, position caching, interpolation, and periodic corrections.

    This class encapsulates all logic related to:
    - Position caching and retrieval
    - Frame-to-frame interpolation
    - Periodic boundary condition corrections
    - Frame averaging

    The FrameManager is designed to work with a Trajectory instance and
    handles the complexity of mapping between scene frames and universe frames.
    """

    def __init__(self, trajectory):
        """
        Initialize the frame manager.

        Args:
            trajectory: Parent Trajectory instance
        """
        self.trajectory = trajectory
        self.cache = PositionCache()

    def _position_at_frame(self, frame: int) -> np.ndarray:
        """
        Get atom positions at a specific universe frame.

        Args:
            frame: Universe frame number

        Returns:
            Scaled atom positions as numpy array
        """
        self.trajectory.uframe = frame
        return self.trajectory.univ_positions

    def adjust_periodic_positions(
        self, pos1: np.ndarray, pos2: np.ndarray
    ) -> np.ndarray:
        """
        Apply periodic boundary correction to positions.

        Corrects for atoms that have crossed periodic boundaries between frames,
        ensuring smooth interpolation across boundaries.

        Args:
            pos1: Reference positions
            pos2: Positions to correct

        Returns:
            Corrected positions (or unchanged if correction not needed)
        """
        if self.trajectory.correct_periodic and self.trajectory.is_orthorhombic:
            return correct_periodic_positions(
                pos1, pos2, self.trajectory.universe.dimensions[:3]
            )
        else:
            return pos2

    def _frame_range(self, frame: int) -> npt.NDArray[np.int64]:
        """
        Get frame numbers to average over.

        Args:
            frame: Center frame number

        Returns:
            Array of frame numbers to include in average
        """
        return frames_to_average(
            frame, self.trajectory.n_frames, average=self.trajectory.average
        )

    def update_position_cache(self, frame: int, cache_ahead: bool = True) -> None:
        """
        Update the position cache for the current frame.

        This method intelligently caches positions based on averaging and
        interpolation settings, prefetching the next frame when needed.

        Args:
            frame: Current frame number
            cache_ahead: If True, cache next frame for interpolation
        """
        frames_to_cache = self._frame_range(frame)

        # If interpolating, ensure we cache 1 frame ahead
        if (
            len(frames_to_cache) == 1
            and frames_to_cache[0] != (self.trajectory.n_frames - 1)
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
        """
        Get mean position from cached frames.

        Computes the average position over the frames specified by the
        averaging window, applying periodic boundary corrections if needed.

        Args:
            frame: Center frame number

        Returns:
            Mean positions (or single frame if not averaging)
        """
        self.update_position_cache(frame)

        if self.trajectory.average == 0:
            return self.cache[frame]

        array = self.cache.get_ordered_array()
        if self.trajectory.correct_periodic and self.trajectory.is_orthorhombic:
            # Correct periodic boundary crossing relative to the first frame
            for i, pos in enumerate(array):
                if i == 0:
                    continue
                array[i] = self.adjust_periodic_positions(array[0], pos)

        return np.mean(array, axis=0)

    def get_positions_at_frame(self, frame: int) -> np.ndarray:
        """
        Get positions for a given frame with all processing applied.

        This is the main entry point for position retrieval, handling:
        - Frame mapping (offset, subframes)
        - Interpolation between frames
        - Averaging over multiple frames
        - Periodic boundary corrections

        Args:
            frame: Scene frame number

        Returns:
            Processed atom positions
        """
        if not self.trajectory.update_with_scene:
            # Just return positions at the frame without any special handling
            return self._position_at_frame(frame)

        # Map scene frame to universe frame
        uframe_current = self.trajectory.frame_mapper(frame)
        uframe_next = uframe_current + 1
        last_frame = self.trajectory.n_frames - 1

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
                and self.trajectory.is_orthorhombic
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


# ============================================================================
# Validation Functions
# ============================================================================


def _validate_non_negative(value: int) -> None:
    """
    Validation function for non-negative integers.

    Args:
        value: Value to validate

    Raises:
        ValueError: If value is negative
    """
    if value < 0:
        raise ValueError(f"Value must be non-negative, got {value}")


def _validate_frame(value: int) -> None:
    """
    Validation function for frame numbers.

    Args:
        value: Frame number to validate

    Raises:
        ValueError: If frame is negative
    """
    if value < 0:
        raise ValueError(f"Frame must be non-negative, got {value}")
