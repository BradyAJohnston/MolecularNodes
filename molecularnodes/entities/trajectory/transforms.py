"""
Custom MDAnalysis Transform classes for subframe interpolation.

This module provides transform classes that enable stepping between frames of existing
trajectories by generating interpolated subframes using different interpolation methods.
"""

import math
from typing import Callable, Literal
import numpy as np
import numpy.typing as npt
try:
    from MDAnalysis.transformations.base import TransformationBase
except ImportError:
    # Fallback for older MDAnalysis versions
    class TransformationBase:
        def __init__(self, max_threads=1, parallelizable=True):
            self.max_threads = max_threads
            self.parallelizable = parallelizable
        
        def __call__(self, ts):
            return self._transform(ts)


InterpolationType = Literal["linear", "sinusoidal", "ease_in_out", "cubic"]


def linear_interpolation(t: float) -> float:
    """Linear interpolation: f(t) = t"""
    return t


def sinusoidal_interpolation(t: float) -> float:
    """Sinusoidal interpolation: f(t) = sin(t * π/2)"""
    return math.sin(t * math.pi / 2)


def ease_in_out_interpolation(t: float) -> float:
    """Ease in-out interpolation: f(t) = (1 - cos(t * π)) / 2"""
    return (1 - math.cos(t * math.pi)) / 2


def cubic_interpolation(t: float) -> float:
    """Cubic interpolation: f(t) = 3t² - 2t³"""
    return 3 * t * t - 2 * t * t * t


class SubframeInterpolator(TransformationBase):
    """
    MDAnalysis Transform class for creating interpolated subframes between existing trajectory frames.
    
    This transform allows you to specify a number of "subframes" that enables stepping between
    frames of the existing trajectory. The subframes are generated using various interpolation
    methods including linear and sinusoidal in-out.
    
    Parameters
    ----------
    subframes : int
        Number of subframes to create between each pair of existing frames.
        For example, subframes=3 will create 3 intermediate frames between each pair.
    interpolation : str, optional
        Type of interpolation to use. Options: "linear", "sinusoidal", "ease_in_out", "cubic".
        Default is "linear".
    ag : MDAnalysis.AtomGroup, optional
        Atom group to apply transformation to. If None, applies to all atoms in timestep.
    periodic_correction : bool, optional
        Whether to apply periodic boundary condition corrections during interpolation.
        Default is True.
    max_threads : int, optional
        Maximum number of threads for parallel processing. Default is 1.
    parallelizable : bool, optional
        Whether this transformation can be parallelized. Default is True.
        
    Attributes
    ----------
    current_frame : int
        The current logical frame number (including subframes)
    next_frame_positions : np.ndarray or None
        Cached positions for the next frame for interpolation
    interpolation_func : Callable
        The interpolation function being used
        
    Examples
    --------
    >>> import MDAnalysis as mda
    >>> from molecularnodes.entities.trajectory.transforms import SubframeInterpolator
    >>> 
    >>> # Create universe
    >>> u = mda.Universe("topology.pdb", "trajectory.xtc")
    >>> 
    >>> # Create transform with 4 subframes and sinusoidal interpolation
    >>> transform = SubframeInterpolator(subframes=4, interpolation="sinusoidal")
    >>> 
    >>> # Add to trajectory
    >>> u.trajectory.add_transformations(transform)
    >>> 
    >>> # Now trajectory will have interpolated frames
    >>> for ts in u.trajectory:
    >>>     print(f"Frame {ts.frame}: {ts.positions[:5]}")  # First 5 atom positions
    """
    
    def __init__(
        self,
        subframes: int,
        interpolation: InterpolationType = "linear",
        ag=None,
        periodic_correction: bool = True,
        max_threads: int = 1,
        parallelizable: bool = True,
    ):
        super().__init__(max_threads=max_threads, parallelizable=parallelizable)
        
        if subframes < 0:
            raise ValueError("subframes must be non-negative")
        
        self.subframes = subframes
        self.ag = ag
        self.periodic_correction = periodic_correction
        
        # Set up interpolation function
        self.interpolation_type = interpolation
        self.interpolation_func = self._get_interpolation_func(interpolation)
        
        # State variables
        self.current_frame = 0
        self.next_frame_positions = None
        self._last_original_frame = -1
        self._cached_positions = {}
        
    def _get_interpolation_func(self, interpolation: InterpolationType) -> Callable[[float], float]:
        """Get the interpolation function for the specified type."""
        interpolation_map = {
            "linear": linear_interpolation,
            "sinusoidal": sinusoidal_interpolation,
            "ease_in_out": ease_in_out_interpolation,
            "cubic": cubic_interpolation,
        }
        
        if interpolation not in interpolation_map:
            raise ValueError(
                f"Unknown interpolation type: {interpolation}. "
                f"Available options: {list(interpolation_map.keys())}"
            )
            
        return interpolation_map[interpolation]
    
    def _get_frame_info(self, frame_idx: int) -> tuple[int, int, float]:
        """
        Get original frame indices and interpolation factor for a given subframe index.
        
        Parameters
        ----------
        frame_idx : int
            The current frame index (including subframes)
            
        Returns
        -------
        tuple[int, int, float]
            (original_frame1, original_frame2, interpolation_factor)
        """
        if self.subframes == 0:
            return frame_idx, frame_idx, 0.0
            
        # Calculate which original frames we're interpolating between
        frame1 = frame_idx // (self.subframes + 1)
        subframe_idx = frame_idx % (self.subframes + 1)
        frame2 = frame1 + 1
        
        # Calculate interpolation factor (0.0 to 1.0)
        if subframe_idx == 0:
            t = 0.0
        else:
            t = subframe_idx / (self.subframes + 1)
            
        return frame1, frame2, t
    
    def _apply_periodic_correction(
        self, 
        pos1: npt.NDArray[np.float64], 
        pos2: npt.NDArray[np.float64], 
        dimensions: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """
        Apply periodic boundary condition corrections to positions.
        
        Parameters
        ----------
        pos1 : np.ndarray
            Reference positions
        pos2 : np.ndarray  
            Positions to correct
        dimensions : np.ndarray
            Box dimensions
            
        Returns
        -------
        np.ndarray
            Corrected positions
        """
        if dimensions is None or len(dimensions) < 3:
            return pos2
            
        corrected_pos = pos2.copy()
        
        for i in range(3):
            diff = pos2[:, i] - pos1[:, i]
            half_box = dimensions[i] / 2
            
            # Apply minimum image convention
            corrected_pos[:, i] -= dimensions[i] * np.round(diff / dimensions[i])
            
        return corrected_pos
    
    def _interpolate_positions(
        self,
        pos1: npt.NDArray[np.float64],
        pos2: npt.NDArray[np.float64], 
        t: float,
        dimensions: npt.NDArray[np.float64] = None
    ) -> npt.NDArray[np.float64]:
        """
        Interpolate between two sets of positions.
        
        Parameters
        ----------
        pos1 : np.ndarray
            Starting positions
        pos2 : np.ndarray
            Ending positions  
        t : float
            Interpolation factor (0.0 to 1.0)
        dimensions : np.ndarray, optional
            Box dimensions for periodic correction
            
        Returns
        -------
        np.ndarray
            Interpolated positions
        """
        if t == 0.0:
            return pos1
        elif t == 1.0:
            return pos2
            
        # Apply periodic correction if requested
        if self.periodic_correction and dimensions is not None:
            pos2 = self._apply_periodic_correction(pos1, pos2, dimensions)
        
        # Apply interpolation function to get smooth transition
        smooth_t = self.interpolation_func(t)
        
        # Linear interpolation with smooth factor
        return pos1 + smooth_t * (pos2 - pos1)
    
    def _transform(self, ts):
        """
        Apply the subframe interpolation transformation to the timestep.
        
        Parameters
        ----------
        ts : MDAnalysis.coordinates.base.Timestep
            The current timestep
            
        Returns
        -------
        MDAnalysis.coordinates.base.Timestep
            The modified timestep with interpolated positions
        """
        # Get frame information
        frame1_idx, frame2_idx, t = self._get_frame_info(ts.frame)
        
        # If no interpolation needed (t=0 or subframes=0), return as-is
        if t == 0.0 or self.subframes == 0:
            return ts
            
        # Get trajectory reference
        trajectory = ts.trajectory
        max_frame = len(trajectory) - 1
        
        # Clamp frame indices to valid range
        frame1_idx = min(frame1_idx, max_frame)
        frame2_idx = min(frame2_idx, max_frame)
        
        # If we're at the last frame, no interpolation possible
        if frame1_idx == frame2_idx:
            return ts
            
        # Get positions for both frames
        # Store current frame to restore later
        current_frame = trajectory.frame
        
        # Get positions from frame1 (should be current positions)
        trajectory[frame1_idx]
        pos1 = trajectory.ts.positions.copy()
        if self.ag is not None:
            pos1 = pos1[self.ag.ix]
            
        # Get positions from frame2  
        trajectory[frame2_idx]
        pos2 = trajectory.ts.positions.copy()
        if self.ag is not None:
            pos2 = pos2[self.ag.ix]
            
        # Restore original frame
        trajectory[current_frame]
        
        # Get box dimensions for periodic correction
        dimensions = getattr(ts, 'dimensions', None)
        
        # Interpolate positions
        interpolated_pos = self._interpolate_positions(pos1, pos2, t, dimensions)
        
        # Apply interpolated positions to timestep
        if self.ag is not None:
            ts.positions[self.ag.ix] = interpolated_pos
        else:
            ts.positions = interpolated_pos
            
        return ts
    
    def __repr__(self):
        return (
            f"SubframeInterpolator(subframes={self.subframes}, "
            f"interpolation='{self.interpolation_type}', "
            f"periodic_correction={self.periodic_correction})"
        )


class AdaptiveSubframeInterpolator(SubframeInterpolator):
    """
    Advanced subframe interpolator that adaptively adjusts interpolation based on motion.
    
    This class extends SubframeInterpolator to provide adaptive interpolation that
    can adjust the number of subframes based on the amount of motion between frames.
    
    Parameters
    ----------
    max_subframes : int
        Maximum number of subframes to create
    motion_threshold : float, optional
        Threshold for motion detection (Angstroms). Default is 2.0.
    min_subframes : int, optional
        Minimum number of subframes to create. Default is 0.
    **kwargs
        Additional arguments passed to SubframeInterpolator
        
    Examples
    --------
    >>> transform = AdaptiveSubframeInterpolator(
    ...     max_subframes=8, 
    ...     motion_threshold=1.5,
    ...     interpolation="ease_in_out"
    ... )
    """
    
    def __init__(
        self,
        max_subframes: int,
        motion_threshold: float = 2.0,
        min_subframes: int = 0,
        **kwargs
    ):
        # Initialize with max_subframes initially
        super().__init__(subframes=max_subframes, **kwargs)
        
        self.max_subframes = max_subframes
        self.min_subframes = min_subframes
        self.motion_threshold = motion_threshold
        self._motion_cache = {}
        
    def _calculate_motion(self, pos1: npt.NDArray[np.float64], pos2: npt.NDArray[np.float64]) -> float:
        """Calculate the average motion between two position sets."""
        if self.periodic_correction and hasattr(self, '_last_dimensions'):
            pos2 = self._apply_periodic_correction(pos1, pos2, self._last_dimensions)
        
        distances = np.linalg.norm(pos2 - pos1, axis=1)
        return np.mean(distances)
    
    def _adapt_subframes(self, motion: float) -> int:
        """Adapt the number of subframes based on motion."""
        if motion <= self.motion_threshold:
            return self.min_subframes
        
        # Scale subframes based on motion
        scale_factor = min(motion / self.motion_threshold, 3.0)  # Cap at 3x threshold
        adapted_subframes = int(self.min_subframes + 
                              (self.max_subframes - self.min_subframes) * 
                              (scale_factor - 1) / 2)
        
        return min(adapted_subframes, self.max_subframes)
    
    def _transform(self, ts):
        """Apply adaptive subframe interpolation."""
        frame1_idx, frame2_idx, t = self._get_frame_info(ts.frame)
        
        # Calculate motion if we haven't cached it
        cache_key = (frame1_idx, frame2_idx)
        if cache_key not in self._motion_cache and frame1_idx != frame2_idx:
            trajectory = ts.trajectory
            current_frame = trajectory.frame
            
            trajectory[frame1_idx]
            pos1 = trajectory.ts.positions.copy()
            trajectory[frame2_idx] 
            pos2 = trajectory.ts.positions.copy()
            trajectory[current_frame]
            
            if self.ag is not None:
                pos1 = pos1[self.ag.ix]
                pos2 = pos2[self.ag.ix]
                
            motion = self._calculate_motion(pos1, pos2)
            self._motion_cache[cache_key] = motion
            
            # Adapt subframes based on motion
            self.subframes = self._adapt_subframes(motion)
        
        # Store dimensions for periodic correction
        self._last_dimensions = getattr(ts, 'dimensions', None)
        
        return super()._transform(ts)