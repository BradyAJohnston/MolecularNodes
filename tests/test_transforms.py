"""
Tests for subframe interpolation transforms.
"""

import math
import numpy as np
import pytest
import MDAnalysis as mda
from molecularnodes.entities.trajectory.transforms import (
    SubframeInterpolator,
    AdaptiveSubframeInterpolator,
    linear_interpolation,
    sinusoidal_interpolation,
    ease_in_out_interpolation,
    cubic_interpolation,
)
from .constants import data_dir


class TestInterpolationFunctions:
    """Test the mathematical interpolation functions."""
    
    def test_linear_interpolation(self):
        """Test linear interpolation function."""
        assert linear_interpolation(0.0) == 0.0
        assert linear_interpolation(0.5) == 0.5
        assert linear_interpolation(1.0) == 1.0
        
    def test_sinusoidal_interpolation(self):
        """Test sinusoidal interpolation function."""
        assert sinusoidal_interpolation(0.0) == 0.0
        assert abs(sinusoidal_interpolation(1.0) - 1.0) < 1e-10
        # Check that it's smooth at 0.5
        mid_val = sinusoidal_interpolation(0.5)
        assert 0.6 < mid_val < 0.8  # sin(π/4) ≈ 0.707
        
    def test_ease_in_out_interpolation(self):
        """Test ease in-out interpolation function."""
        assert ease_in_out_interpolation(0.0) == 0.0
        assert abs(ease_in_out_interpolation(1.0) - 1.0) < 1e-10
        assert ease_in_out_interpolation(0.5) == 0.5
        
    def test_cubic_interpolation(self):
        """Test cubic interpolation function."""
        assert cubic_interpolation(0.0) == 0.0
        assert cubic_interpolation(1.0) == 1.0
        assert cubic_interpolation(0.5) == 0.5


class TestSubframeInterpolator:
    """Test the SubframeInterpolator class."""
    
    @pytest.fixture(scope="class")
    def universe(self):
        """Create a test universe."""
        top = data_dir / "md_ppr/box.gro"
        traj = data_dir / "md_ppr/first_5_frames.xtc"
        u = mda.Universe(top, traj)
        return u
    
    def test_init_valid_parameters(self):
        """Test initialization with valid parameters."""
        transform = SubframeInterpolator(subframes=3)
        assert transform.subframes == 3
        assert transform.interpolation_type == "linear"
        assert transform.periodic_correction is True
        
    def test_init_invalid_subframes(self):
        """Test initialization with invalid subframes."""
        with pytest.raises(ValueError, match="subframes must be non-negative"):
            SubframeInterpolator(subframes=-1)
            
    def test_init_invalid_interpolation(self):
        """Test initialization with invalid interpolation type."""
        with pytest.raises(ValueError, match="Unknown interpolation type"):
            SubframeInterpolator(subframes=1, interpolation="invalid")
            
    def test_all_interpolation_types(self):
        """Test that all interpolation types can be initialized."""
        for interp_type in ["linear", "sinusoidal", "ease_in_out", "cubic"]:
            transform = SubframeInterpolator(subframes=2, interpolation=interp_type)
            assert transform.interpolation_type == interp_type
            
    def test_get_frame_info_no_subframes(self):
        """Test frame info calculation with no subframes."""
        transform = SubframeInterpolator(subframes=0)
        frame1, frame2, t = transform._get_frame_info(5)
        assert frame1 == 5
        assert frame2 == 5
        assert t == 0.0
        
    def test_get_frame_info_with_subframes(self):
        """Test frame info calculation with subframes."""
        transform = SubframeInterpolator(subframes=3)  # 4 total frames per original frame
        
        # Test various frame indices
        test_cases = [
            (0, 0, 1, 0.0),      # Original frame 0
            (1, 0, 1, 0.25),     # First subframe
            (2, 0, 1, 0.5),      # Second subframe  
            (3, 0, 1, 0.75),     # Third subframe
            (4, 1, 2, 0.0),      # Original frame 1
            (5, 1, 2, 0.25),     # First subframe of frame 1
        ]
        
        for frame_idx, expected_f1, expected_f2, expected_t in test_cases:
            frame1, frame2, t = transform._get_frame_info(frame_idx)
            assert frame1 == expected_f1, f"Frame {frame_idx}: expected frame1={expected_f1}, got {frame1}"
            assert frame2 == expected_f2, f"Frame {frame_idx}: expected frame2={expected_f2}, got {frame2}"
            assert abs(t - expected_t) < 1e-10, f"Frame {frame_idx}: expected t={expected_t}, got {t}"
            
    def test_periodic_correction(self):
        """Test periodic boundary condition correction."""
        transform = SubframeInterpolator(subframes=1, periodic_correction=True)
        
        # Create test positions that cross periodic boundary
        pos1 = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        pos2 = np.array([[9.0, 0.0, 0.0], [1.0, 1.0, 1.0]])  # Wrapped around boundary
        dimensions = np.array([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
        
        corrected = transform._apply_periodic_correction(pos1, pos2, dimensions)
        
        # The 9.0 should be corrected to -1.0 (closer to 0.0)
        expected = np.array([[-1.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        np.testing.assert_allclose(corrected, expected, atol=1e-10)
        
    def test_interpolate_positions_linear(self):
        """Test position interpolation with linear method."""
        transform = SubframeInterpolator(subframes=1, interpolation="linear")
        
        pos1 = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        pos2 = np.array([[2.0, 0.0, 0.0], [3.0, 1.0, 1.0]])
        
        # Test t=0.0 (should return pos1)
        result = transform._interpolate_positions(pos1, pos2, 0.0)
        np.testing.assert_allclose(result, pos1)
        
        # Test t=1.0 (should return pos2)
        result = transform._interpolate_positions(pos1, pos2, 1.0)
        np.testing.assert_allclose(result, pos2)
        
        # Test t=0.5 (should return midpoint)
        result = transform._interpolate_positions(pos1, pos2, 0.5)
        expected = np.array([[1.0, 0.0, 0.0], [2.0, 1.0, 1.0]])
        np.testing.assert_allclose(result, expected)
        
    def test_interpolate_positions_nonlinear(self):
        """Test position interpolation with non-linear methods."""
        for interp_type in ["sinusoidal", "ease_in_out", "cubic"]:
            transform = SubframeInterpolator(subframes=1, interpolation=interp_type)
            
            pos1 = np.array([[0.0, 0.0, 0.0]])
            pos2 = np.array([[1.0, 0.0, 0.0]])
            
            # Test endpoints
            result = transform._interpolate_positions(pos1, pos2, 0.0)
            np.testing.assert_allclose(result, pos1)
            
            result = transform._interpolate_positions(pos1, pos2, 1.0)
            np.testing.assert_allclose(result, pos2)
            
            # Test midpoint - should be different from linear for non-linear methods
            result = transform._interpolate_positions(pos1, pos2, 0.5)
            if interp_type == "ease_in_out":
                # Ease in-out should give exactly 0.5 at t=0.5
                expected = np.array([[0.5, 0.0, 0.0]])
                np.testing.assert_allclose(result, expected)
            else:
                # Other methods should be different from linear
                linear_midpoint = np.array([[0.5, 0.0, 0.0]])
                assert not np.allclose(result, linear_midpoint)
                
    def test_transform_no_subframes(self, universe):
        """Test transform with no subframes (should be pass-through)."""
        transform = SubframeInterpolator(subframes=0)
        
        # Get original positions
        universe.trajectory[0]
        original_positions = universe.trajectory.ts.positions.copy()
        
        # Apply transform
        ts = transform._transform(universe.trajectory.ts)
        
        # Should be unchanged
        np.testing.assert_allclose(ts.positions, original_positions)
        
    def test_transform_with_subframes(self, universe):
        """Test transform with subframes."""
        transform = SubframeInterpolator(subframes=2, interpolation="linear")
        
        # Get positions for frame 0 and 1
        universe.trajectory[0]
        pos0 = universe.trajectory.ts.positions.copy()
        universe.trajectory[1] 
        pos1 = universe.trajectory.ts.positions.copy()
        
        # Test interpolation at t=0.5 (should be midpoint)
        universe.trajectory[0]
        ts = universe.trajectory.ts
        
        # Manually set frame to simulate subframe 1 (t=1/3)
        ts.frame = 1  # This represents the first subframe
        
        # Apply transform
        transformed_ts = transform._transform(ts)
        
        # Should be interpolated between pos0 and pos1
        # Note: exact verification depends on the frame calculation logic
        assert not np.allclose(transformed_ts.positions, pos0)
        assert not np.allclose(transformed_ts.positions, pos1)
        
    def test_transform_with_atom_group(self, universe):
        """Test transform applied to specific atom group."""
        # Select first 10 atoms
        ag = universe.atoms[:10]
        transform = SubframeInterpolator(subframes=1, ag=ag)
        
        universe.trajectory[0]
        original_positions = universe.trajectory.ts.positions.copy()
        
        # Apply transform (this is a simplified test)
        ts = transform._transform(universe.trajectory.ts)
        
        # Only the selected atoms should potentially be modified
        # The rest should remain unchanged
        # Note: actual verification would require proper subframe setup
        assert ts.positions.shape == original_positions.shape
        
    def test_repr(self):
        """Test string representation."""
        transform = SubframeInterpolator(
            subframes=3, 
            interpolation="sinusoidal",
            periodic_correction=False
        )
        
        repr_str = repr(transform)
        assert "SubframeInterpolator" in repr_str
        assert "subframes=3" in repr_str
        assert "sinusoidal" in repr_str
        assert "periodic_correction=False" in repr_str


class TestAdaptiveSubframeInterpolator:
    """Test the AdaptiveSubframeInterpolator class."""
    
    @pytest.fixture(scope="class")
    def universe(self):
        """Create a test universe."""
        top = data_dir / "md_ppr/box.gro"
        traj = data_dir / "md_ppr/first_5_frames.xtc"
        u = mda.Universe(top, traj)
        return u
    
    def test_init(self):
        """Test initialization."""
        transform = AdaptiveSubframeInterpolator(
            max_subframes=8,
            motion_threshold=1.5,
            min_subframes=1
        )
        
        assert transform.max_subframes == 8
        assert transform.motion_threshold == 1.5
        assert transform.min_subframes == 1
        assert transform.subframes == 8  # Should start with max
        
    def test_calculate_motion(self):
        """Test motion calculation."""
        transform = AdaptiveSubframeInterpolator(max_subframes=4)
        
        # Create test positions
        pos1 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        pos2 = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        
        motion = transform._calculate_motion(pos1, pos2)
        expected_motion = 1.0  # Average distance is 1.0 Angstrom
        
        assert abs(motion - expected_motion) < 1e-10
        
    def test_adapt_subframes(self):
        """Test subframe adaptation logic."""
        transform = AdaptiveSubframeInterpolator(
            max_subframes=8,
            min_subframes=2,
            motion_threshold=2.0
        )
        
        # Low motion should give minimum subframes
        adapted = transform._adapt_subframes(1.0)  # Below threshold
        assert adapted == 2
        
        # Motion at threshold should give minimum subframes
        adapted = transform._adapt_subframes(2.0)  # At threshold
        assert adapted == 2
        
        # High motion should give more subframes
        adapted = transform._adapt_subframes(4.0)  # 2x threshold
        assert adapted > 2
        assert adapted <= 8
        
        # Very high motion should be capped at maximum
        adapted = transform._adapt_subframes(10.0)  # Very high
        assert adapted <= 8


class TestIntegrationWithMDAnalysis:
    """Integration tests with actual MDAnalysis trajectories."""
    
    @pytest.fixture(scope="class") 
    def universe(self):
        """Create a test universe."""
        top = data_dir / "md_ppr/box.gro"
        traj = data_dir / "md_ppr/first_5_frames.xtc"
        u = mda.Universe(top, traj)
        return u
    
    def test_add_transformation_to_trajectory(self, universe):
        """Test adding the transform to a trajectory."""
        transform = SubframeInterpolator(subframes=2)
        
        # Add transformation to trajectory
        universe.trajectory.add_transformations(transform)
        
        # Check that transformation was added
        assert len(universe.trajectory.transformations) > 0
        
    def test_trajectory_iteration_with_transform(self, universe):
        """Test iterating through trajectory with transform applied."""
        transform = SubframeInterpolator(subframes=1, interpolation="linear")
        universe.trajectory.add_transformations(transform)
        
        positions_list = []
        frame_count = 0
        
        # Iterate through first few frames
        for ts in universe.trajectory[:3]:
            positions_list.append(ts.positions.copy())
            frame_count += 1
            if frame_count >= 3:  # Limit to avoid long test
                break
                
        # Should have collected some positions
        assert len(positions_list) > 0
        
        # Positions should be different between frames
        if len(positions_list) > 1:
            assert not np.allclose(positions_list[0], positions_list[1])


@pytest.mark.parametrize("interpolation_type", ["linear", "sinusoidal", "ease_in_out", "cubic"])
@pytest.mark.parametrize("subframes", [0, 1, 3, 5])
def test_interpolation_consistency(interpolation_type, subframes):
    """Test that all interpolation types work consistently."""
    # Create simple test data
    transform = SubframeInterpolator(
        subframes=subframes,
        interpolation=interpolation_type
    )
    
    # Test basic functionality
    assert transform.interpolation_type == interpolation_type
    assert transform.subframes == subframes
    
    # Test interpolation function
    func = transform.interpolation_func
    assert func(0.0) == 0.0
    assert abs(func(1.0) - 1.0) < 1e-10
    
    # Test position interpolation
    pos1 = np.array([[0.0, 0.0, 0.0]])
    pos2 = np.array([[1.0, 1.0, 1.0]])
    
    result = transform._interpolate_positions(pos1, pos2, 0.5)
    assert result.shape == pos1.shape
    
    # Result should be between pos1 and pos2
    assert np.all(result >= np.minimum(pos1, pos2))
    assert np.all(result <= np.maximum(pos1, pos2))