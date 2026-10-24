import numpy as np
import pytest
from molecularnodes import viewpoint


def _cluster(spacing=0.1):
    grid = np.arange(-1, 2) * spacing
    return np.array(np.meshgrid(grid, grid, grid)).reshape(3, -1).T


def _shell_with_hole(hole=(-1.0, 0.0, 0.0), radius=1.2, n=4000, cap_degrees=30):
    points = viewpoint.fibonacci_sphere(n) * radius
    keep = points @ np.asarray(hole) / radius < np.cos(np.radians(cap_degrees))
    return points[keep]


def test_fibonacci_sphere_is_unit_and_balanced():
    dirs = viewpoint.fibonacci_sphere(500)
    assert np.allclose(np.linalg.norm(dirs, axis=1), 1.0)
    assert np.allclose(dirs.mean(axis=0), 0.0, atol=0.01)


@pytest.mark.parametrize(
    "forward", [(0, 1, 0), (1, 0, 0), (0.3, -0.5, 0.2), (0, 0, -1), (0, 0, 1)]
)
def test_look_along_is_a_proper_upright_rotation(forward):
    right, up, fwd = viewpoint.look_along(forward)
    assert np.allclose(fwd, np.asarray(forward) / np.linalg.norm(forward))
    # Blender's camera looks down its own -Z: columns (right, up, -forward)
    assert np.isclose(np.linalg.det(np.column_stack([right, up, -fwd])), 1.0)
    if abs(fwd[2]) < 0.95:
        assert up[2] > 0
        assert abs(right[2]) < 1e-9


def test_pca_view_looks_along_the_thinnest_axis():
    rng = np.random.default_rng(0)
    points = rng.normal(size=(2000, 3)) * [5.0, 0.5, 2.0]
    right, up, forward = viewpoint.pca_view(points)
    assert abs(forward @ [0, 1, 0]) > 0.99
    assert abs(right @ [1, 0, 0]) > 0.99
    assert up[2] > 0.99
    # deterministic sign: the camera looks in from -Y
    assert forward[1] > 0


@pytest.mark.parametrize("method", ["unobstructed", "visible_area"])
def test_search_looks_through_the_opening_in_a_shell(method):
    target = _cluster()
    shell = _shell_with_hole(hole=(-1, 0, 0))
    _, _, forward = viewpoint.best_view(
        target, method=method, occluders=shell, occluder_radii=0.1, resolution=128
    )
    # camera sits on the -X side, looking along +X through the hole
    assert forward @ [1, 0, 0] > np.cos(np.radians(30))


def test_unobstructed_scores_blocked_side_lower():
    target = _cluster()
    shell = _shell_with_hole(hole=(-1, 0, 0))
    score = viewpoint.unobstructed_scorer(target, shell, occluder_radii=0.1)
    through, blocked = score(np.array([[1.0, 0, 0], [-1.0, 0, 0]]))
    assert through > blocked + 5


def test_entropy_prefers_a_flat_plate_face_on():
    grid = np.arange(-5, 6) * 0.3
    plate = np.array(np.meshgrid(grid, grid, [0.0])).reshape(3, -1).T
    _, _, forward = viewpoint.best_view(plate, method="entropy", resolution=64)
    assert abs(forward[2]) > 0.95


def test_rasterize_single_sphere_area():
    winner = viewpoint.rasterize(
        np.zeros((1, 3)), np.array([1.0]), (0, 1, 0), (0, 0, 0), 2.0, 200
    )
    pixel = 4.0 / 200
    assert np.isclose((winner == 0).sum() * pixel**2, np.pi, rtol=0.02)


def test_rasterize_nearest_sphere_wins():
    spheres = np.array([[0.0, -1.0, 0.0], [0.0, 1.0, 0.0]])
    winner = viewpoint.rasterize(
        spheres, np.array([0.5, 0.5]), (0, 1, 0), (0, 0, 0), 1.0, 32
    )
    # looking along +Y, the sphere at y = -1 is nearer
    assert set(winner[winner >= 0]) == {0}


def test_thin_keeps_one_point_per_cell():
    points = np.random.default_rng(1).uniform(0, 1, size=(5000, 3))
    thinned = viewpoint.thin(points, 0.25)
    assert len(thinned) == 64


def test_ray_scorer_counts_clear_lines_of_sight():
    class Wall:
        """Blocks every ray heading towards +X."""

        def ray_cast(self, origin, direction):
            return (origin, None, None, None) if direction[0] > 0 else (None,) * 4

    score = viewpoint.ray_scorer([Wall()], _cluster())
    # looking along -X means rays go back towards +X and hit the wall
    blocked, clear = score(np.array([[-1.0, 0, 0], [1.0, 0, 0]]))
    # the integer part counts clear atoms; the fraction only breaks ties
    assert int(blocked) == 0
    assert int(clear) == 27


def test_unknown_method_raises():
    with pytest.raises(ValueError, match="Unknown method"):
        viewpoint.best_view(_cluster(), method="nope")


def test_without_near_drops_the_targets_own_geometry():
    target = _cluster()
    own = target + 0.02
    far = np.array([[5.0, 0, 0]])
    kept = viewpoint.without_near(np.concatenate([own, far]), target, 0.2)
    assert np.allclose(kept, far)
