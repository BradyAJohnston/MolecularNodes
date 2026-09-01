import numpy as np
import pytest
from molecularnodes import framing


def random_basis(rng) -> np.ndarray:
    "A random orthonormal (right, up, forward) basis."
    basis, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(basis) < 0:
        basis[0] *= -1
    return basis


def frustum_slack(points, location, basis, frame, margin=0.0):
    """How much room is left between the points and each side of the frame.

    Zero means the points touch that side; negative means they spill out of it.
    """
    left, right, bottom, top = framing.apply_margin(frame, margin)
    offsets = np.asarray(points, dtype=np.float64) - location
    depth = offsets @ basis[2]
    x = offsets @ basis[0]
    y = offsets @ basis[1]
    return np.min(
        [
            np.min(right * depth - x),
            np.min(x - left * depth),
            np.min(top * depth - y),
            np.min(y - bottom * depth),
        ]
    ), depth


@pytest.mark.parametrize("shifted", [False, True])
def test_fit_is_exact_over_random_scenes(shifted):
    """Every point in frame, and at least one of them touching the edge.

    Anything less than touching means the camera pulled back further than the
    points required, which is what framing a box around them used to do.
    """
    rng = np.random.default_rng(0 if shifted else 1)
    for _ in range(200):
        points = rng.normal(size=(int(rng.integers(2, 400)), 3)) * rng.uniform(0.1, 50)
        points += rng.normal(size=3) * 20
        basis = random_basis(rng)
        half_x, half_y = rng.uniform(0.2, 1.5), rng.uniform(0.2, 1.5)
        shift = rng.uniform(-0.2, 0.2, size=2) if shifted else np.zeros(2)
        frame = (
            -half_x + shift[0],
            half_x + shift[0],
            -half_y + shift[1],
            half_y + shift[1],
        )

        location = framing.fit_camera_to_points(points, basis, frame)
        slack, depth = frustum_slack(points, location, basis, frame)

        assert (depth > 0).all(), "points ended up behind the camera"
        # the near-distance guard can legitimately hold the camera back off a
        # point set with no extent across the frame
        clamped = np.isclose(depth.min(), 1e-3)
        assert slack >= -1e-9, "points spilled out of frame"
        assert clamped or slack < 1e-6, "camera sat further back than needed"


def test_shifted_frustum_is_not_framed_as_if_centred():
    "A lens shift moves the frame, so the camera has to move with it."
    points = np.array([[-1.0, -1, 0], [1, 1, 0], [0, 0, 1]])
    basis = np.eye(3)
    basis[2] *= -1  # look down -Z
    centred = framing.fit_camera_to_points(points, basis, (-1, 1, -1, 1))
    shifted = framing.fit_camera_to_points(points, basis, (-0.6, 1.4, -1, 1))
    assert not np.allclose(centred, shifted)
    slack, _ = frustum_slack(points, shifted, basis, (-0.6, 1.4, -1, 1))
    assert slack >= -1e-9


def test_margin_pulls_back_and_negative_margin_crops():
    points = np.random.default_rng(2).normal(size=(50, 3))
    basis = random_basis(np.random.default_rng(3))
    frame = (-0.5, 0.5, -0.4, 0.4)
    centre = points.mean(axis=0)

    distances = [
        np.linalg.norm(framing.fit_camera_to_points(points, basis, frame, m) - centre)
        for m in (-0.2, 0.0, 0.2, 0.5)
    ]
    assert distances == sorted(distances), "more margin should mean further back"


def test_margin_leaves_the_frame_it_says_it_does():
    "A margin of 0.25 should leave a quarter of each half-frame empty."
    points = np.random.default_rng(4).normal(size=(100, 3))
    basis = np.array([[1.0, 0, 0], [0, 1, 0], [0, 0, -1]])
    frame = (-1.0, 1.0, -1.0, 1.0)

    location = framing.fit_camera_to_points(points, basis, frame, margin=0.25)
    offsets = points - location
    depth = offsets @ basis[2]
    # whichever axis ends up binding is the one pushed out to the margin
    used = max(
        np.max(np.abs(offsets @ basis[0]) / depth),
        np.max(np.abs(offsets @ basis[1]) / depth),
    )
    assert used == pytest.approx(0.75, rel=1e-6)


def test_single_and_identical_points_do_not_land_on_the_subject():
    for points in ([[1.0, 2.0, 3.0]], np.zeros((5, 3))):
        location = framing.fit_camera_to_points(
            points, np.array([[1.0, 0, 0], [0, 1, 0], [0, 0, -1]]), (-1, 1, -1, 1)
        )
        assert np.isfinite(location).all()
        assert np.linalg.norm(location - np.asarray(points)[0]) > 0


def test_framing_points_beats_framing_their_bounding_box():
    """The reason for solving on the points at all.

    An axis-aligned box projects larger than the points inside it from any
    direction that isn't straight down an axis.
    """
    rng = np.random.default_rng(5)
    points = rng.normal(size=(2000, 3)) @ np.diag([20.0, 3.0, 3.0])
    basis = random_basis(rng)
    frame = (-0.5, 0.5, -0.35, 0.35)
    low, high = points.min(axis=0), points.max(axis=0)
    corners = np.array(
        [
            [x, y, z]
            for x in (low[0], high[0])
            for y in (low[1], high[1])
            for z in (low[2], high[2])
        ]
    )

    centre = points.mean(axis=0)
    on_points = np.linalg.norm(framing.fit_camera_to_points(points, basis, frame) - centre)
    on_box = np.linalg.norm(framing.fit_camera_to_points(corners, basis, frame) - centre)
    assert on_points < on_box


@pytest.mark.parametrize(
    "bad,match",
    [
        ([(0, 0, 0), (1, 2)], "could not be read as numbers"),
        ([(1, 2)], r"shape \(N, 3\)"),
        ([], r"shape \(N, 3\)"),
        ([(0, 0, float("nan"))], "NaN or infinite"),
    ],
)
def test_as_points_rejects_what_it_cannot_frame(bad, match):
    with pytest.raises(ValueError, match=match):
        framing.as_points(bad)


def test_fit_rejects_a_bad_basis_and_an_impossible_margin():
    points = [[0.0, 0, 0], [1, 1, 1]]
    with pytest.raises(ValueError, match=r"\(3, 3\) matrix"):
        framing.fit_camera_to_points(points, np.eye(2), (-1, 1, -1, 1))
    with pytest.raises(ValueError, match="margin must be greater than -1"):
        framing.fit_camera_to_points(points, np.eye(3), (-1, 1, -1, 1), margin=-1.0)
    with pytest.raises(ValueError, match="degenerate"):
        framing.fit_camera_to_points(points, np.eye(3), (-1, 1, -1, 1), margin=1.0)


def test_enclosing_sphere_encloses():
    rng = np.random.default_rng(6)
    for shape in ([1.0, 1.0, 1.0], [30.0, 1.0, 1.0], [5.0, 5.0, 0.01]):
        points = rng.normal(size=(500, 3)) @ np.diag(shape)
        centre, radius = framing.enclosing_sphere(points)
        assert (np.linalg.norm(points - centre, axis=1) <= radius + 1e-9).all()
        # half the widest pair of points is a lower bound on the smallest
        # enclosing sphere, so this bounds how loose the approximation can be
        widest = np.max(np.linalg.norm(points[:, None] - points[None], axis=-1))
        assert radius <= 1.15 * widest / 2


def test_enclosing_sphere_of_one_point():
    centre, radius = framing.enclosing_sphere([[1.0, 2.0, 3.0]])
    assert radius == 0.0
    assert np.allclose(centre, [1.0, 2.0, 3.0])


def test_orthographic_fit_covers_the_points():
    rng = np.random.default_rng(7)
    points = rng.normal(size=(200, 3)) @ np.diag([4.0, 2.0, 1.0])
    basis = random_basis(rng)
    frame = (-1.0, 1.0, -0.75, 0.75)

    location, factor = framing.fit_orthographic_to_points(points, basis, frame)
    offsets = points - location
    assert (offsets @ basis[2] > 0).all(), "points ended up behind the camera"

    left, right, bottom, top = [bound * factor for bound in frame]
    x, y = offsets @ basis[0], offsets @ basis[1]
    assert (x >= left - 1e-9).all() and (x <= right + 1e-9).all()
    assert (y >= bottom - 1e-9).all() and (y <= top + 1e-9).all()
