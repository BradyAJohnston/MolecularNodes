"""
Fitting a camera to a set of points.

The camera is placed so that every point falls inside the view frustum, as close
to the subject as that allows. This is solved directly on the points rather than
by building a box around them and asking Blender to frame that: an axis-aligned
box projects larger than the points it contains, and by a different amount from
each direction, so a box-framed camera pulls back further than it needs to and
by an inconsistent amount as the viewpoint changes.

The solve is exact and runs in a single pass over the points. For a fixed camera
orientation, each side of the frustum gives a linear constraint on the camera
position, and the closest camera satisfying all four is read straight off the
extremes of the point set - see :func:`fit_camera_to_points`.

Notes
-----
Choosing *where* to look from is a separate problem from how close to be, and a
well studied one - Vazquez et al.'s viewpoint entropy [1]_ and Secord et al.'s
perceptual model of viewpoint preference [2]_ both score candidate directions.
Nothing here picks a direction: the orientation is given, and only the position
is solved for.

An alternative to the exact solve is to fit the smallest enclosing sphere of the
points [3]_ and frame that. It gives identical framing from every direction,
which is useful for turntables, at the cost of wasting the frame on anything not
spherical. :func:`enclosing_sphere` provides it for that purpose.

References
----------
.. [1] Vazquez, Feixas, Sbert and Heidrich (2003). "Automatic View Selection
   Using Viewpoint Entropy and its Application to Image-Based Modelling".
   Computer Graphics Forum 22(4).
.. [2] Secord, Lu, Finkelstein, Singh and Nealen (2011). "Perceptual Models of
   Viewpoint Preference". ACM Transactions on Graphics 30(5).
.. [3] Welzl (1991). "Smallest enclosing disks (balls and ellipsoids)".
   New Results and New Trends in Computer Science, LNCS 555.
"""

import numpy as np
import numpy.typing as npt

# a camera solved onto a single point, or onto points with no extent across the
# frame, would otherwise land exactly on top of them
_MIN_DISTANCE = 1e-3

Points = npt.NDArray[np.float64]


def as_points(value) -> Points:
    """
    Coerce a target into an ``(N, 3)`` array of positions.

    Parameters
    ----------
    value : array_like
        Anything shaped like a sequence of 3D positions - a list of tuples, an
        ``(N, 3)`` array, or a single ``(3,)`` position.

    Returns
    -------
    ndarray
        An ``(N, 3)`` array of float positions.

    Raises
    ------
    ValueError
        If the values are not 3D positions, or there are none of them.
    """
    try:
        points = np.asarray(value, dtype=np.float64)
    except (ValueError, TypeError) as e:
        # a ragged sequence, or one holding something that isn't a number
        raise ValueError(
            "Expected a sequence of 3D positions with shape (N, 3), but the "
            f"values given could not be read as numbers: {e}"
        ) from e
    if points.ndim == 1:
        points = points.reshape(1, -1)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(
            "Expected a sequence of 3D positions with shape (N, 3), got an array "
            f"of shape {tuple(np.shape(value))}."
        )
    if len(points) == 0:
        raise ValueError("Cannot frame an empty set of points.")
    if not np.isfinite(points).all():
        raise ValueError("Positions to frame contain NaN or infinite values.")
    return points


def apply_margin(
    frame: tuple[float, float, float, float], margin: float
) -> tuple[float, float, float, float]:
    """
    Shrink a camera frame towards its centre to leave a margin around the subject.

    Parameters
    ----------
    frame : tuple[float, float, float, float]
        ``(left, right, bottom, top)`` frustum bounds, as ratios of x (or y) to
        depth, so that a point at depth ``d`` is in frame when
        ``left * d <= x <= right * d``.
    margin : float
        Fraction of the frame to leave empty on each side. ``0.1`` leaves a ten
        percent border; a negative value crops in past the subject's edges.

    Returns
    -------
    tuple[float, float, float, float]
        The bounds scaled about the frame centre.
    """
    if margin <= -1.0:
        raise ValueError(f"margin must be greater than -1, got {margin}.")
    left, right, bottom, top = frame
    mid_x, mid_y = (left + right) / 2, (bottom + top) / 2
    half_x = (right - left) / 2 * (1.0 - margin)
    half_y = (top - bottom) / 2 * (1.0 - margin)
    return (mid_x - half_x, mid_x + half_x, mid_y - half_y, mid_y + half_y)


def fit_camera_to_points(
    points: npt.ArrayLike,
    basis: npt.ArrayLike,
    frame: tuple[float, float, float, float],
    margin: float = 0.0,
) -> Points:
    """
    Find the closest camera position that keeps every point in frame.

    Parameters
    ----------
    points : array_like
        ``(N, 3)`` world-space positions to fit into the frame.
    basis : array_like
        ``(3, 3)`` orthonormal camera basis as rows ``(right, up, forward)``,
        where ``forward`` is the direction the camera looks along.
    frame : tuple[float, float, float, float]
        ``(left, right, bottom, top)`` frustum bounds as ratios to depth. Need
        not be symmetric, so a shifted (off-axis) frustum is handled.
    margin : float, default 0.0
        Fraction of the frame to leave empty around the subject.

    Returns
    -------
    ndarray
        The ``(3,)`` world-space camera position. Orientation is untouched -
        only where the camera sits along its own axes is solved for.

    Notes
    -----
    With the orientation fixed, a point ``p`` is inside the right-hand side of
    the frustum when ``x_p <= right * depth_p``. Writing the camera position in
    its own basis as ``(c_r, c_u, c_f)`` and the point's as ``(a, b, z)``, that
    is ``a - c_r <= right * (z - c_f)``, or

        ``a - right * z <= c_r - right * c_f``

    The left-hand side depends only on the points, so one pass of ``max`` over
    them gives the binding constraint for that side of the frustum. Four sides
    give four such constraints; the two horizontal ones fix how far back the
    camera must be for width, the two vertical ones for height, and the further
    of the two wins. The lateral position is then the midpoint of whatever
    freedom is left, which centres the subject in frame.

    Being linear in the point coordinates, this is exact and needs no search -
    unlike fitting a box and framing that, it never pulls back further than the
    points require.
    """
    points = as_points(points)
    basis = np.asarray(basis, dtype=np.float64)
    if basis.shape != (3, 3):
        raise ValueError(f"basis must be a (3, 3) matrix, got {basis.shape}.")
    right_axis, up_axis, forward_axis = basis
    left, right, bottom, top = apply_margin(frame, margin)
    if right <= left or top <= bottom:
        raise ValueError(f"Camera frame {frame} with margin {margin} is degenerate.")

    # the point set in the camera's own basis
    a = points @ right_axis
    b = points @ up_axis
    z = points @ forward_axis

    # the binding point for each side of the frustum
    bound_right = np.max(a - right * z)
    bound_left = np.max(left * z - a)
    bound_top = np.max(b - top * z)
    bound_bottom = np.max(bottom * z - b)

    # how far along the view axis the camera has to sit to satisfy each pair of
    # opposing sides; the smaller value is further back, so it is the binding one
    depth_x = (bound_right + bound_left) / (left - right)
    depth_y = (bound_top + bound_bottom) / (bottom - top)
    depth = min(depth_x, depth_y)

    # a point set with no extent across the frame solves to the point itself
    nearest = float(np.min(z))
    depth = min(depth, nearest - _MIN_DISTANCE)

    # centre the subject in whatever freedom the binding axis left over
    offset_r = (bound_right - bound_left + (left + right) * depth) / 2
    offset_u = (bound_top - bound_bottom + (bottom + top) * depth) / 2
    return offset_r * right_axis + offset_u * up_axis + depth * forward_axis


def fit_orthographic_to_points(
    points: npt.ArrayLike,
    basis: npt.ArrayLike,
    frame: tuple[float, float, float, float],
    margin: float = 0.0,
) -> tuple[Points, float]:
    """
    Fit an orthographic camera to a set of points.

    An orthographic frustum is a box, so distance does not change the framing -
    the scale does. The camera is centred on the points and pulled back clear of
    them, and the factor its scale needs to change by is returned alongside.

    Parameters
    ----------
    points : array_like
        ``(N, 3)`` world-space positions to fit into the frame.
    basis : array_like
        ``(3, 3)`` orthonormal camera basis as rows ``(right, up, forward)``.
    frame : tuple[float, float, float, float]
        ``(left, right, bottom, top)`` bounds of the current frame in world
        units, as read off the camera at its current scale.
    margin : float, default 0.0
        Fraction of the frame to leave empty around the subject.

    Returns
    -------
    location : ndarray
        The ``(3,)`` world-space camera position.
    scale_factor : float
        What to multiply the camera's ``ortho_scale`` by.
    """
    points = as_points(points)
    basis = np.asarray(basis, dtype=np.float64)
    right_axis, up_axis, forward_axis = basis
    left, right, bottom, top = apply_margin(frame, margin)
    if right <= left or top <= bottom:
        raise ValueError(f"Camera frame {frame} with margin {margin} is degenerate.")

    a = points @ right_axis
    b = points @ up_axis
    z = points @ forward_axis

    # how much wider the frame has to be to hold the points, on the tighter axis
    span_x = float(np.ptp(a)) or _MIN_DISTANCE
    span_y = float(np.ptp(b)) or _MIN_DISTANCE
    scale_factor = max(span_x / (right - left), span_y / (top - bottom))

    # sit clear of the points, since nothing in front of the camera is drawn
    depth = float(np.min(z)) - max(float(np.ptp(z)), _MIN_DISTANCE)
    centre_r = (float(np.max(a)) + float(np.min(a))) / 2
    centre_u = (float(np.max(b)) + float(np.min(b))) / 2
    return (
        centre_r * right_axis + centre_u * up_axis + depth * forward_axis,
        scale_factor,
    )


def enclosing_sphere(points: npt.ArrayLike) -> tuple[Points, float]:
    """
    A sphere containing every point, for framing that does not change with angle.

    Ritter's two-pass approximation, refined until every point is enclosed. It
    comes out within a few percent of the true smallest enclosing sphere, which
    Welzl's algorithm solves exactly in expected linear time [3]_ at rather more
    complexity than framing a camera warrants.

    Framing this sphere instead of the points themselves gives the same framing
    from every direction - worth having for a turntable, where solving each
    frame exactly makes an elongated subject appear to breathe as it rotates.

    Parameters
    ----------
    points : array_like
        ``(N, 3)`` positions to enclose.

    Returns
    -------
    centre : ndarray
        The ``(3,)`` centre of the sphere.
    radius : float
        Its radius.
    """
    points = as_points(points)
    if len(points) == 1:
        return points[0].copy(), 0.0

    # Ritter: start from a point pair found by walking the extremes, then grow
    first = points[np.argmax(np.sum((points - points[0]) ** 2, axis=1))]
    second = points[np.argmax(np.sum((points - first) ** 2, axis=1))]
    centre = (first + second) / 2
    radius = float(np.linalg.norm(second - first)) / 2

    for _ in range(2):
        offsets = points - centre
        distances = np.linalg.norm(offsets, axis=1)
        outside = distances > radius
        if not outside.any():
            break
        # grow just enough to take in the furthest stray point, moving the
        # centre half way towards it so the far side stays covered
        furthest = int(np.argmax(distances))
        grown = (radius + distances[furthest]) / 2
        centre = centre + offsets[furthest] * (1 - radius / distances[furthest]) / 2
        radius = float(grown)

    # a final sweep guarantees enclosure however the approximation landed
    radius = float(np.max(np.linalg.norm(points - centre, axis=1)))

    # Ritter's two passes can come out worse than simply centring on the mean,
    # which happens for point sets that are already close to spherical
    centroid = points.mean(axis=0)
    centroid_radius = float(np.max(np.linalg.norm(points - centroid, axis=1)))
    if centroid_radius < radius:
        return centroid, centroid_radius
    return centre, radius
