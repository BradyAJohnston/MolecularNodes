"""
Choosing which direction to look at a subject from.

:mod:`molecularnodes.framing` solves how close the camera has to be for a given
orientation. This module picks the orientation. Every method returns a camera
basis as rows ``(right, up, forward)``, the convention ``framing`` uses, so the
two compose: pick a direction here, then fit the distance there.

Methods
-------
``"pca"``
    Look along the axis of least variance, with the longest axis across the
    frame. Closed form and instant, and blind to occlusion.
``"unobstructed"``
    The direction from which the most target atoms are directly visible: an atom
    counts when no nearer sphere covers its projected centre, whether that
    sphere belongs to the surroundings or to the target itself. A port of
    molar_vis's unobstructed view [1]_.
``"visible_area"``
    The same idea scored on projected area instead of atom count, using a small
    depth raster at a fixed world scale, so a view that shows the target broad
    side on beats one that shows the same atoms foreshortened.
``"entropy"``
    Viewpoint entropy [2]_: the Shannon entropy of how the visible area divides
    between labelled parts (atoms, residues, chains), the background included.
    Highest when many parts are visible at similar sizes. Needs no target, so it
    suits framing a whole structure.

The searched methods share :func:`search`: a Fibonacci sphere of candidate
directions, then rings of refinement around the best few, as in [1]_.

Radii and positions are in world units throughout, so 1 Angstrom is 0.1.

References
----------
.. [1] molar_vis, ``crates/molar_vis_core/src/unobstructed.rs``.
   https://github.com/yesint/molar_vis
.. [2] Vazquez, Feixas, Sbert and Heidrich (2003). "Automatic View Selection
   Using Viewpoint Entropy and its Application to Image-Based Modelling".
   Computer Graphics Forum 22(4).
"""

from typing import Callable, Sequence
import numpy as np
import numpy.typing as npt
from .framing import as_points, enclosing_sphere

Points = npt.NDArray[np.float64]
Scorer = Callable[[Points], npt.NDArray[np.float64]]

METHODS = ("pca", "unobstructed", "visible_area", "entropy")

# a carbon's van der Waals radius in world units, used when none are given
DEFAULT_RADIUS = 0.17

_WORLD_UP = np.array([0.0, 0.0, 1.0])


def fibonacci_sphere(n: int) -> Points:
    """
    ``n`` near-uniform unit directions on the sphere (a Fibonacci lattice).
    """
    i = np.arange(n) + 0.5
    z = 1.0 - 2.0 * i / n
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, None))
    theta = np.pi * (3.0 - np.sqrt(5.0)) * i
    return np.column_stack((r * np.cos(theta), r * np.sin(theta), z))


def tangent_basis(direction: npt.ArrayLike) -> tuple[Points, Points]:
    """
    Two unit vectors spanning the plane perpendicular to ``direction``.
    """
    d = _unit(direction)
    helper = np.array([1.0, 0.0, 0.0]) if abs(d[0]) < 0.9 else np.array([0, 1.0, 0])
    t1 = _unit(np.cross(helper, d))
    return t1, np.cross(d, t1)


def look_along(forward: npt.ArrayLike, up: npt.ArrayLike = _WORLD_UP) -> Points:
    """
    A camera basis looking along ``forward`` and kept upright against ``up``.

    Parameters
    ----------
    forward : array_like
        Direction the camera looks along.
    up : array_like, default (0, 0, 1)
        Which way is up. When ``forward`` is nearly parallel to it, world Y is
        used instead, so looking straight down still gives a defined roll.

    Returns
    -------
    ndarray
        ``(3, 3)`` rows ``(right, up, forward)``.
    """
    forward = _unit(forward)
    up = _unit(up)
    if abs(forward @ up) > 0.95:
        up = np.array([0.0, 1.0, 0.0])
    right = _unit(np.cross(forward, up))
    return np.array([right, np.cross(right, forward), forward])


def thin(points: npt.ArrayLike, spacing: float) -> Points:
    """
    Keep one point per cube of side ``spacing``.

    Rendered geometry can run to millions of vertices, far denser than the
    spheres standing in for it as occluders need to be.
    """
    points = as_points(points)
    if spacing <= 0:
        return points
    _, keep = np.unique(np.floor(points / spacing), axis=0, return_index=True)
    return points[np.sort(keep)]


def without_near(
    points: npt.ArrayLike, target: npt.ArrayLike, distance: float
) -> Points:
    """
    Drop the points within ``distance`` of any target point.

    Occluders read off rendered geometry include the target's own rendering,
    whose surface sits just in front of every atom centre and would hide all
    of them from every direction.
    """
    from scipy.spatial import cKDTree

    points = as_points(points)
    near = cKDTree(as_points(target)).query(points, distance_upper_bound=distance)[0]
    return points[~np.isfinite(near)]


def search(
    score: Scorer,
    resolution: int = 256,
    seeds: int = 3,
    rings: Sequence[float] = (10.0, 4.0, 1.5),
    samples: int = 16,
) -> Points:
    """
    The direction that maximises ``score``.

    Scores ``resolution`` directions spread over the sphere, keeps the best
    ``seeds`` that are at least 15 degrees apart, and refines each by scoring
    ``samples`` directions on rings of decreasing angle around it.

    Parameters
    ----------
    score : callable
        Maps an ``(M, 3)`` array of view directions (where the camera looks) to
        ``M`` scores, higher being better.
    resolution : int, default 256
        Directions in the coarse pass.
    seeds : int, default 3
        Separate maxima to refine; more guards against a narrow best view that
        the coarse pass only caught the edge of.
    rings : sequence of float, default (10, 4, 1.5)
        Refinement ring angles in degrees.
    samples : int, default 16
        Directions per ring.

    Returns
    -------
    ndarray
        The best ``(3,)`` unit view direction.
    """
    directions = fibonacci_sphere(max(int(resolution), 8))
    scores = np.asarray(score(directions), dtype=np.float64)
    separation = np.cos(np.radians(15.0))
    best: list[tuple[Points, float]] = []
    for i in np.argsort(-scores, kind="stable"):
        if all(directions[i] @ d < separation for d, _ in best):
            best.append((directions[i], float(scores[i])))
            if len(best) >= seeds:
                break

    angles = np.linspace(0.0, 2 * np.pi, samples, endpoint=False)
    for ring in rings:
        tilt = np.radians(ring)
        for k, (base, base_score) in enumerate(best):
            t1, t2 = tangent_basis(base)
            around = np.outer(np.cos(angles), t1) + np.outer(np.sin(angles), t2)
            candidates = base * np.cos(tilt) + around * np.sin(tilt)
            candidates /= np.linalg.norm(candidates, axis=1, keepdims=True)
            ring_scores = np.asarray(score(candidates), dtype=np.float64)
            i = int(np.argmax(ring_scores))
            if ring_scores[i] > base_score:
                best[k] = (candidates[i], float(ring_scores[i]))
    return max(best, key=lambda pair: pair[1])[0]


def pca_view(points: npt.ArrayLike) -> Points:
    """
    A basis looking along the axis of least variance, longest axis across.

    The sign of each axis is fixed by convention - up points up, and the camera
    looks in from the side facing -Y, as Blender's front view does - so the
    same structure gives the same view. Near-degenerate spreads (a sphere, a
    flat disc seen for its in-plane axes) have no well defined axes, and any
    rotation among them is as good as another; expect those to be unstable
    across platforms.
    """
    points = as_points(points)
    centred = points - points.mean(axis=0)
    _, vectors = np.linalg.eigh(centred.T @ centred)
    # eigh sorts ascending: least variance first
    forward, up, _ = vectors.T
    if forward[1] < 0:
        forward = -forward
    if up[2] < 0:
        up = -up
    up = _unit(up - (up @ forward) * forward)
    return np.array([np.cross(forward, up), up, forward])


def unobstructed_scorer(
    target: npt.ArrayLike,
    occluders: npt.ArrayLike | None = None,
    target_radii: npt.ArrayLike | float = DEFAULT_RADIUS,
    occluder_radii: npt.ArrayLike | float = DEFAULT_RADIUS,
) -> Scorer:
    """
    Score directions by how many target atoms are directly visible.

    An atom is hidden when a nearer sphere - an occluder or another target
    atom - covers its projected centre. Ties, common when the target is in the
    open and every view sees all of it, are broken towards the view where the
    target's projected bounding box is largest.
    """
    target = as_points(target)
    blockers, radii = _spheres(target, target_radii, occluders, occluder_radii)
    n_target = len(target)
    eps = 1e-6

    def score(directions: Points) -> npt.NDArray[np.float64]:
        out = np.empty(len(directions))
        for k, forward in enumerate(directions):
            t1, t2 = tangent_basis(forward)
            u, v, z = blockers @ t1, blockers @ t2, blockers @ forward
            tu, tv, tz = u[:n_target], v[:n_target], z[:n_target]
            # only spheres in front of the furthest target atom and overlapping
            # the target's footprint can hide anything
            near = (
                (z < tz.max() - eps)
                & (u + radii > tu.min())
                & (u - radii < tu.max())
                & (v + radii > tv.min())
                & (v - radii < tv.max())
            )
            bu, bv, bz, br2 = u[near], v[near], z[near], radii[near] ** 2
            hidden = np.zeros(n_target, dtype=bool)
            for chunk in _chunks(n_target, len(bu)):
                covered = (tu[chunk, None] - bu) ** 2 + (tv[chunk, None] - bv) ** 2
                covered = (covered < br2) & (bz < tz[chunk, None] - eps)
                hidden[chunk] = covered.any(axis=1)
            spread = np.ptp(tu) * np.ptp(tv)
            out[k] = n_target - hidden.sum() + 0.5 * spread / (spread + 1.0)
        return out

    return score


def visible_area_scorer(
    target: npt.ArrayLike,
    occluders: npt.ArrayLike | None = None,
    target_radii: npt.ArrayLike | float = DEFAULT_RADIUS,
    occluder_radii: npt.ArrayLike | float = DEFAULT_RADIUS,
    resolution: int = 64,
) -> Scorer:
    """
    Score directions by the projected area of the target left visible.

    Spheres are rasterised into a depth buffer covering the target's enclosing
    sphere, the same world size from every direction, so the pixel count of
    target spheres winning the depth test measures visible area.
    """
    target = as_points(target)
    spheres, radii = _spheres(target, target_radii, occluders, occluder_radii)
    centre, extent = enclosing_sphere(target)
    half = extent + float(np.max(radii[: len(target)]))
    n_target = len(target)

    def score(directions: Points) -> npt.NDArray[np.float64]:
        out = np.empty(len(directions))
        for k, forward in enumerate(directions):
            winner = rasterize(spheres, radii, forward, centre, half, resolution)
            out[k] = np.count_nonzero((winner >= 0) & (winner < n_target))
        return out

    return score


def entropy_scorer(
    points: npt.ArrayLike,
    radii: npt.ArrayLike | float = DEFAULT_RADIUS,
    labels: npt.ArrayLike | None = None,
    resolution: int = 64,
) -> Scorer:
    """
    Score directions by viewpoint entropy over labelled parts.

    Parameters
    ----------
    points : array_like
        ``(N, 3)`` sphere centres.
    radii : array_like or float
        Sphere radii.
    labels : array_like, optional
        ``(N,)`` integer part for each sphere - residue or chain IDs, say. Each
        sphere is its own part when left out, the analogue of faces in [2]_.
    resolution : int, default 64
        Raster size across the structure's enclosing sphere.
    """
    points = as_points(points)
    radii = _radii(radii, len(points))
    if labels is None:
        labels = np.arange(len(points))
    _, labels = np.unique(np.asarray(labels), return_inverse=True)
    centre, extent = enclosing_sphere(points)
    half = extent + float(radii.max())
    total = resolution * resolution

    def score(directions: Points) -> npt.NDArray[np.float64]:
        out = np.empty(len(directions))
        for k, forward in enumerate(directions):
            winner = rasterize(points, radii, forward, centre, half, resolution)
            hit = winner[winner >= 0]
            areas = np.bincount(labels[hit])
            areas = np.append(areas[areas > 0], total - len(hit))
            p = areas[areas > 0] / total
            out[k] = -np.sum(p * np.log2(p))
        return out

    return score


def ray_scorer(
    trees: Sequence, target: npt.ArrayLike, offset: float = DEFAULT_RADIUS
) -> Scorer:
    """
    Score directions by target atoms with a clear line of sight to the camera.

    Casts a ray from each target atom towards the camera against the actual
    rendered geometry, so a cartoon occludes as a ribbon rather than as the
    atoms it was built from. Ties are broken as in :func:`unobstructed_scorer`.

    Parameters
    ----------
    trees : sequence
        Objects with a ``ray_cast(origin, direction)`` method returning a hit
        location or ``None`` first - ``mathutils.bvhtree.BVHTree`` built from
        the evaluated meshes.
    target : array_like
        ``(N, 3)`` atom positions.
    offset : float
        How far along the ray to start, to clear the atom's own geometry.
    """
    target = as_points(target)

    def score(directions: Points) -> npt.NDArray[np.float64]:
        out = np.empty(len(directions))
        for k, forward in enumerate(directions):
            back = -forward
            clear = 0
            for p in target + back * offset:
                if all(tree.ray_cast(p, back)[0] is None for tree in trees):
                    clear += 1
            t1, t2 = tangent_basis(forward)
            spread = np.ptp(target @ t1) * np.ptp(target @ t2)
            out[k] = clear + 0.5 * spread / (spread + 1.0)
        return out

    return score


def rasterize(
    points: Points,
    radii: npt.NDArray[np.float64],
    forward: npt.ArrayLike,
    centre: npt.ArrayLike,
    half: float,
    resolution: int,
) -> npt.NDArray[np.int64]:
    """
    Orthographic depth raster of spheres: the nearest sphere at each pixel.

    Parameters
    ----------
    points, radii
        Sphere centres and radii.
    forward : array_like
        View direction.
    centre : array_like
        World position at the middle of the raster.
    half : float
        Half the raster's width in world units.
    resolution : int
        Pixels along each side.

    Returns
    -------
    ndarray
        ``(resolution * resolution,)`` index of the front sphere at each pixel,
        ``-1`` where there is none.
    """
    forward = _unit(forward)
    t1, t2 = tangent_basis(forward)
    rel = points - np.asarray(centre)
    u, v, z = rel @ t1, rel @ t2, rel @ forward
    size = 2.0 * half / resolution
    inside = (np.abs(u) < half + radii) & (np.abs(v) < half + radii)
    index = np.flatnonzero(inside)
    u, v, z, r = u[index], v[index], z[index], radii[index]

    reach = int(np.ceil(r.max() / size)) if len(r) else 0
    steps = np.arange(-reach, reach + 1)
    gx = (
        np.floor((u + half) / size).astype(np.int64)[:, None, None]
        + steps[None, :, None]
    )
    gy = (
        np.floor((v + half) / size).astype(np.int64)[:, None, None]
        + steps[None, None, :]
    )
    du = (gx + 0.5) * size - half - u[:, None, None]
    dv = (gy + 0.5) * size - half - v[:, None, None]
    d2 = du * du + dv * dv
    r2 = (r * r)[:, None, None]
    hit = (d2 < r2) & (gx >= 0) & (gx < resolution) & (gy >= 0) & (gy < resolution)

    owner = np.broadcast_to(index[:, None, None], hit.shape)[hit]
    depth = (z[:, None, None] - np.sqrt(np.where(hit, r2 - d2, 0.0)))[hit]
    pixel = (gy * resolution + gx)[hit]

    winner = np.full(resolution * resolution, -1, dtype=np.int64)
    if len(pixel):
        order = np.lexsort((depth, pixel))
        pixel, owner = pixel[order], owner[order]
        first = np.flatnonzero(np.r_[True, pixel[1:] != pixel[:-1]])
        winner[pixel[first]] = owner[first]
    return winner


def best_view(
    target: npt.ArrayLike,
    method: str = "unobstructed",
    occluders: npt.ArrayLike | None = None,
    target_radii: npt.ArrayLike | float = DEFAULT_RADIUS,
    occluder_radii: npt.ArrayLike | float = DEFAULT_RADIUS,
    labels: npt.ArrayLike | None = None,
    up: npt.ArrayLike = _WORLD_UP,
    resolution: int = 256,
) -> Points:
    """
    A camera basis looking at ``target`` from the best direction for ``method``.

    Parameters
    ----------
    target : array_like
        ``(N, 3)`` positions of what should be seen.
    method : str, default "unobstructed"
        One of :data:`METHODS`; see the module docs.
    occluders : array_like, optional
        ``(M, 3)`` positions of everything else that could get in the way.
        Ignored by ``"pca"`` and ``"entropy"``.
    target_radii, occluder_radii : array_like or float
        Sphere radii in world units.
    labels : array_like, optional
        Parts for ``"entropy"``.
    up : array_like, default (0, 0, 1)
        Which way is up in the result. ``"pca"`` sets its own.
    resolution : int, default 256
        Directions in the coarse search.

    Returns
    -------
    ndarray
        ``(3, 3)`` rows ``(right, up, forward)``.
    """
    if method == "pca":
        return pca_view(target)
    if method == "unobstructed":
        scorer = unobstructed_scorer(target, occluders, target_radii, occluder_radii)
    elif method == "visible_area":
        scorer = visible_area_scorer(target, occluders, target_radii, occluder_radii)
    elif method == "entropy":
        scorer = entropy_scorer(target, target_radii, labels)
    else:
        raise ValueError(f"Unknown method {method!r}, expected one of {METHODS}.")
    return look_along(search(scorer, resolution), up)


def _unit(vector: npt.ArrayLike) -> Points:
    vector = np.asarray(vector, dtype=np.float64)
    norm = np.linalg.norm(vector)
    if norm == 0:
        raise ValueError("Cannot normalise a zero-length direction.")
    return vector / norm


def _radii(radii: npt.ArrayLike | float, n: int) -> npt.NDArray[np.float64]:
    radii = np.asarray(radii, dtype=np.float64)
    if radii.ndim == 0:
        return np.full(n, float(radii))
    if radii.shape != (n,):
        raise ValueError(f"Expected {n} radii, got an array of shape {radii.shape}.")
    return radii


def _spheres(target, target_radii, occluders, occluder_radii):
    """Target spheres first, then occluders, as one set of centres and radii."""
    radii = _radii(target_radii, len(target))
    if occluders is None or len(occluders) == 0:
        return target, radii
    occluders = as_points(occluders)
    return (
        np.concatenate([target, occluders]),
        np.concatenate([radii, _radii(occluder_radii, len(occluders))]),
    )


def _chunks(n_rows: int, n_cols: int, budget: int = 2_000_000):
    """Row slices keeping an ``n_rows x n_cols`` comparison within ``budget``."""
    step = max(1, budget // max(n_cols, 1))
    for start in range(0, n_rows, step):
        yield slice(start, min(start + step, n_rows))
