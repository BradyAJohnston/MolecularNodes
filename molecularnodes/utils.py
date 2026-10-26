import io
import json
import os
import sys
from contextlib import ExitStack
from typing import List, cast
import bpy
import numpy as np

_INCREASED_CLIP_END = 1e4
_DEFAULT_CAMERA_CLIP_END = 1e2
_DEFAULT_VIEWPORT_CLIP_END = 1e3


class Unset:
    """Sentinel for an argument that was not passed (``None`` is a meaningful value)."""

    def __repr__(self) -> str:
        return "<unset>"


_UNSET = Unset()


def _increase_view_distance():
    scene = bpy.context.scene
    assert scene
    if scene.camera is not None:
        camera = cast(bpy.types.Camera, scene.camera.data)
        if camera.clip_end == _DEFAULT_CAMERA_CLIP_END:
            camera.clip_end = _INCREASED_CLIP_END

    # viewport clip end is stored per 3D Viewport space, so update every
    # VIEW_3D area across all screens (covers currently inactive workspaces)
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type != "VIEW_3D":
                continue
            for space in area.spaces:
                if space.type == "VIEW_3D":
                    space = cast(bpy.types.SpaceView3D, space)
                    if space.clip_end == _DEFAULT_VIEWPORT_CLIP_END:
                        space.clip_end = _INCREASED_CLIP_END


def fraction(x, y):
    return x % y / y


def correct_periodic_1d(
    value1: np.ndarray, value2: np.ndarray, boundary: float
) -> np.ndarray:
    # shift value2 by whole box lengths to the periodic image nearest value1.
    # A pure function of its inputs: mutating value2 in place would poison the
    # frame position cache, chaining each frame's correction onto the previous
    # already-shifted frame so positions unwrap further out every crossing
    diff = value2 - value1
    return value2 - boundary * np.round(diff / boundary)


def correct_periodic_positions(
    positions_1: np.ndarray, positions_2: np.ndarray, dimensions: np.ndarray
) -> np.ndarray:
    if not np.allclose(dimensions[3:], 90.0):
        raise ValueError(
            f"Only works with orthorhombic unitcells, and not dimensions={dimensions}"
        )
    final_positions = positions_2.copy()
    for i in range(3):
        final_positions[:, i] = correct_periodic_1d(
            positions_1[:, i], positions_2[:, i], dimensions[i]
        )
    return final_positions


def is_periodic_box(dimensions: np.ndarray | None) -> bool:
    """Whether unit cell dimensions describe a usable periodic box.

    Parameters
    ----------
    dimensions : np.ndarray or None
        Unit cell as ``[lx, ly, lz, alpha, beta, gamma]``, as returned by
        `MDAnalysis.coordinates.timestep.Timestep.dimensions`.

    Returns
    -------
    bool
        False when there is no box, when a length is not positive, when the
        angles do not form a valid cell, or for the ``[1, 1, 1, 90, 90, 90]``
        placeholder that some PDB writers use (which MDAnalysis's DSSP also
        ignores).
    """
    from MDAnalysis.lib.mdamath import triclinic_vectors

    if dimensions is None:
        return False
    dimensions = np.asarray(dimensions, dtype=np.float64)
    if dimensions.shape != (6,) or not np.all(np.isfinite(dimensions)):
        return False
    if np.any(dimensions[:3] <= 0) or np.allclose(dimensions, [1, 1, 1, 90, 90, 90]):
        return False
    # MDAnalysis returns all-zero vectors for angles that cannot form a cell
    return bool(np.any(triclinic_vectors(dimensions)))


def minimum_image_positions(
    positions_1: np.ndarray, positions_2: np.ndarray, dimensions: np.ndarray
) -> np.ndarray:
    """Shift each of positions_2 to its periodic image nearest positions_1.

    Triclinic boxes use the minimum image convention of
    `MDAnalysis.lib.distances.minimize_vectors`. Orthorhombic boxes give the
    same result with `correct_periodic_positions`, which is several times
    faster.

    Parameters
    ----------
    positions_1 : np.ndarray
        Reference positions, shape (n, 3).
    positions_2 : np.ndarray
        Positions to shift, shape (n, 3).
    dimensions : np.ndarray
        Unit cell as ``[lx, ly, lz, alpha, beta, gamma]``, with lengths in the
        same units as the positions.

    Returns
    -------
    np.ndarray
        A new array: positions_2 shifted by whole box vectors.
    """
    from MDAnalysis.lib.distances import minimize_vectors

    if np.allclose(dimensions[3:], 90.0):
        return correct_periodic_positions(positions_1, positions_2, dimensions[:3])
    vectors = positions_2 - positions_1
    # shift positions_2 rather than rebuilding it from positions_1, so atoms
    # that did not cross a boundary keep their positions (up to rounding of
    # the zero shift)
    return positions_2 + (minimize_vectors(vectors, dimensions) - vectors)


def frame_mapper(
    frame: int,
    subframes: int = 0,
    offset: int = 0,
    mapping: np.ndarray | List | None = None,
) -> int:
    frame = max(frame - offset, 0)

    if mapping is not None:
        if isinstance(mapping, list):
            mapping = np.array(mapping)
        if not isinstance(mapping, np.ndarray):
            raise ValueError(
                "Frame mapping must be an array of values to map frames to"
            )
        # add the subframes to the frame mapping
        frame_map = np.repeat(mapping, subframes + 1)
        # get the current and next frames
        frame_a = frame_map[frame]
    else:
        frame_a = frame

        if subframes > 0:
            frame_a = int(frame / (subframes + 1))

    return frame_a


def frames_to_average(
    frame: int, upper_bound: int, average: int = 0, lower_bound: int = 0
) -> np.ndarray:
    length = average * 2 + 1
    frames = np.arange(length) + frame - average
    frames = frames[frames >= lower_bound]
    frames = frames[frames < upper_bound]
    return frames


def pca_dejitter(frames: np.ndarray, index: int) -> np.ndarray:
    """Smooth one frame of a window by PCA across the window's frames.

    Each frame of the window, flattened, is one sample. The frame at `index` is
    projected onto the window's principal components, the weights of the
    components are faded out linearly from the first (kept in full) to the last
    (dropped), and the frame is rebuilt from them. The dominant motion of the
    window survives while frame-to-frame jitter, spread over the later
    components, is removed. With every weight at zero this is the window mean.

    Adapted from the PCA De-Jitter tutorial for Houdini by SideFX:
    https://www.sidefx.com/tutorials/pca-dejitter/

    Parameters
    ----------
    frames : np.ndarray
        Values for each frame of the window, shape `(n_frames, ...)`.
    index : int
        Index into `frames` of the frame to smooth.

    Returns
    -------
    np.ndarray
        The smoothed frame, with the shape and dtype of `frames[index]`.
    """
    samples = frames.reshape(len(frames), -1).astype(np.float64)
    mean = samples.mean(axis=0)
    _, _, components = np.linalg.svd(samples - mean, full_matrices=False)
    # centring leaves at most n_frames - 1 components that carry anything
    components = components[: len(frames) - 1]
    weights = components @ (samples[index] - mean)
    weights *= np.linspace(1.0, 0.0, len(components))
    smoothed = mean + weights @ components
    return smoothed.reshape(frames.shape[1:]).astype(frames.dtype)


def rigid_group_fit(
    current: np.ndarray,
    target: np.ndarray,
    groups: np.ndarray,
    reference: np.ndarray,
    extra_reference: tuple[np.ndarray, np.ndarray] | None = None,
) -> np.ndarray:
    """Move each group of `current` rigidly onto `target`.

    For every group the rotation and translation that best superimpose its
    `reference` points in `current` onto the same points in `target` (Kabsch)
    are applied to all of the group's points. The shape of each group is kept
    from `current` while its placement follows `target`, so averaging `target`
    over frames can't shorten bonds or pass atoms of a group through each other.

    Parameters
    ----------
    current : np.ndarray
        Positions to move, shape `(n, 3)`.
    target : np.ndarray
        Positions to fit onto, shape `(n, 3)`.
    groups : np.ndarray
        Group of each point as integers `0..n_groups-1`, shape `(n,)`.
    reference : np.ndarray
        Boolean mask of the points each group is fitted on, shape `(n,)`. Every
        group needs at least one; with one the group is only translated.
    extra_reference : tuple[np.ndarray, np.ndarray] | None, optional
        Points of other groups to also fit a group on, as `(point indices,
        groups)`, such as the next residue's N so neighbouring residues agree
        on the bond between them.

    Returns
    -------
    np.ndarray
        The moved positions, with the shape and dtype of `current`.
    """
    n_groups = int(groups.max()) + 1
    ref_index = np.flatnonzero(reference)
    ref_groups = groups[ref_index]
    if extra_reference is not None:
        ref_index = np.concatenate([ref_index, extra_reference[0]])
        ref_groups = np.concatenate([ref_groups, extra_reference[1]])
    counts = np.bincount(ref_groups, minlength=n_groups)[:, None]
    p = current[ref_index].astype(np.float64)
    q = target[ref_index].astype(np.float64)

    def centroid(values):
        return np.stack(
            [np.bincount(ref_groups, values[:, i], n_groups) for i in range(3)],
            axis=1,
        ) / np.maximum(counts, 1)

    p_centre = centroid(p)
    q_centre = centroid(q)
    p = p - p_centre[ref_groups]
    q = q - q_centre[ref_groups]
    covariance = np.empty((n_groups, 3, 3))
    for i in range(3):
        for j in range(3):
            covariance[:, i, j] = np.bincount(ref_groups, p[:, i] * q[:, j], n_groups)
    u, _, vt = np.linalg.svd(covariance)
    # flip the last axis where needed so each fit is a rotation, not a reflection
    sign = np.sign(np.linalg.det(np.einsum("gji,gkj->gik", vt, u)))
    vt[:, 2, :] *= np.where(sign == 0, 1.0, sign)[:, None]
    rotation = np.einsum("gji,gkj->gik", vt, u)
    moved = (
        np.einsum("nij,nj->ni", rotation[groups], current - p_centre[groups])
        + q_centre[groups]
    )
    return moved.astype(current.dtype)


# data types for the np.array that will store per-chain symmetry operations


def array_transforms_from_dict(transforms_dict):
    if isinstance(transforms_dict, str):
        transforms_dict = json.loads(transforms_dict.replace("nan", "0.0"))

    dtype = [
        ("assembly_id", int),
        ("sym_id", int),
        ("chain_id", "U10"),
        ("transform", float, (4, 4)),
        ("pdb_model_num", int),
    ]

    transforms = []
    for i, assembly in enumerate(transforms_dict.values()):
        for j, transform in enumerate(assembly):
            chains = transform["chain_ids"]
            arr = np.zeros((len(chains)), dtype=dtype)
            arr["assembly_id"] = i + 1
            arr["sym_id"] = j
            arr["chain_id"] = chains
            arr["transform"] = transform["matrix"]
            arr["pdb_model_num"] = transform["pdb_model_num"]
            transforms.append(arr)

    return np.hstack(transforms)


def count_value_changes(arr1, arr2):
    """
    Counts the number of times either array changes value (increases or decreases).

    Used for computing the `ures_id` attribute for imported structures.

    Parameters:
    -----------
    arr1 : numpy.ndarray
        First integer array
    arr2 : numpy.ndarray
        Second integer array

    Returns:
    --------
    numpy.ndarray
        Array of cumulative changes count
    """
    diff1 = np.diff(arr1) != 0
    diff2 = np.diff(arr2) != 0

    combined_changes = np.logical_or(diff1, diff2)
    result = np.zeros(len(arr1), dtype=int)
    result[1:] = np.cumsum(combined_changes)

    return result


class suppress_stdout(object):
    """
    A context manager that blocks stdout for its scope, usage:
        with suppress_stdout():
            do_something()

    Suppression is skipped when the ``MN_VERBOSE`` environment variable is
    set, letting Blender's render output through for debugging.

    From: https://stackoverflow.com/a/14797594

    """

    def __init__(self, *args, **kw):
        self._suppress = not os.environ.get("MN_VERBOSE")
        if not self._suppress:
            return
        sys.stdout.flush()
        self._origstdout = sys.stdout
        try:
            stdout_fno = sys.stdout.fileno()
        except (AttributeError, io.UnsupportedOperation, ValueError):
            # sys.stdout has been replaced with a stream that has no real
            # file descriptor (e.g. ipykernel in Jupyter / Colab), but the
            # C-level stdout that Blender writes to is still fd 1
            stdout_fno = 1
        try:
            self._oldstdout_fno = os.dup(stdout_fno)
            self._devnull = os.open(os.devnull, os.O_WRONLY)
        except OSError:
            # no usable stdout fd at all; skip suppression rather than fail
            self._suppress = False

    def __enter__(self):
        if not self._suppress:
            return
        self._newstdout = os.dup(1)
        os.dup2(self._devnull, 1)
        os.close(self._devnull)
        self._fdopen_stdout = os.fdopen(self._newstdout, "w")
        sys.stdout = self._fdopen_stdout

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self._suppress:
            return
        sys.stdout = self._origstdout
        sys.stdout.flush()
        self._fdopen_stdout.close()
        os.dup2(self._oldstdout_fno, 1)
        os.close(self._oldstdout_fno)


class temp_override_property:
    """
    A context manager to temporarily override an object property

    """

    def __init__(self, obj, prop_name, value):
        self.obj = obj
        self.prop_name = prop_name
        self.value = value
        self.orig_value = _UNSET

    def __enter__(self):
        if not hasattr(self.obj, self.prop_name):
            raise ValueError(f"Property {self.prop_name} not found")
        self.orig_value = getattr(self.obj, self.prop_name)
        setattr(self.obj, self.prop_name, self.value)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.orig_value is not _UNSET:
            setattr(self.obj, self.prop_name, self.orig_value)


def temp_override_properties(stack: ExitStack, properties: list) -> None:
    """
    Add multiple property overrides to an ExitStack

    """
    for property in properties:
        stack.enter_context(temp_override_property(*property))
