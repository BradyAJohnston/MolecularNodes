import os
from pathlib import Path

import numpy as np
from mathutils import Matrix

ADDON_DIR = Path(__file__).resolve().parent
MN_DATA_FILE = os.path.join(ADDON_DIR, "assets", "MN_data_file_4.2.blend")


def correct_periodic_1d(
    value1: np.ndarray, value2: np.ndarray, boundary: float
) -> np.ndarray:
    diff = value2 - value1
    half = boundary / 2
    value2[diff > half] -= boundary
    value2[diff < -half] += boundary
    return value2


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


def lerp(a: np.ndarray, b: np.ndarray, t: float = 0.5) -> np.ndarray:
    """
    Linearly interpolate between two values.

    Parameters
    ----------
    a : array_like
        The starting value.
    b : array_like
        The ending value.
    t : float, optional
        The interpolation parameter. Default is 0.5.

    Returns
    -------
    array_like
        The interpolated value(s).

    Notes
    -----
    This function performs linear interpolation between `a` and `b` using the
    interpolation parameter `t` such that the result lies between `a` and `b`.

    Examples
    --------
    >>> lerp(1, 2, 0.5)
    1.5

    >>> lerp(3, 7, 0.2)
    3.8

    >>> lerp([1, 2, 3], [4, 5, 6], 0.5)
    array([2.5, 3.5, 4.5])

    """
    return np.add(a, np.multiply(np.subtract(b, a), t))


# data types for the np.array that will store per-chain symmetry operations
dtype = [
    ("assembly_id", int),
    ("transform_id", int),
    ("chain_id", "U10"),
    ("rotation", float, 4),  # quaternion form
    ("translation", float, 3),
]


def array_quaternions_from_dict(transforms_dict):
    n_transforms = 0
    for assembly in transforms_dict.values():
        for transform in assembly:
            n_transforms += len(transform[0])

    arr = np.array((n_transforms), dtype=dtype)

    transforms = []
    for i, assembly in enumerate(transforms_dict.values()):
        for j, transform in enumerate(assembly):
            chains = transform[0]
            matrix = transform[1]
            arr = np.zeros((len(chains)), dtype=dtype)
            translation, rotation, scale = Matrix(matrix).decompose()
            arr["assembly_id"] = i + 1
            arr["transform_id"] = j
            arr["chain_id"] = chains
            arr["rotation"] = rotation
            arr["translation"] = translation
            transforms.append(arr)

    return np.hstack(transforms)
