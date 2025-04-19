import json
import sys
from typing import List
import numpy as np
from mathutils import Matrix
from .assets import ADDON_DIR


def add_current_module_to_path():
    path = str(ADDON_DIR.parent)
    sys.path.append(path)


def fraction(x, y):
    return x % y / y


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


# data types for the np.array that will store per-chain symmetry operations


def array_quaternions_from_dict(transforms_dict):
    n_transforms = 0

    if isinstance(transforms_dict, str):
        transforms_dict = json.loads(transforms_dict.replace("nan", "0.0"))

    for assembly in transforms_dict.values():
        # add the number of chains for each transform together, this gives us the total
        # number of transforms we need to initialise
        n_transforms += sum(len(transform["chain_ids"]) for transform in assembly)

    dtype = [
        ("assembly_id", int),
        ("transform_id", int),
        ("chain_id", "U10"),
        ("rotation", float, 4),  # quaternion form
        ("translation", float, 3),
        ("pdb_model_num", int),
    ]

    arr = np.array((n_transforms), dtype=dtype)

    transforms = []
    for i, assembly in enumerate(transforms_dict.values()):
        for j, transform in enumerate(assembly):
            chains = transform["chain_ids"]
            arr = np.zeros((len(chains)), dtype=dtype)
            matrix = transform["matrix"]
            translation, rotation, scale = Matrix(matrix).decompose()
            arr["assembly_id"] = i + 1
            arr["transform_id"] = j
            arr["chain_id"] = chains
            arr["rotation"] = rotation
            arr["translation"] = translation
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
