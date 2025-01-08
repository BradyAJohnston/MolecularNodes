import os
import sys
import numpy as np

from pathlib import Path
from math import floor
from mathutils import Matrix

ADDON_DIR = Path(__file__).resolve().parent
MN_DATA_FILE = os.path.join(ADDON_DIR, "assets", "MN_data_file_4.2.blend")


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
    mapping: np.ndarray | None = None,
) -> int:
    frame = max(frame - offset, 0)

    if mapping is not None:
        if not isinstance(mapping, np.ndarray):
            raise ValueError(
                "Frame mapping must be an array of values to map frames to"
            )
        # add the subframes to the frame mapping
        frame_map = np.repeat(mapping, subframes + 1)
        # get the current and next frames
        frame_a = frame_map[frame]

    frame_a = frame

    if subframes > 0:
        frame_a = int(frame / (subframes + 1))

    return frame_a


def frames_to_average(frame: int, average: int = 0, lower_bound: int = 0) -> np.ndarray:
    length = average * 2 + 1
    frames = np.arange(length) + frame - average
    frames = frames[frames >= lower_bound]
    return frames


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
