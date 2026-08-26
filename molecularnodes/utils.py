import json
import os
import sys
from contextlib import ExitStack
from typing import List, cast
import addon_utils
import bpy
import numpy as np

_INCREASED_CLIP_END = 1e4
_DEFAULT_CAMERA_CLIP_END = 1e2
_DEFAULT_VIEWPORT_CLIP_END = 1e3


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


def load_extension_module():
    # check enabled addons
    for addon in bpy.context.preferences.addons.keys():
        if addon.endswith(".molecularnodes"):
            is_enabled, is_loaded = addon_utils.check(addon)
            if is_enabled and is_loaded:
                return sys.modules[addon]
    # return this parent module
    mn_module_name = __name__.rsplit(".", 1)[0]
    return sys.modules[mn_module_name]


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
        self._oldstdout_fno = os.dup(sys.stdout.fileno())
        self._devnull = os.open(os.devnull, os.O_WRONLY)

    def __enter__(self):
        if not self._suppress:
            return
        self._newstdout = os.dup(1)
        os.dup2(self._devnull, 1)
        os.close(self._devnull)
        sys.stdout = os.fdopen(self._newstdout, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self._suppress:
            return
        sys.stdout = self._origstdout
        sys.stdout.flush()
        os.dup2(self._oldstdout_fno, 1)


class temp_override_property:
    """
    A context manager to temporarily override an object property

    """

    _UNSET = object()

    def __init__(self, obj, prop_name, value):
        self.obj = obj
        self.prop_name = prop_name
        self.value = value
        self.orig_value = self._UNSET

    def __enter__(self):
        if not hasattr(self.obj, self.prop_name):
            raise ValueError(f"Property {self.prop_name} not found")
        self.orig_value = getattr(self.obj, self.prop_name)
        setattr(self.obj, self.prop_name, self.value)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.orig_value is not self._UNSET:
            setattr(self.obj, self.prop_name, self.orig_value)


def temp_override_properties(stack: ExitStack, properties: list) -> None:
    """
    Add multiple property overrides to an ExitStack

    """
    for property in properties:
        stack.enter_context(temp_override_property(*property))
