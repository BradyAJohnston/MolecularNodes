import bpy
from pathlib import Path
from databpy.material import append_from_blend
from ..assets import MN_DATA_FILE
from typing import TypeVar
import numpy as np

MATERIAL_NAMES = [
    "MN Default",
    "MN Flat Outline",
    "MN Squishy",
    "MN Transparent Outline",
    "MN Ambient Occlusion",
]


def path_resolve(path: str | Path) -> Path:
    if isinstance(path, str):
        return Path(bpy.path.abspath(path))
    elif isinstance(path, Path):
        return Path(bpy.path.abspath(str(path)))
    else:
        raise ValueError(f"Unable to resolve path: {path}")


def append_material(name: str) -> bpy.types.Material:
    "Append a material from the MN_DATA_FILE."
    return append_from_blend(name, str(MN_DATA_FILE))


def add_all_materials() -> None:
    "Append all pre-defined materials from the MN_DATA_FILE."
    for name in MATERIAL_NAMES:
        append_material(name)


T = TypeVar(
    "T",
    float,
    int,
    bool,
    str,
    np.ndarray,
    tuple[float, float, float],
    list[float, float, float],
)


def socket(node_name: str, socket: str | int, type_: type[T]):
    def getter(self) -> float | int | bool | str | np.ndarray:
        value = self.nodes[node_name].inputs[socket].default_value
        if isinstance(value, (int, float, bool)):
            return value
        elif all(isinstance(x, float) for x in value):
            return np.array(value, dtype=float)
        elif all(isinstance(x, int) for x in value):
            return np.array(value, dtype=int)
        else:
            raise ValueError(f"Unable to convert {value} to {type_}")

    def setter(self, value: T) -> None:
        self.nodes[node_name].inputs[socket].default_value = value

    return property(getter, setter)


def option(node_name: str, input: str | int, type_: type[T]):
    def getter(self) -> T:
        return getattr(self.nodes[node_name], input)

    def setter(self, value: T) -> None:
        setattr(self.nodes[node_name], input, value)

    return property(getter, setter)
