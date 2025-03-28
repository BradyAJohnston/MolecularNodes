import bpy
from pathlib import Path
from typing import TypeVar
import numpy as np
from .nodes import input_named_attribute
from .arrange import arrange_tree


def path_resolve(path: str | Path) -> Path:
    if isinstance(path, str):
        return Path(bpy.path.abspath(path))
    elif isinstance(path, Path):
        return Path(bpy.path.abspath(str(path)))
    else:
        raise ValueError(f"Unable to resolve path: {path}")


class TreeInterface:
    @property
    def node_tree(self):
        raise NotImplementedError("Must be implemented by subclass")

    @property
    def nodes(self):
        return self.node_tree.nodes

    @property
    def links(self):
        return self.node_tree.links


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


class SocketLinkedError(Exception):
    def __init__(self, message="Socket is linked, default value is not used."):
        self.message = message
        super().__init__(self.message)


def _socket_bool(socket: bpy.types.NodeSocket):
    def getter(self) -> bool:
        if socket.is_linked:
            raise SocketLinkedError()
        else:
            return getattr(socket, "default_value")

    def setter(self, value: bool | str) -> None:
        if socket.is_linked:
            socket.node.id_data.links.remove(socket.links[0])
            arrange_tree(socket.node.id_data)
        if isinstance(value, bool):
            setattr(socket, "default_value", value)

        if isinstance(value, str):
            interface_item = socket.node.node_tree.interface.items_tree[socket.name]
            if interface_item.force_non_field:
                raise ValueError(
                    f"Cannot use a named attribute for this input: '{socket}' as it accepts only single values"
                )
            input_named_attribute(socket, value, "BOOLEAN")
            arrange_tree(socket.node.id_data)

    return property(getter, setter)


def socket(socket: bpy.types.NodeSocket, type_: type[T]):
    if type_ is bool:
        return _socket_bool(socket)

    def getter(self) -> float | int | bool | str | np.ndarray:
        value = getattr(socket, "default_value")
        if isinstance(value, (int, float, bool)):
            return value
        elif isinstance(value, str):
            return value
        elif all(isinstance(x, float) for x in value):
            return np.array(value, dtype=float)
        elif all(isinstance(x, int) for x in value):
            return np.array(value, dtype=int)
        else:
            raise ValueError(f"Unable to convert {value} to {type_}")

    def setter(self, value: T) -> None:
        setattr(socket, "default_value", value)

    return property(getter, setter)


def option(node_name: str, input: str | int, type_: type[T]):
    def getter(self) -> T:
        return getattr(self.nodes[node_name], input)

    def setter(self, value: T) -> None:
        setattr(self.nodes[node_name], input, value)

    return property(getter, setter)
