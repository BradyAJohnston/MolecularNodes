import bpy
from typing import TypeVar
import numpy as np
from .arrange import arrange_tree
from .nodes import input_named_attribute


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
            return socket.default_value

    def setter(self, value: bool | str) -> None:
        if socket.is_linked:
            socket.node.id_data.links.remove(socket.links[0])
            arrange_tree(socket.node.id_data)
        if isinstance(value, bool):
            socket.default_value = value

        if isinstance(value, str):
            interface_item = socket.node.node_tree.interface.items_tree[socket.name]
            if interface_item.force_non_field:
                raise ValueError(
                    f"Cannot use a named attribute for this input: '{socket}' as it accepts only single values"
                )
            input_named_attribute(socket, value, "BOOLEAN")
            arrange_tree(socket.node.id_data)

    return property(getter, setter)


def socket_lookup(node_name: str, socket_name: str, type_: type[T]):
    def getter(self) -> type[T]:
        return self.nodes[node_name].inputs[socket_name].default_value

    def setter(self, value: T) -> None:
        self.nodes[node_name].inputs[socket_name].default_value = value

    return property(getter, setter)


def check_linked(socket: bpy.types.NodeSocket):
    if socket.is_linked:
        raise SocketLinkedError()


def remove_linked(socket: bpy.types.NodeSocket):
    if socket.is_linked:
        socket.node.id_data.links.remove(socket.links[0])
        arrange_tree(socket.node.id_data)


def getset_float(socket: bpy.types.NodeSocket):
    def getter(self) -> float:
        check_linked(socket)
        return getattr(socket, "default_value")

    def setter(self, value: float | str) -> None:
        remove_linked(socket)
        if isinstance(value, str):
            input_named_attribute(socket, value, "FLOAT")
        elif isinstance(value, float):
            setattr(socket, "default_value", value)
        else:
            raise ValueError(f"Value must be a float, not a {type(value)}")

    return property(getter, setter)


def getset_bool(socket: bpy.types.NodeSocket):
    def getter(self) -> bool:
        check_linked(socket)
        return getattr(socket, "default_value")

    def setter(self, value: bool) -> None:
        remove_linked(socket)
        if isinstance(value, bool):
            setattr(socket, "default_value", value)
        if isinstance(value, str):
            input_named_attribute(socket, value, "BOOLEAN")
        else:
            raise ValueError(f"Value must be a bool, not a {type(value)}")

    return property(getter, setter)


def getset_vector(socket: bpy.types.NodeSocket):
    def getter(self) -> np.ndarray:
        check_linked(socket)
        return getattr(socket, "default_value")

    def setter(self, value: np.ndarray) -> None:
        remove_linked(socket)
        if isinstance(value, (np.ndarray, list, tuple)):
            if len(value) != 3:
                raise ValueError(
                    f"Value must be a numpy array of shape (3,), not a {value.shape}"
                )
            setattr(socket, "default_value", value)
        elif isinstance(value, str):
            input_named_attribute(socket, value, "FLOAT_VECTOR")
        else:
            raise ValueError(f"Value must be a numpy array, not a {type(value)}")


def socket(socket: bpy.types.NodeSocket):
    if socket.type == "BOOLEAN":
        return _socket_bool(socket)

    def getter(self) -> float | int | bool | str | np.ndarray:
        value = socket.default_value
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
        socket.default_value = value

    return property(getter, setter)


def option(node_name: str, input: str | int, type_: type[T]):
    def getter(self) -> T:
        return getattr(self.nodes[node_name], input)

    def setter(self, value: T) -> None:
        setattr(self.nodes[node_name], input, value)

    return property(getter, setter)
