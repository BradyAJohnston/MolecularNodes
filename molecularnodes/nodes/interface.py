import bpy
from typing import TypeVar
import numpy as np
from .arrange import arrange_tree
from mathutils import Vector

NODE_SPACING = 250


import bpy
from typing import TypeVar, Dict, Any
import numpy as np
from .arrange import arrange_tree
from mathutils import Vector

NODE_SPACING = 250


class TreeInterface:
    def __init__(self):
        self._dynamic_properties = set()

    @property
    def node_tree(self):
        raise NotImplementedError("Must be implemented by subclass")

    @property
    def nodes(self):
        return self.node_tree.nodes

    @property
    def links(self):
        return self.node_tree.links

    def __getattr__(self, name):
        # This method is only called when the attribute doesn't exist
        # No need to check if it exists in _dynamic_properties
        raise AttributeError(
            f"'{self.__class__.__name__}' has no attribute '{name}'"
            f"potentially a dynamic property in {self._dynamic_properties}"
        )

    def __setattr__(self, name, value):
        # For setting attributes, we need to check if it's allowed
        if (
            name != "_dynamic_properties"
            and not name.startswith("_")
            and hasattr(self, "_dynamic_properties")
            and name not in self._dynamic_properties
            and name
            not in ("node_tree", "nodes", "links", "remove", "tree", "material")
        ):
            raise AttributeError(
                f"Cannot set non-existent property '{name}' on '{self.__class__.__name__}'"
                f"\nAvailable dynamic properties: {self._dynamic_properties}"
            )
        # Use the normal attribute setting mechanism
        super().__setattr__(name, value)

    def _register_property(self, name):
        """Register a dynamic property that can be set on this interface"""
        self._dynamic_properties.add(name)

    def _register_properties(self, names):
        """Register multiple dynamic properties at once"""
        self._dynamic_properties.update(names)


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


def input_named_attribute(
    socket: bpy.types.NodeSocket, name: str, data_type: str | None
) -> bpy.types.GeometryNode:
    """
    Add a named attribute node to the tree and connect it to the given socket
    """
    tree = socket.node.id_data
    node_na = tree.nodes.new("GeometryNodeInputNamedAttribute")

    if data_type is not None:
        node_na.data_type = data_type
    node_na.inputs["Name"].default_value = name
    node_na.location = socket.node.location - Vector([NODE_SPACING, 0])

    tree.links.new(node_na.outputs["Attribute"], socket)

    return node_na


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

    return property(getter, setter, doc=socket.description)


def getset_int(
    socket: bpy.types.NodeSocketInt
    | bpy.types.NodeSocketIntFactor
    | bpy.types.NodeSocketIntPercentage
    | bpy.types.NodeSocketIntUnsigned,
):
    def getter(self) -> int:
        check_linked(socket)
        return getattr(socket, "default_value")

    def setter(self, value: int | str) -> None:
        remove_linked(socket)
        if isinstance(value, str):
            input_named_attribute(socket, value, "INT")
        elif isinstance(value, int):
            setattr(socket, "default_value", value)
        else:
            raise ValueError(f"Value must be an int, not a {type(value)}")

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
            remove_linked(socket)
            input_named_attribute(socket, value, "FLOAT_VECTOR")
        else:
            raise ValueError(f"Value must be a numpy array, not a {type(value)}")

    return property(getter, setter)


def getset_string(socket: bpy.types.NodeSocket):
    def getter(self) -> str:
        check_linked(socket)
        return getattr(socket, "default_value")

    def setter(self, value: str) -> None:
        remove_linked(socket)
        if isinstance(value, str):
            setattr(socket, "default_value", value)
        else:
            raise ValueError(f"Value must be a string, not a {type(value)}")

    return property(getter, setter)


def getset_color(socket: bpy.types.NodeSocketColor):
    def getter(self) -> np.ndarray:
        check_linked(socket)
        return getattr(socket, "default_value")

    def setter(self, value: np.ndarray) -> None:
        remove_linked(socket)
        if isinstance(value, (np.ndarray, list, tuple)):
            if len(value) != 4:
                raise ValueError(
                    f"Value must be a numpy array of shape (4,), not a {value.shape}"
                )
            setattr(socket, "default_value", value)
        elif isinstance(value, str):
            remove_linked(socket)
            input_named_attribute(socket, value, "FLOAT_VECTOR")
        else:
            raise ValueError(f"Value must be a numpy array, not a {type(value)}")

    return property(getter, setter)


def getset_rotation(socket: bpy.types.NodeSocketRotation):
    def getter(self) -> np.ndarray:
        check_linked(socket)
        return getattr(socket, "default_value")

    def setter(self, value: np.ndarray) -> None:
        remove_linked(socket)
        if isinstance(value, (np.ndarray, list, tuple)):
            if len(value) != 4:
                raise ValueError(
                    f"Value must be a numpy array of shape (4,), not a {value.shape}"
                )
            setattr(socket, "default_value", value)
        elif isinstance(value, str):
            remove_linked(socket)
            input_named_attribute(socket, value, "QUATERNION")
        else:
            raise ValueError(f"Value must be a numpy array, not a {type(value)}")

    return property(getter, setter)


def socket(socket: bpy.types.NodeSocket):
    if socket.type == "BOOLEAN":
        return _socket_bool(socket)

    match socket.type:
        case "FLOAT":
            return getset_float(socket)
        case "VALUE":
            return getset_float(socket)
        case "INT":
            return getset_int(socket)
        case "BOOLEAN":
            return getset_bool(socket)
        case "ROTATION":
            return getset_rotation(socket)
        case "VECTOR":
            return getset_vector(socket)
        case "STRING":
            return getset_string(socket)
        case "RGBA":
            return getset_color(socket)
        case _:
            raise ValueError(f"Unknown socket type: {socket.type}")

    # def getter(self) -> float | int | bool | str | np.ndarray:
    #     value = getattr(socket, "default_value")

    #     if isinstance(value, (int, float, bool)):
    #         return value
    #     elif isinstance(value, str):
    #         return value
    #     elif all(isinstance(x, float) for x in value):
    #         return np.array(value, dtype=float)
    #     elif all(isinstance(x, int) for x in value):
    #         return np.array(value, dtype=int)
    #     else:
    #         raise ValueError(f"Unable to convert {value} to ")

    # def setter(self, value: T) -> None:
    #     setattr(socket, "default_value", value)

    # return property(getter, setter)


def option(node_name: str, input: str, type_: type[T]):
    def getter(self) -> T:
        return self.nodes[node_name][input]

    def setter(self, value: T) -> None:
        setattr(self.nodes[node_name], input, value)

    return property(getter, setter)
