from typing import Type, TypeVar
import bpy
import numpy as np
from mathutils import Vector
from .arrange import arrange_tree

NODE_SPACING = 250


class TreeInterface:
    def __init__(self):
        self._dynamic_properties = set()
        self.tree: bpy.types.NodeTree
        self._allowable_properties = {"tree", "remove"}

    def __getattr__(self, name):
        # This method is only called when the attribute doesn't exist
        # No need to check if it exists in _dynamic_properties
        raise AttributeError(
            f"'{self.__class__.__name__}' has no attribute '{name}'.\n"
            f"Potentially a dynamic property in {self._dynamic_properties}"
        )

    def __setattr__(self, name, value):
        # For setting attributes, we need to check if it's allowed
        if (
            name != "_dynamic_properties"
            and not name.startswith("_")
            and hasattr(self, "_dynamic_properties")
            and name not in self._dynamic_properties
            and name not in self._allowable_properties
        ):
            raise AttributeError(
                f"Cannot set non-existent property '{name}' on '{self.__class__.__name__}' "
                f"\nAvailable dynamic properties: {self._dynamic_properties}"
            )
        # Use the normal attribute setting mechanism
        super().__setattr__(name, value)

    def _register_property(self, name):
        """Register a dynamic property that can be set on this interface"""
        self._dynamic_properties.add(name)


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
    def __init__(
        self, message="Socket is linked, default value has no bearing on the result"
    ):
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
        return self.tree.nodes[node_name].inputs[socket_name].default_value

    def setter(self, value: T) -> None:
        self.tree.nodes[node_name].inputs[socket_name].default_value = value

    return property(getter, setter)


def check_linked(socket: bpy.types.NodeSocket):
    if socket.is_linked:
        raise SocketLinkedError()


def remove_linked(socket: bpy.types.NodeSocket):
    if socket.is_linked:
        socket.node.id_data.links.remove(socket.links[0])
        arrange_tree(socket.node.id_data)


def create_socket_property(
    socket: bpy.types.NodeSocket,
    value_type: Type,
    data_type: str | None = None,
    length_check: int | None = None,
):
    """
    Generic function to create socket properties with consistent behavior.

    Parameters:
    - socket: The socket to create a property for
    - value_type: The expected Python type(s) for the value
    - data_type: The Blender data type for named attributes
    - length_check: If not None, validates that array-like values have this length
    """

    def getter(self):
        check_linked(socket)
        return getattr(socket, "default_value")

    def setter(self, value):
        remove_linked(socket)

        if isinstance(value, str):
            input_named_attribute(socket, value, data_type)
        elif isinstance(value, value_type):
            # For array-like types, check length if specified
            if length_check is not None and hasattr(value, "__len__"):
                if len(value) != length_check:
                    raise ValueError(
                        f"Value must have length {length_check}, not {len(value)}"
                    )
            setattr(socket, "default_value", value)
        else:
            raise ValueError(
                f"Value must be a {value_type.__name__}, not a {type(value)}"
            )

    return property(getter, setter, doc=socket.description)


def getset_float(socket: bpy.types.NodeSocketFloat):
    return create_socket_property(socket, (float, int), "FLOAT")


def getset_int(socket: bpy.types.NodeSocketInt):
    return create_socket_property(socket, int, "INT")


def getset_bool(socket: bpy.types.NodeSocketBool):
    return create_socket_property(socket, bool, "BOOLEAN")


def getset_vector(socket: bpy.types.NodeSocket):
    return create_socket_property(socket, (np.ndarray, list, tuple), "FLOAT_VECTOR", 3)


def getset_menu(socket: bpy.types.NodeSocketMenu):
    def getter(self) -> str:
        check_linked(socket)
        return getattr(socket, "default_value")

    def setter(self, value: str):
        if isinstance(value, str):
            remove_linked(socket)
            # try first to use title case so we can supply lowercase values
            # but if that fails then try and fail with the orignal value
            try:
                setattr(socket, "default_value", value.title())
            except TypeError:
                setattr(socket, "default_value", value)
        else:
            raise ValueError(f"Value must be a string, not a {type(value)}")

    return property(getter, setter, doc=socket.description)


def getset_string(socket: bpy.types.NodeSocket):
    return create_socket_property(socket, str, "STRING")


def getset_color(socket: bpy.types.NodeSocketColor):
    return create_socket_property(socket, (np.ndarray, list, tuple), "FLOAT_VECTOR", 4)


def getset_rotation(socket: bpy.types.NodeSocketRotation):
    return create_socket_property(socket, (np.ndarray, list, tuple), "QUATERNION", 4)


def socket(socket: bpy.types.NodeSocket):
    match socket.type:
        case "FLOAT" | "VALUE":
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
        case "MENU":
            return getset_menu(socket)
        case _:
            raise ValueError(f"Unknown socket type: {socket.type}")


def option(node_name: str, input: str, type_: type[T]):
    def getter(self) -> T:
        return self.tree.nodes[node_name][input]

    def setter(self, value: T) -> None:
        setattr(self.tree.nodes[node_name], input, value)

    return property(getter, setter)
