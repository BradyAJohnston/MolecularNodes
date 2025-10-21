"""
Socket wrapper providing fluent linking API for Blender node sockets.
"""

from typing import TYPE_CHECKING, Union

import bpy

if TYPE_CHECKING:
    from .node import NodeWrapper


class SocketWrapper:
    """
    Wrapper around bpy.types.NodeSocket providing fluent linking capabilities.

    Examples:
        # Link sockets together
        output.link_to(input)

        # Operator overload syntax
        input << output

        # Chain links
        output.link_to(node1.input(0)).link_to(node2.input(0))
    """

    def __init__(self, socket: bpy.types.NodeSocket):
        """
        Initialize socket wrapper.

        Args:
            socket: The Blender socket to wrap
        """
        self._socket = socket

    @property
    def socket(self) -> bpy.types.NodeSocket:
        """Get the underlying Blender socket."""
        return self._socket

    @property
    def node(self) -> bpy.types.Node:
        """Get the node this socket belongs to."""
        return self._socket.node

    @property
    def tree(self) -> bpy.types.NodeTree:
        """Get the node tree this socket belongs to."""
        return self._socket.node.id_data

    @property
    def name(self) -> str:
        """Get the socket name."""
        return self._socket.name

    @property
    def identifier(self) -> str:
        """Get the socket identifier."""
        return self._socket.identifier

    @property
    def is_output(self) -> bool:
        """Check if this is an output socket."""
        return self._socket.is_output

    @property
    def is_linked(self) -> bool:
        """Check if this socket is linked."""
        return self._socket.is_linked

    @property
    def default_value(self):
        """Get/set the socket's default value."""
        return self._socket.default_value

    @default_value.setter
    def default_value(self, value):
        self._socket.default_value = value

    def link_to(self, target: Union["SocketWrapper", bpy.types.NodeSocket]) -> "SocketWrapper":
        """
        Create a link from this socket to a target socket.

        Args:
            target: Target socket to link to (can be SocketWrapper or raw socket)

        Returns:
            The target socket wrapper for chaining

        Raises:
            ValueError: If trying to link incompatible socket directions
        """
        # Unwrap if needed
        target_socket = target.socket if isinstance(target, SocketWrapper) else target

        # Determine link direction
        if self.is_output and not target_socket.is_output:
            # Output -> Input
            self.tree.links.new(target_socket, self._socket)
        elif not self.is_output and target_socket.is_output:
            # Input <- Output
            self.tree.links.new(self._socket, target_socket)
        else:
            raise ValueError(
                f"Cannot link {self.name} ({'output' if self.is_output else 'input'}) "
                f"to {target_socket.name} ({'output' if target_socket.is_output else 'input'}). "
                "Must link output to input or input to output."
            )

        # Return target wrapper for chaining
        return target if isinstance(target, SocketWrapper) else SocketWrapper(target_socket)

    def __lshift__(self, other: Union["SocketWrapper", bpy.types.NodeSocket]) -> "SocketWrapper":
        """
        Operator overload for linking: input << output

        Args:
            other: Output socket to link from

        Returns:
            Self for chaining
        """
        other_wrapped = other if isinstance(other, SocketWrapper) else SocketWrapper(other)
        other_wrapped.link_to(self)
        return self

    def unlink(self) -> "SocketWrapper":
        """
        Remove all links from this socket.

        Returns:
            Self for chaining
        """
        if self.is_linked:
            for link in list(self._socket.links):
                self.tree.links.remove(link)
        return self

    def __repr__(self) -> str:
        direction = "output" if self.is_output else "input"
        return f"<SocketWrapper {direction}:{self.name} on {self.node.name}>"


class SocketCollection:
    """
    Collection of sockets providing dict-like and attribute access.

    Examples:
        # Access by name
        sockets.Geometry
        sockets["Geometry"]

        # Access by index
        sockets[0]
    """

    def __init__(self, sockets: bpy.types.NodeOutputs | bpy.types.NodeInputs):
        """
        Initialize socket collection.

        Args:
            sockets: Blender socket collection (inputs or outputs)
        """
        self._sockets = sockets

    def __getattr__(self, name: str) -> SocketWrapper:
        """Get socket by name as attribute."""
        if name.startswith("_"):
            # Don't interfere with private attributes
            return object.__getattribute__(self, name)

        try:
            return SocketWrapper(self._sockets[name])
        except KeyError:
            raise AttributeError(f"No socket named '{name}'")

    def __getitem__(self, key: Union[str, int]) -> SocketWrapper:
        """Get socket by name or index."""
        return SocketWrapper(self._sockets[key])

    def __len__(self) -> int:
        """Get number of sockets."""
        return len(self._sockets)

    def __iter__(self):
        """Iterate over socket wrappers."""
        for socket in self._sockets:
            yield SocketWrapper(socket)

    def __repr__(self) -> str:
        socket_names = [s.name for s in self._sockets]
        return f"<SocketCollection {socket_names}>"
