"""Socket type definitions for node group interfaces.

These dataclasses define the properties for node group input/output sockets.
Each socket type provides full IDE autocomplete and type checking.

Socket classes are prefixed with 'Socket' to distinguish them from node classes.
For example: SocketVector (interface socket) vs Vector (input node).
"""

from __future__ import annotations
from dataclasses import dataclass


@dataclass
class SocketBase:
    """Base class for all socket definitions."""

    name: str
    bl_socket_type: str = ""
    description: str = ""


@dataclass
class SocketGeometry(SocketBase):
    """Geometry socket - holds mesh, curve, point cloud, or volume data."""

    bl_socket_type: str = "NodeSocketGeometry"


@dataclass
class SocketBoolean(SocketBase):
    """Boolean socket - true/false value."""

    default: bool = False
    bl_socket_type: str = "NodeSocketBool"


@dataclass
class SocketVector(SocketBase):
    """Vector socket - 3D vector (X, Y, Z)."""

    default: tuple[float, float, float] = (0.0, 0.0, 0.0)
    min_value: tuple[float, float, float] | None = None
    max_value: tuple[float, float, float] | None = None
    bl_socket_type: str = "NodeSocketVector"


@dataclass
class SocketFloat(SocketBase):
    """Float socket - single floating point value."""

    default: float = 0.0
    min_value: float | None = None
    max_value: float | None = None
    bl_socket_type: str = "NodeSocketFloat"


@dataclass
class SocketInt(SocketBase):
    """Integer socket - whole number value."""

    default: int = 0
    min_value: int | None = None
    max_value: int | None = None
    bl_socket_type: str = "NodeSocketInt"


@dataclass
class SocketString(SocketBase):
    """String socket - text value."""

    default: str = ""
    bl_socket_type: str = "NodeSocketString"


@dataclass
class SocketColor(SocketBase):
    """Color socket - RGBA color value."""

    default: tuple[float, float, float, float] = (1.0, 1.0, 1.0, 1.0)
    bl_socket_type: str = "NodeSocketColor"


@dataclass
class SocketMaterial(SocketBase):
    """Material socket - Blender material reference."""

    bl_socket_type: str = "NodeSocketMaterial"


@dataclass
class SocketImage(SocketBase):
    """Image socket - Blender image datablock reference."""

    bl_socket_type: str = "NodeSocketImage"


@dataclass
class SocketObject(SocketBase):
    """Object socket - Blender object reference."""

    bl_socket_type: str = "NodeSocketObject"


@dataclass
class SocketCollection(SocketBase):
    """Collection socket - Blender collection reference."""

    bl_socket_type: str = "NodeSocketCollection"


@dataclass
class SocketRotation(SocketBase):
    """Rotation socket - rotation value (Euler or Quaternion)."""

    default: tuple[float, float, float, float] = (1.0, 0.0, 0.0, 0.0)  # Quaternion
    bl_socket_type: str = "NodeSocketRotation"


@dataclass
class SocketMatrix(SocketBase):
    """Matrix socket - 4x4 transformation matrix."""

    bl_socket_type: str = "NodeSocketMatrix"


@dataclass
class SocketTexture(SocketBase):
    """Texture socket - Blender texture reference."""

    bl_socket_type: str = "NodeSocketTexture"
