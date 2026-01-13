"""
Some of the nodes need to be manually specified because they are a bit tricky to generate automatically.
"""

from __future__ import annotations
import bpy
from bpy.types import NodeSocketFloat
from typing_extensions import Literal
from ..builder import NodeBuilder, NodeSocket
from . import types
from .types import LINKABLE, TYPE_INPUT_BOOLEAN, TYPE_INPUT_VECTOR

_RANDOM_VALUE_DATA_TYPES = Literal["FLOAT", "INT", "BOOLEAN", "FLOAT_VECTOR"]


class RandomValue(NodeBuilder):
    """Random Value node"""

    name = "FunctionNodeRandomValue"
    node: bpy.types.FunctionNodeRandomValue
    _default_input_id = "ID"

    def __init__(
        self,
        data_type: _RANDOM_VALUE_DATA_TYPES,
        id: int | LINKABLE | None = None,
        seed: int | LINKABLE | None = None,
        **kwargs,
    ):
        super().__init__()
        self.node.data_type = data_type
        key_args = {
            "ID": id,
            "Seed": seed,
        }
        key_args.update(kwargs)

        self._establish_links(**key_args)

    @property
    def data_type(self) -> _RANDOM_VALUE_DATA_TYPES:
        return self.node.data_type  # type: ignore

    @data_type.setter
    def data_type(self, value: _RANDOM_VALUE_DATA_TYPES):
        self.node.data_type = value

    @property
    def o_value(self) -> NodeSocket:
        """Output socket: Value"""
        match self.data_type:
            case "FLOAT":
                return self._output("Value_001")
            case "INT":
                return self._output("Value_002")
            case "BOOLEAN":
                return self._output("Value_003")
            case "FLOAT_VECTOR":
                return self._output("Value")

    def i_min(self) -> NodeSocket:
        """Input socket: Minimum"""
        match self.data_type:
            case "FLOAT":
                return self._input("Min_001")
            case "INT":
                return self._input("Min_002")
            case "BOOLEAN":
                raise ValueError(
                    f"Boolean data type does not support minimum value, use 'Probability'"
                )
            case "FLOAT_VECTOR":
                return self._input("Min")

    def i_max(self) -> NodeSocket:
        """Input socket: Maximum"""
        match self.data_type:
            case "FLOAT":
                return self._input("Max_001")
            case "INT":
                return self._input("Max_002")
            case "BOOLEAN":
                raise ValueError(
                    f"Boolean data type does not support maximum value, use 'Probability'"
                )
            case "FLOAT_VECTOR":
                return self._input("Max")

    def i_probability(self) -> NodeSocket:
        """Input socket: Probability"""
        if self.data_type != "BOOLEAN":
            raise ValueError(
                f"Probability socket is only supported for boolean data types, not for data type: {self.data_type}"
            )

        return self._input("Probability")

    @classmethod
    def float(
        cls,
        min: float | LINKABLE = 0.0,
        max: float | LINKABLE = 1.0,
        id: int | LINKABLE | None = None,
        seed: int | LINKABLE = 1,
    ) -> NodeBuilder:
        buidler = cls(Min_001=min, Max_001=max, id=id, seed=seed, data_type="FLOAT")
        buidler._default_output_id = "Value_001"
        return buidler

    @classmethod
    def integer(
        cls,
        min: int | LINKABLE = 0,
        max: int | LINKABLE = 1,
        id: int | LINKABLE | None = None,
        seed: int | LINKABLE = 1,
    ) -> NodeBuilder:
        buidler = cls(Min_002=min, Max_002=max, id=id, seed=seed, data_type="INT")
        buidler._default_output_id = "Value_002"
        return buidler

    @classmethod
    def boolean(
        cls,
        probability: float | LINKABLE = 0.5,
        id: int | LINKABLE | None = None,
        seed: int | LINKABLE = 1,
    ) -> NodeBuilder:
        builder = cls(Probability=probability, id=id, seed=seed, data_type="BOOLEAN")
        builder._default_output_id = "Value_003"
        return builder

    @classmethod
    def vector(
        cls,
        min: tuple[float, float, float] | list[float] | LINKABLE = (0.0, 0.0, 0.0),
        max: tuple[float, float, float] | list[float] | LINKABLE = (1.0, 1.0, 1.0),
        id: int | LINKABLE | None = None,
        seed: int | LINKABLE = 1,
    ) -> NodeBuilder:
        buidler = cls(Min=min, Max=max, id=id, seed=seed, data_type="FLOAT_VECTOR")
        buidler._default_output_id = "Value"
        return buidler


class SeparateXYZ(NodeBuilder):
    """Split a vector into its X, Y, and Z components"""

    name = "ShaderNodeSeparateXYZ"
    node: bpy.types.ShaderNodeSeparateXYZ  # type: ignore

    def __init__(self, vector: TYPE_INPUT_VECTOR | None = None):
        super().__init__()
        self._establish_links(**{"Vector": vector})

    @property
    def i_vector(self) -> NodeSocket:
        """Input socket: Vector"""
        return self._input("Vector")

    @property
    def o_x(self) -> NodeSocket:
        """Output socket: X"""
        return self._output("X")

    @property
    def o_y(self) -> NodeSocket:
        """Output socket: Y"""
        return self._output("Y")

    @property
    def o_z(self) -> NodeSocket:
        """Output socket: Z"""
        return self._output("Z")


class CombineXYZ(NodeBuilder):
    """Create a vector from X, Y, and Z components"""

    name = "ShaderNodeCombineXYZ"
    node: bpy.types.ShaderNodeCombineXYZ  # type: ignore

    def __init__(
        self,
        x: float | LINKABLE = 0.0,
        y: float | LINKABLE = 0.0,
        z: float | LINKABLE = 0.0,
    ):
        super().__init__()
        self._establish_links(**{"X": x, "Y": y, "Z": z})

    @property
    def o_vector(self) -> NodeSocket:
        """Output socket: Vector"""
        return self._output("Vector")

    @property
    def i_x(self) -> NodeSocket:
        """Input socket: X"""
        return self._input("X")

    @property
    def i_y(self) -> NodeSocket:
        """Input socket: Y"""
        return self._input("Y")

    @property
    def i_z(self) -> NodeSocket:
        """Input socket: Z"""
        return self._input("Z")


_MIX_VALUE_DATA_TYPES = Literal["FLOAT", "VECTOR", "COLOR", "ROTATION"]


class Mix(NodeBuilder):
    """Mix values by a factor"""

    name = "ShaderNodeMix"
    node: bpy.types.ShaderNodeMix  # type: ignore

    def __init__(
        self,
        data_type: _MIX_VALUE_DATA_TYPES = "FLOAT",
        **kwargs,
    ):
        super().__init__()
        self._default_input_id = f"A_{data_type.title()}"
        self._default_output_id = f"Result_{data_type.title()}"
        self.node.data_type = "RGBA" if data_type == "COLOR" else data_type
        key_args = {}
        key_args.update(kwargs)
        self._establish_links(**key_args)

    @property
    def data_type(self) -> str:
        return self.node.data_type

    @data_type.setter
    def data_type(self, value: _MIX_VALUE_DATA_TYPES):
        self.node.data_type = value  # type: ignore

    @property
    def factor_mode(self) -> Literal["UNIFORM", "NON_UNIFORM"]:
        return self.node.factor_mode

    @factor_mode.setter
    def factor_mode(self, value: Literal["NON_UNIFORM", "UNIFORM"]):
        self.node.factor_mode = value

    @property
    def o_result(self) -> NodeSocket:
        """Output socket: Result"""
        return self._default_output_socket

    @property
    def i_factor(self) -> NodeSocket:
        """Input socket: Factor"""
        match self.data_type:
            case "FLOAT":
                name = "Factor_Float"
            case "VECTOR":
                name = (
                    "Factor_Float" if self.factor_mode == "UNIFORM" else "Factor_Vector"
                )
            case "RGBA":
                name = "Factor_Color"
            case "ROTATION":
                name = "Factor_Rotation"
            case _:
                raise ValueError(f"Unsupported data type: {self.data_type}")

        idx = self._input_idx(name)
        return self.node.inputs[idx]

    @property
    def i_value_a(self) -> NodeSocket:
        """Input socket: Value A"""
        type_name = "Color" if self.data_type == "RGBA" else self.data_type
        name = f"A_{type_name}"
        idx = self._input_idx(name)
        return self.node.inputs[idx]

    @property
    def i_value_b(self) -> NodeSocket:
        """Input socket: Value B"""
        type_name = "Color" if self.data_type == "RGBA" else self.data_type
        name = f"B_{type_name}"
        idx = self._input_idx(name)
        return self.node.inputs[idx]

    @classmethod
    def float(
        cls,
        factor: float | LINKABLE = 0.5,
        a: float | LINKABLE = 0.0,
        b: float | LINKABLE = 0.0,
        clamp_factor: bool | LINKABLE = True,
    ) -> Mix:
        builder = cls(
            Factor_Float=factor,
            A_Float=a,
            B_Float=b,
            data_type="COLOR",
        )
        builder.node.clamp_factor = clamp_factor
        return builder

    @classmethod
    def vector(
        cls,
        factor: float | LINKABLE = 0.5,
        a: tuple[float, float, float] | list[float] | LINKABLE = (0.0, 0.0, 0.0),
        b: tuple[float, float, float] | list[float] | LINKABLE = (1.0, 1.0, 1.0),
        clamp_factor: bool = True,
        factor_mode: Literal["UNIFORM", "NON_UNIFORM"] = "UNIFORM",
    ) -> Mix:
        match factor_mode:
            case "UNIFORM":
                builder = cls(
                    Factor_Float=factor,
                    A_Vector=a,
                    B_Vector=b,
                    data_type="VECTOR",
                )
            case "NON_UNIFORM":
                builder = cls(
                    Factor_Vector=factor,
                    A_Vector=a,
                    B_Vector=b,
                    data_type="VECTOR",
                )

        builder.node.clamp_factor = clamp_factor
        return builder

    @classmethod
    def color(
        cls,
        factor: float | LINKABLE = 0.5,
        a: tuple[float, float, float] | list[float] | LINKABLE = (0.0, 0.0, 0.0),
        b: tuple[float, float, float] | list[float] | LINKABLE = (1.0, 1.0, 1.0),
        blend_type: str = "add",
        clamp_factor: bool = True,
        clamp_result: bool = True,
    ) -> Mix:
        builder = cls(
            Factor_Float=factor,
            A_Color=a,
            B_Color=b,
            data_type="COLOR",
        )
        builder.node.blend_type = blend_type.capitalize()
        builder.node.clamp_factor = clamp_factor
        builder.node.clamp_result = clamp_result
        return builder

    @classmethod
    def rotation(
        cls,
        a: tuple[float, float, float, float] | list[float] | LINKABLE | None = None,
        b: tuple[float, float, float, float] | list[float] | LINKABLE | None = None,
        factor: float | LINKABLE = 0.5,
        clamp_factor: bool = True,
    ) -> Mix:
        builder = cls(
            Factor_Float=factor,
            A_Rotation=a,
            B_Rotation=b,
            data_type="ROTATION",
        )
        builder.node.clamp_factor = clamp_factor
        return builder


class Math(NodeBuilder):
    """Perform math operations"""

    name = "ShaderNodeMath"
    node: bpy.types.ShaderNodeMath  # type: ignore

    def __init__(
        self,
        operation: types.NodeMathItems = "ADD",
        use_clamp: bool = False,
        **kwargs,
    ):
        super().__init__()
        self.operation = operation
        self.use_clamp = use_clamp
        self._establish_links(**kwargs)

    @property
    def operation(self) -> types.NodeMathItems:
        return self.node.operation

    @operation.setter
    def operation(self, value: types.NodeMathItems):
        self.node.operation = value

    @property
    def use_clamp(self) -> bool:
        return self.node.use_clamp

    @use_clamp.setter
    def use_clamp(self, value: bool):
        self.node.use_clamp = value

    @property
    def o_value(self) -> NodeSocketFloat:
        return self._output("Value")  # type: ignore

    def _input(self, identifier: str) -> NodeSocketFloat:
        return self._input(identifier)

    @property
    def i_value(self) -> NodeSocketFloat:
        return self._input("Value")

    @property
    def i_value_001(self) -> NodeSocketFloat:
        return self._input("Value_001")

    @property
    def i_value_002(self) -> NodeSocketFloat:
        return self._input("Value_002")

    @classmethod
    def add(
        cls,
        a: float | LINKABLE = 0.5,
        b: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation of `a + b`."""
        return cls(operation="ADD", Value=a, Value_001=b)

    @classmethod
    def subtract(
        cls,
        a: float | LINKABLE = 0.5,
        b: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation of `a - b`."""
        return cls(operation="SUBTRACT", Value=a, Value_001=b)

    @classmethod
    def multiply(
        cls,
        a: float | LINKABLE = 0.5,
        b: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation of `a * b`."""
        return cls(operation="MULTIPLY", Value=a, Value_001=b)

    @classmethod
    def divide(
        cls,
        a: float | LINKABLE = 0.5,
        b: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation of `a / b`."""
        return cls(operation="DIVIDE", Value=a, Value_001=b)

    @classmethod
    def multiply_add(
        cls,
        a: float | LINKABLE = 0.5,
        b: float | LINKABLE = 0.5,
        c: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `a * b + c`."""
        return cls(operation="MULTIPLY_ADD", Value=a, Value_001=b, value_002=c)

    @classmethod
    def power(
        cls,
        base: float | LINKABLE = 0.5,
        exponent: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `base ** exponent`."""
        return cls(operation="POWER", Value=base, Value_001=exponent)

    @classmethod
    def logarithm(
        cls,
        value: float | LINKABLE = 0.5,
        base: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `log(value, base)`."""
        return cls(operation="LOGARITHM", Value=value, Value_001=base)

    @classmethod
    def sqrt(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `sqrt(value)`."""
        return cls(operation="SQRT", Value=value)

    @classmethod
    def inverse_sqrt(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `inverse_sqrt(value)`."""
        return cls(operation="INVERSE_SQRT", Value=value)

    @classmethod
    def absolute(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `abs(value)`."""
        return cls(operation="ABSOLUTE", Value=value)

    @classmethod
    def exponent(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `exp(value)`."""
        return cls(operation="EXPONENT", Value=value)

    @classmethod
    def minimum(
        cls,
        a: float | LINKABLE = 0.5,
        b: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `min(a, b)`."""
        return cls(operation="MINIMUM", Value=a, Value_001=b)

    @classmethod
    def maximum(
        cls,
        a: float | LINKABLE = 0.5,
        b: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `max(a, b)`."""
        return cls(operation="MAXIMUM", Value=a, Value_001=b)

    @classmethod
    def less_than(
        cls,
        value: float | LINKABLE = 0.5,
        threshold: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `value < threshold` returning 1 or 0."""
        return cls(operation="LESS_THAN", Value=value, Value_001=threshold)

    @classmethod
    def greater_than(
        cls,
        value: float | LINKABLE = 0.5,
        threshold: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `value > threshold` returning 1 or 0."""
        return cls(operation="GREATER_THAN", Value=value, Value_001=threshold)

    @classmethod
    def sign(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `sign(value)` returning -1, 0, or 1."""
        return cls(operation="SIGN", Value=value)

    @classmethod
    def compare(
        cls,
        a: float | LINKABLE = 0.5,
        b: float | LINKABLE = 0.5,
        epsilon: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `compare(a, b, epsilon)` returning -1, 0, or 1."""
        return cls(operation="COMPARE", Value=a, Value_001=b, value_002=epsilon)

    @classmethod
    def smooth_min(
        cls,
        a: float | LINKABLE = 0.5,
        b: float | LINKABLE = 0.5,
        distance: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `smooth_min(a, b, distance)`."""
        return cls(operation="SMOOTH_MIN", Value=a, Value_001=b, value_002=distance)

    @classmethod
    def smooth_max_(
        cls,
        a: float | LINKABLE = 0.5,
        b: float | LINKABLE = 0.5,
        distance: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `smooth_max(a, b, distance)`."""
        return cls(operation="SMOOTH_MAX", Value=a, Value_001=b, value_002=distance)

    @classmethod
    def round(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Round A to the nearest integer. Round up if 0.5 or greater."""
        return cls(operation="ROUND", Value=value)

    @classmethod
    def floor(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """The largest integer smaller than or equal to `value`"""
        return cls(operation="FLOOR", Value=value)

    @classmethod
    def ceil(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """The smallest integer greater than or equal to `value`"""
        return cls(operation="CEIL", Value=value)

    @classmethod
    def truncate(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """The integer part of `value` removing the fractional part"""
        return cls(operation="TRUNC", Value=value)

    @classmethod
    def fraction(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """The fractional part of `value`"""
        return cls(operation="FRACT", Value=value)

    @classmethod
    def truncated_modulo(
        cls,
        a: float | LINKABLE = 0.5,
        b: float | LINKABLE = 0.5,
    ) -> "Math":
        """The remained of truncated division using fmod(a, b)"""
        return cls(operation="MODULO", Value=a, Value_001=b)

    @classmethod
    def floored_modulo(
        cls,
        a: float | LINKABLE = 0.5,
        b: float | LINKABLE = 0.5,
    ) -> "Math":
        """The remained of floored division"""
        return cls(operation="FLOORED_MODULO", Value=a, Value_001=b)

    @classmethod
    def wrap(
        cls,
        value: float | LINKABLE = 0.5,
        max: float | LINKABLE = 0.5,
        min: float | LINKABLE = 0.5,
    ) -> "Math":
        """Wrap value to range, wrap(value, max, min)"""
        return cls(operation="WRAP", Value=value, Value_001=max, value_002=min)

    @classmethod
    def snap(
        cls,
        value: float | LINKABLE = 0.5,
        increment: float | LINKABLE = 0.5,
    ) -> "Math":
        """Snap to increment of `snap(value, increment)`"""
        return cls(operation="SNAP", Value=value, Value_001=increment)

    @classmethod
    def ping_pong(
        cls,
        value: float | LINKABLE = 0.5,
        scale: float | LINKABLE = 0.5,
    ) -> "Math":
        """Wraps a value and reverses every other cycle"""
        return cls(operation="PINGPONG", Value=value, Value_001=scale)

    @classmethod
    def sine(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation 'sin(value)'."""
        return cls(operation="SINE", Value=value)

    @classmethod
    def cosine(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation 'cos(value)'."""
        return cls(operation="COSINE", Value=value)

    @classmethod
    def tangent(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation 'tan(value)'."""
        return cls(operation="TANGENT", Value=value)

    @classmethod
    def arcsine(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `arcsin(value)'."""
        return cls(operation="ARCSINE", Value=value)

    @classmethod
    def arccosine(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation 'arccos(value)'."""
        return cls(operation="ARCCOSINE", Value=value)

    @classmethod
    def arctangent(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation 'arctan(value)'."""
        return cls(operation="ARCTANGENT", Value=value)

    @classmethod
    def arctan2(
        cls,
        a: float | LINKABLE = 0.5,
        b: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation 'arctan(a / b)'."""
        return cls(operation="ARCTAN2", Value=a, Value_001=b)

    @classmethod
    def sinh(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `sinh(value)`."""
        return cls(operation="SINH", Value=value)

    @classmethod
    def cosh(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `cosh(value)`."""
        return cls(operation="COSH", Value=value)

    @classmethod
    def tanh(
        cls,
        value: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `tanh(value)`."""
        return cls(operation="TANH", Value=value)

    @classmethod
    def radians(
        cls,
        degrees: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation `radians(degrees)`."""
        return cls(operation="RADIANS", Value=degrees)

    @classmethod
    def degrees(
        cls,
        radians: float | LINKABLE = 0.5,
    ) -> "Math":
        """Create Math with operation 'To Degrees'."""
        return cls(operation="DEGREES", Value=radians)


class BooleanMath(NodeBuilder):
    """Boolean Math node"""

    name = "FunctionNodeBooleanMath"
    node: bpy.types.FunctionNodeBooleanMath

    def __init__(self, operation: types.NodeBooleanMathItems = "AND", **kwargs):
        super().__init__()
        self.operator = operation
        self._establish_links(**kwargs)

    @property
    def operation(self) -> types.NodeBooleanMathItems:
        return self.node.operation

    @operation.setter
    def operation(self, value: types.NodeBooleanMathItems):
        self.node.operation = value

    @property
    def i_boolean(self) -> bpy.types.NodeSocketBool:
        return self._input("Boolean")  # type: ignore

    @property
    def i_boolean_001(self) -> bpy.types.NodeSocketBool:
        return self._input("Boolean_001")  # type: ignore

    @property
    def o_boolean(self) -> bpy.types.NodeSocketBool:
        return self._output("Boolean")  # type: ignore

    @classmethod
    def l_and(
        cls,
        boolean: TYPE_INPUT_BOOLEAN = False,
        boolean_001: TYPE_INPUT_BOOLEAN = False,
    ) -> "BooleanMath":
        """Create Boolean Math with operation 'AND'."""
        return cls(operation="AND", Boolean=boolean, Boolean_001=boolean_001)

    @classmethod
    def l_or(
        cls,
        boolean: TYPE_INPUT_BOOLEAN = False,
        boolean_001: TYPE_INPUT_BOOLEAN = False,
    ) -> "BooleanMath":
        """Create Boolean Math with operation 'OR'."""
        return cls(operation="OR", Boolean=boolean, Boolean_001=boolean_001)

    @classmethod
    def l_not(cls, boolean: TYPE_INPUT_BOOLEAN = False) -> "BooleanMath":
        """Create Boolean Math with operation 'NOT'."""
        return cls(operation="NOT", Boolean=boolean)

    @classmethod
    def l_not_and(
        cls,
        boolean: TYPE_INPUT_BOOLEAN = False,
        boolean_001: TYPE_INPUT_BOOLEAN = False,
    ) -> "BooleanMath":
        """Create Boolean Math with operation 'NAND'."""
        return cls(operation="NAND", Boolean=boolean, Boolean_001=boolean_001)

    @classmethod
    def l_nor(
        cls,
        boolean: TYPE_INPUT_BOOLEAN = False,
        boolean_001: TYPE_INPUT_BOOLEAN = False,
    ) -> "BooleanMath":
        """Create Boolean Math with operation 'NOR'."""
        return cls(operation="NOR", Boolean=boolean, Boolean_001=boolean_001)

    @classmethod
    def l_equal(
        cls,
        boolean: TYPE_INPUT_BOOLEAN = False,
        boolean_001: TYPE_INPUT_BOOLEAN = False,
    ) -> "BooleanMath":
        """Create Boolean Math with operation 'XNOR'."""
        return cls(operation="XNOR", Boolean=boolean, Boolean_001=boolean_001)

    @classmethod
    def l_not_equal(
        cls,
        boolean: TYPE_INPUT_BOOLEAN = False,
        boolean_001: TYPE_INPUT_BOOLEAN = False,
    ) -> "BooleanMath":
        """Create Boolean Math with operation 'XOR'."""
        return cls(operation="XOR", Boolean=boolean, Boolean_001=boolean_001)

    @classmethod
    def l_imply(
        cls,
        boolean: TYPE_INPUT_BOOLEAN = False,
        boolean_001: TYPE_INPUT_BOOLEAN = False,
    ) -> "BooleanMath":
        """Create Boolean Math with operation 'IMPLY'."""
        return cls(operation="IMPLY", Boolean=boolean, Boolean_001=boolean_001)

    @classmethod
    def l_subtract(
        cls,
        boolean: TYPE_INPUT_BOOLEAN = False,
        boolean_001: TYPE_INPUT_BOOLEAN = False,
    ) -> "BooleanMath":
        """Create Boolean Math with operation 'NIMPLY'."""
        return cls(operation="NIMPLY", Boolean=boolean, Boolean_001=boolean_001)
