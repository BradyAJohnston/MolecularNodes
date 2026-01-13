"""
Some of the nodes need to be manually specified because they are a bit tricky to generate automatically.
"""

from __future__ import annotations
import bpy
from typing_extensions import Literal
from ..builder import NodeBuilder, NodeSocket
from . import types
from .types import (
    LINKABLE,
    TYPE_INPUT_VECTOR,
    DataTypes,
)


class SocketAccessor:
    def __init__(self, node_builder: NodeBuilder):
        self._builder = node_builder

    @property
    def _inputs(self):
        return self._builder.node.inputs

    @property
    def _outputs(self):
        return self._builder.node.outputs


class RandomValue(NodeBuilder):
    """Random Value node"""

    name = "FunctionNodeRandomValue"

    _default_input_id = "ID"

    def __init__(
        self,
        data_type: Literal["FLOAT", "INT", "BOOL", "FLOAT_VECTOR"],
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
    def value(self) -> NodeSocket:
        """Output socket: Value"""
        return self._default_output_socket

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
        builder = cls(Probability=probability, id=id, seed=seed, data_type="BOOL")
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
    def inputs(self) -> SocketAccessor:
        class SeparateXYZInputsAccessor(SocketAccessor):
            _builder: SeparateXYZ

            @property
            def vector(self) -> NodeSocket:
                """Input socket: Vector"""
                return self._inputs["Vector"]

        return SeparateXYZInputsAccessor(self)

    @property
    def outputs(self) -> SocketAccessor:
        class SeparateXYZOutputsAccessor(SocketAccessor):
            _builder: SeparateXYZ

            @property
            def x(self) -> NodeSocket:
                """Output socket: X"""
                return self._outputs["X"]

            @property
            def y(self) -> NodeSocket:
                """Output socket: Y"""
                return self._outputs["Y"]

            @property
            def z(self) -> NodeSocket:
                """Output socket: Z"""
                return self._outputs["Z"]

        return SeparateXYZOutputsAccessor(self)

    @property
    def x(self) -> NodeSocket:
        """Output socket: X"""
        return self.node.outputs["X"]

    @property
    def y(self) -> NodeSocket:
        """Output socket: Y"""
        return self.node.outputs["Y"]

    @property
    def z(self) -> NodeSocket:
        """Output socket: Z"""
        return self.node.outputs["Z"]


class Mix(NodeBuilder):
    """Mix values by a factor"""

    name = "ShaderNodeMix"
    node: bpy.types.ShaderNodeMix  # type: ignore

    def __init__(
        self,
        data_type: Literal["FLOAT", "VECTOR", "COLOR", "ROTATION"] = "FLOAT",
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
    def data_type(self, value: Literal["FLOAT", "VECTOR", "RGBA", "ROTATION"]):
        self.node.data_type = value

    @property
    def factor_mode(self) -> str:
        return self.node.factor_mode

    @factor_mode.setter
    def factor_mode(self, value: Literal["NON_UNIFORM", "UNIFORM"]):
        self.node.factor_mode = value

    @property
    def outputs(self) -> SocketAccessor:
        class MixOutputAccessor(SocketAccessor):
            _builder: Mix

            @property
            def result(self) -> NodeSocket:
                """Output socket: Result"""
                return self._builder._default_output_socket

        return MixOutputAccessor(self)

    @property
    def inputs(self) -> SocketAccessor:
        class MixInputAccessor(SocketAccessor):
            _builder: Mix

            @property
            def factor(self) -> NodeSocket:
                """Input socket: Factor"""
                match self._builder.data_type:
                    case "FLOAT":
                        name = "Factor_Float"
                    case "VECTOR":
                        name = (
                            "Factor_Vector"
                            if self._builder.factor_mode == "NON_UNIFORM"
                            else "Factor_Float"
                        )
                    case "RGBA":
                        name = "Factor_Color"
                    case "ROTATION":
                        name = "Factor_Rotation"
                    case _:
                        raise ValueError(
                            f"Unsupported data type: {self._builder.data_type}"
                        )

                name = self._builder._input_idx(name)
                return self._builder.node.inputs[f"Factor_{name}"]

        return MixInputAccessor(self)

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
            data_type=DataTypes.FLOAT,
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
                    data_type=DataTypes.VECTOR,
                )
            case "NON_UNIFORM":
                builder = cls(
                    Factor_Vector=factor,
                    A_Vector=a,
                    B_Vector=b,
                    data_type=DataTypes.VECTOR,
                )
            case _:
                raise ValueError(f"Invalid factor mode: {factor_mode}")

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
            data_type=DataTypes.COLOR,
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
            data_type=DataTypes.ROTATION,
        )
        builder.node.clamp_factor = clamp_factor
        return builder


class Math(NodeBuilder):
    """Perform math operations"""

    name = "ShaderNodeMath"
    node: bpy.types.ShaderNodeMath  # type: ignore

    def __init__(
        self, operation: types.NodeMathItems, use_clamp: bool | None = None, **kwargs
    ):
        super().__init__()
        if operation is not None:
            self.node.operation = operation
        if use_clamp is not None:
            self.node.use_clamp = use_clamp
        self._establish_links(**kwargs)

    @property
    def outputs(self) -> SocketAccessor:
        class MathOutputs(SocketAccessor):
            @property
            def value(self):
                return self._outputs["Value"]

        return MathOutputs(self)

    @property
    def inputs(self) -> SocketAccessor:
        class MathInputs(SocketAccessor):
            @property
            def value(self):
                return self._inputs["Value"]

            @property
            def value_001(self):
                return self._inputs["Value_001"]

            @property
            def value_002(self):
                return self._inputs["Value_002"]

        return MathInputs(self)

    @property
    def use_clamp(self):
        return self.node.use_clamp

    @use_clamp.setter
    def use_clamp(self, value: bool):
        self.node.use_clamp = value

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
