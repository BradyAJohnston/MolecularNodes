"""
Some of the nodes need to be manually specified because they are a bit tricky to generate automatically.
"""

from __future__ import annotations
from asyncio import sleep
from types import ClassMethodDescriptorType
from mathutils import Vector
from ..builder import NodeBuilder, NodeSocket
from .types import (
    LINKABLE,
    TYPE_INPUT_BOOLEAN,
    TYPE_INPUT_ROTATION,
    TYPE_INPUT_VECTOR,
    DataTypes,
)


class RandomValue(NodeBuilder):
    """Random Value node"""

    name = "FunctionNodeRandomValue"

    _default_input_id = "ID"

    def __init__(
        self,
        data_type: str,
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
    def float_(
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
    def integer_(
        cls,
        min: int | LINKABLE = 0,
        max: int | LINKABLE = 1,
        id: int | LINKABLE | None = None,
        seed: int | LINKABLE = 1,
    ) -> NodeBuilder:
        buidler = cls(
            Min_002=min, Max_002=max, id=id, seed=seed, data_type=DataTypes.INT
        )
        buidler._default_output_id = "Value_002"
        return buidler

    @classmethod
    def boolean_(
        cls,
        probability: float | LINKABLE = 0.5,
        id: int | LINKABLE | None = None,
        seed: int | LINKABLE = 1,
    ) -> NodeBuilder:
        builder = cls(
            Probability=probability, id=id, seed=seed, data_type=DataTypes.BOOL
        )
        builder._default_output_id = "Value_003"
        return builder

    @classmethod
    def vector_(
        cls,
        min: tuple[float, float, float] | list[float] | LINKABLE = (0.0, 0.0, 0.0),
        max: tuple[float, float, float] | list[float] | LINKABLE = (1.0, 1.0, 1.0),
        id: int | LINKABLE | None = None,
        seed: int | LINKABLE = 1,
    ) -> NodeBuilder:
        buidler = cls(Min=min, Max=max, id=id, seed=seed, data_type="FLOAT_VECTOR")
        buidler._default_output_id = "Value"
        return buidler


class Mix(NodeBuilder):
    """Mix values by a factor"""

    name = "ShaderNodeMix"

    def __init__(
        self,
        default_input_id: str,
        default_output_id: str,
        data_type: str | None = None,
        **kwargs,
    ):
        super().__init__()
        self._default_input_id = default_input_id
        self._default_output_id = default_output_id
        self.node.data_type = data_type
        key_args = {}
        key_args.update(kwargs)
        self._establish_links(**key_args)

    @property
    def result(self) -> NodeSocket:
        """Output socket: Result"""
        return self._default_output_socket

    @classmethod
    def float_(
        cls,
        factor: float | LINKABLE = 0.5,
        a: float | LINKABLE = 0.0,
        b: float | LINKABLE = 0.0,
        clamp_factor: bool | LINKABLE = True,
    ) -> Mix:
        builder = cls(
            default_input_id="A_Float",
            default_output_id="Result_Float",
            Factor_Float=factor,
            A_Float=a,
            B_Float=b,
            data_type=DataTypes.FLOAT,
        )
        builder.node.clamp_factor = clamp_factor
        return builder

    @classmethod
    def vector_(
        cls,
        factor: float | LINKABLE = 0.5,
        a: tuple[float, float, float] | list[float] | LINKABLE = (0.0, 0.0, 0.0),
        b: tuple[float, float, float] | list[float] | LINKABLE = (1.0, 1.0, 1.0),
        clamp_factor: bool | LINKABLE = True,
        factor_mode: str = "uniform",
    ) -> Mix:
        match factor_mode:
            case "uniform":
                builder = cls(
                    default_input_id="A_Vector",
                    default_output_id="Result_Vector",
                    Factor_Float=factor,
                    A_Vector=a,
                    B_Vector=b,
                    data_type=DataTypes.VECTOR,
                )
            case "non_uniform":
                builder = cls(
                    default_input_id="A_Vector",
                    default_output_id="Result_Vector",
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
    def color_(
        cls,
        factor: float | LINKABLE = 0.5,
        a: tuple[float, float, float] | list[float] | LINKABLE = (0.0, 0.0, 0.0),
        b: tuple[float, float, float] | list[float] | LINKABLE = (1.0, 1.0, 1.0),
        blend_type: str = "add",
        clamp_factor: bool = True,
        clamp_result: bool = True,
    ) -> Mix:
        builder = cls(
            default_input_id="A_Color",
            default_output_id="Result_Color",
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
    def rotation_(
        cls,
        a: tuple[float, float, float, float] | list[float] | LINKABLE | None = None,
        b: tuple[float, float, float, float] | list[float] | LINKABLE | None = None,
        factor: float | LINKABLE = 0.5,
        clamp_factor: bool = True,
    ) -> Mix:
        builder = cls(
            default_input_id="A_Rotation",
            default_output_id="Result_Rotation",
            Factor_Float=factor,
            A_Rotation=a,
            B_Rotation=b,
            data_type=DataTypes.ROTATION,
        )
        builder.node.clamp_factor = clamp_factor
        return builder
