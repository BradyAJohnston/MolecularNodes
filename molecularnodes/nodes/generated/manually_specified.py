"""
Some of the nodes need to be manually specified because they are a bit tricky to generate automatically.
"""

from __future__ import annotations
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
        return self.node.outputs["Value"]

    @classmethod
    def float_(
        cls,
        min: float | LINKABLE = 0.0,
        max: float | LINKABLE = 1.0,
        id: int | LINKABLE | None = None,
        seed: int | LINKABLE = 1,
    ) -> NodeBuilder:
        node = cls(
            Min_001=min, Max_001=max, id=id, seed=seed, data_type=DataTypes.FLOAT
        )
        node._default_output_id = "Value_001"
        return node

    @classmethod
    def integer_(
        cls,
        min: int | LINKABLE = 0,
        max: int | LINKABLE = 1,
        id: int | LINKABLE | None = None,
        seed: int | LINKABLE = 1,
    ) -> NodeBuilder:
        node = cls(Min_002=min, Max_002=max, id=id, seed=seed, data_type=DataTypes.INT)
        node._default_output_id = "Value_002"
        return node

    @classmethod
    def boolean_(
        cls,
        probability: float | LINKABLE = 0.5,
        id: int | LINKABLE | None = None,
        seed: int | LINKABLE = 1,
    ) -> NodeBuilder:
        node = cls(Probability=probability, id=id, seed=seed, data_type=DataTypes.BOOL)
        node._default_output_id = "Value_003"
        return node

    @classmethod
    def vector_(
        cls,
        min: tuple[float, float, float] | list[float] | LINKABLE = (0.0, 0.0, 0.0),
        max: tuple[float, float, float] | list[float] | LINKABLE = (1.0, 1.0, 1.0),
        id: int | LINKABLE | None = None,
        seed: int | LINKABLE = 1,
    ) -> NodeBuilder:
        node = cls(Min=min, Max=max, id=id, seed=seed, data_type=DataTypes.VECTOR)
        node._default_output_id = "Value"
        return node


class Mix(NodeBuilder):
    """Mix values by a factor"""

    name = "ShaderNodeMix"

    def __init__(self, data_type: str | None = None, **kwargs):
        super().__init__()
        self.node.data_type = data_type
        key_args = {}
        key_args.update(kwargs)
        self._establish_links(**key_args)

    @property
    def result(self) -> NodeSocket:
        """Output socket: Result"""
        return self.node.outputs["Result"]

    @classmethod
    def float_(
        cls,
        factor: float | LINKABLE = 0.5,
        a: float | LINKABLE = 0.0,
        b: float | LINKABLE = 0.0,
        clamp_factor: bool | LINKABLE = True,
    ) -> NodeBuilder:
        node = cls(
            Factor_Float=factor, A_Float=a, B_Float=b, id=id, data_type=DataTypes.FLOAT
        )
        node._default_input_id = "A_Float"
        node._default_output_id = "Result_Float"
        node.node.clamp_factor = clamp_factor
        return node

    @classmethod
    def vector_(
        cls,
        factor: float | LINKABLE = 0.5,
        a: tuple[float, float, float] | list[float] | LINKABLE = (0.0, 0.0, 0.0),
        b: tuple[float, float, float] | list[float] | LINKABLE = (1.0, 1.0, 1.0),
        clamp_factor: bool | LINKABLE = True,
        factor_mode: str = "uniform",
    ) -> NodeBuilder:
        match factor_mode:
            case "uniform":
                node = cls(
                    Factor_Float=factor,
                    A_Vector=a,
                    B_Vector=b,
                    id=id,
                    data_type=DataTypes.VECTOR,
                )
            case "non_uniform":
                node = cls(
                    Factor_Vector=factor,
                    A_Vector=a,
                    B_Vector=b,
                    id=id,
                    data_type=DataTypes.VECTOR,
                )
            case _:
                raise ValueError(f"Invalid factor mode: {factor_mode}")

        node._default_input_id = "A_Vector"
        node._default_output_id = "Result_Vector"
        node.node.clamp_factor = clamp_factor
        return node

    @classmethod
    def color_(
        cls,
        factor: float | LINKABLE = 0.5,
        a: tuple[float, float, float] | list[float] | LINKABLE = (0.0, 0.0, 0.0),
        b: tuple[float, float, float] | list[float] | LINKABLE = (1.0, 1.0, 1.0),
        blend_type: str = "add",
        clamp_factor: bool = True,
        clamp_result: bool = True,
    ) -> NodeBuilder:
        node = cls(
            Factor_Float=factor,
            A_Color=a,
            B_Color=b,
            id=id,
            data_type=DataTypes.COLOR,
        )
        node._default_input_id = "A_Color"
        node._default_output_id = "Result_Color"
        node.node.blend_type = blend_type.capitalize()
        node.node.clamp_factor = clamp_factor
        node.node.clamp_result = clamp_result
        return node

    @classmethod
    def rotation_(
        cls,
        a: tuple[float, float, float, float] | list[float] | LINKABLE | None = None,
        b: tuple[float, float, float, float] | list[float] | LINKABLE | None = None,
        factor: float | LINKABLE = 0.5,
        clamp_factor: bool = True,
    ) -> NodeBuilder:
        node = cls(
            Factor_Float=factor,
            A_Rotation=a,
            B_Rotation=b,
            data_type=DataTypes.ROTATION,
        )
        node._default_input_id = "A_Rotation"
        node._default_output_id = "Result_Rotation"
        node.node.clamp_factor = clamp_factor
        return node
