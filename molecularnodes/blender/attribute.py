from dataclasses import dataclass
from enum import Enum
from typing import Type
from functools import reduce

import bpy
import numpy as np
from numpy.ma.core import prod


@dataclass
class AttributeTypeInfo:
    dname: str
    dtype: type
    width: int


@dataclass
class DomainInfo:
    name: str


class AttributeMismatchError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


# https://docs.blender.org/api/current/bpy_types_enum_items/attribute_domain_items.html#rna-enum-attribute-domain-items
class Domains:
    POINT = DomainInfo(name="POINT")
    EDGE = DomainInfo(name="EDGE")
    FACE = DomainInfo(name="FACE")
    CORNER = DomainInfo(name="CORNER")
    CURVE = DomainInfo(name="CURVE")
    INSTANCE = DomainInfo(name="INSTNANCE")
    LAYER = DomainInfo(name="LAYER")


@dataclass
class AttributeType:
    type_name: str
    value_name: str
    dtype: Type
    dimensions: tuple


# https://docs.blender.org/api/current/bpy_types_enum_items/attribute_type_items.html#rna-enum-attribute-type-items
class AttributeTypes(Enum):
    # https://docs.blender.org/api/current/bpy.types.FloatAttribute.html#bpy.types.FloatAttribute
    FLOAT = AttributeType(
        type_name="FLOAT", value_name="value", dtype=float, dimensions=(1,)
    )
    # https://docs.blender.org/api/current/bpy.types.FloatVectorAttribute.html#bpy.types.FloatVectorAttribute
    FLOAT_VECTOR = AttributeType(
        type_name="FLOAT_VECTOR", value_name="vector", dtype=float, dimensions=(3,)
    )
    # https://docs.blender.org/api/current/bpy.types.Float2Attribute.html#bpy.types.Float2Attribute
    FLOAT2 = AttributeType(
        type_name="FLOAT2", value_name="vector", dtype=float, dimensions=(2,)
    )
    # alternatively use color_srgb to get the color info in sRGB color space, otherwise linear color space
    # https://docs.blender.org/api/current/bpy.types.FloatColorAttributeValue.html#bpy.types.FloatColorAttributeValue
    FLOAT_COLOR = AttributeType(
        type_name="FLOAT_COLOR", value_name="color", dtype=float, dimensions=(4,)
    )
    # https://docs.blender.org/api/current/bpy.types.ByteColorAttribute.html#bpy.types.ByteColorAttribute
    # TODO unsure about this, int values are stored but float values are returned
    BYTE_COLOR = AttributeType(
        type_name="BYTE_COLOR", value_name="color", dtype=int, dimensions=(4,)
    )
    # https://docs.blender.org/api/current/bpy.types.QuaternionAttribute.html#bpy.types.QuaternionAttribute
    QUATERNION = AttributeType(
        type_name="QUATERNION", value_name="value", dtype=float, dimensions=(4,)
    )
    # https://docs.blender.org/api/current/bpy.types.IntAttribute.html#bpy.types.IntAttribute
    INT = AttributeType(type_name="INT", value_name="value", dtype=int, dimensions=(1,))
    # https://docs.blender.org/api/current/bpy.types.ByteIntAttributeValue.html#bpy.types.ByteIntAttributeValue
    INT8 = AttributeType(
        type_name="INT8", value_name="value", dtype=int, dimensions=(1,)
    )
    # https://docs.blender.org/api/current/bpy.types.Int2Attribute.html#bpy.types.Int2Attribute
    INT32_2D = AttributeType(
        type_name="INT32_2D", value_name="value", dtype=int, dimensions=(2,)
    )
    # https://docs.blender.org/api/current/bpy.types.Float4x4Attribute.html#bpy.types.Float4x4Attribute
    FLOAT4X4 = AttributeType(
        type_name="FLOAT4X4", value_name="value", dtype=float, dimensions=(4, 4)
    )
    # https://docs.blender.org/api/current/bpy.types.BoolAttribute.html#bpy.types.BoolAttribute
    BOOLEAN = AttributeType(
        type_name="BOOLEAN", value_name="value", dtype=bool, dimensions=(1,)
    )

def guess_atype_from_array(array: np.ndarray) -> AttributeType:
    if not isinstance(array, np.ndarray):
        raise ValueError(f"`array` must be a numpy array, not {type(array)=}")

    dtype = array.dtype
    shape = array.shape
    n_row = shape[0]

    # for 1D arrays we we use the float, int of boolean attribute types
    if shape == (n_row, 1) or shape == (n_row,):
        if np.issubdtype(dtype, np.int_):
            return AttributeTypes.INT.value
        elif np.issubdtype(dtype, np.float_):
            return AttributeTypes.FLOAT.value
        elif np.issubdtype(dtype, np.bool_):
            return AttributeTypes.BOOLEAN.value

    # for 2D arrays we use the float_vector, float_color, float4x4 attribute types
    elif shape == (n_row, 4, 4):
        return AttributeTypes.FLOAT4X4.value
    elif shape == (n_row, 3):
        return AttributeTypes.FLOAT_VECTOR.value
    elif shape == (n_row, 4):
        return AttributeTypes.FLOAT_COLOR.value

    # if we didn't match against anything return float
    return AttributeTypes.FLOAT.value

class Attribute:
    def __init__(self, attribute: bpy.types.Attribute):
        self.attribute = attribute
        self.n_attr = len(attribute.data)

    @property
    def atype(self):
        try:
            atype = AttributeTypes[self.attribute.data_type].value
        except KeyError:
            raise ValueError(f"Unknown attribute type: {self.attribute.data_type}")

        return atype

    @property
    def value_name(self):
        return self.atype.value_name

    @property
    def is_1d(self):
        return self.atype.dimensions == (1,)

    @property
    def type_name(self):
        return self.atype.type_name

    @property
    def shape(self):
        return (self.n_attr, *self.atype.dimensions)

    @property
    def dtype(self) -> Type:
        return self.atype.dtype

    @property
    def n_values(self) -> int:
        return np.prod(self.shape, dtype=int)


    def as_array(self) -> np.ndarray:
            """
            Returns the attribute data as a numpy array
            """
            # initialize empty 1D array that is needed to then be filled with values
            # from the Blender attribute
            array = np.zeros(self.n_values, dtype=self.dtype)
            self.attribute.data.foreach_get(self.value_name, array)

            # if the attribute has more than one dimension reshape the array before returning
            if self.is_1d:
                return array
            else:
                return array.reshape(self.shape)

    def __str__(self):
        return "Attribute: {}, type: {}, size: {}".format(self.attribute.name, self.type_name, self.shape)
