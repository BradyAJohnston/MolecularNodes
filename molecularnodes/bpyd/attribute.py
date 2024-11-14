from dataclasses import dataclass
from enum import Enum
from typing import Type
import bpy
import numpy as np

from pathlib import Path


def path_resolve(path: str | Path) -> Path:
    if isinstance(path, str):
        return Path(bpy.path.abspath(path))
    elif isinstance(path, Path):
        return Path(bpy.path.abspath(str(path)))
    else:
        raise ValueError(f"Unable to resolve path: {path}")


@dataclass
class AttributeTypeInfo:
    dname: str
    dtype: type
    width: int


@dataclass
class DomainType:
    name: str

    def __str__(self):
        return self.name


class AttributeMismatchError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


# https://docs.blender.org/api/current/bpy_types_enum_items/attribute_domain_items.html#rna-enum-attribute-domain-items
class Domains:
    POINT = DomainType(name="POINT")
    EDGE = DomainType(name="EDGE")
    FACE = DomainType(name="FACE")
    CORNER = DomainType(name="CORNER")
    CURVE = DomainType(name="CURVE")
    INSTANCE = DomainType(name="INSTNANCE")
    LAYER = DomainType(name="LAYER")


@dataclass
class AttributeType:
    type_name: str
    value_name: str
    dtype: Type
    dimensions: tuple

    def __str__(self) -> str:
        return self.type_name


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
            return AttributeTypes.INT
        elif np.issubdtype(dtype, np.float_):
            return AttributeTypes.FLOAT
        elif np.issubdtype(dtype, np.bool_):
            return AttributeTypes.BOOLEAN

    # for 2D arrays we use the float_vector, float_color, float4x4 attribute types
    elif shape == (n_row, 4, 4):
        return AttributeTypes.FLOAT4X4
    elif shape == (n_row, 3):
        return AttributeTypes.FLOAT_VECTOR
    elif shape == (n_row, 4):
        return AttributeTypes.FLOAT_COLOR

    # if we didn't match against anything return float
    return AttributeTypes.FLOAT


class Attribute:
    """
    Wrapper around a Blender attribute to provide a more convenient interface with numpy arrays
    """

    def __init__(self, attribute: bpy.types.Attribute):
        self.attribute = attribute
        self.n_attr = len(attribute.data)
        self.atype = AttributeTypes[self.attribute.data_type].value

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

    @classmethod
    def from_object(
        cls,
        obj: bpy.types.Object,
        name: str,
        atype: AttributeType,
        domain: DomainType,
    ):
        att = obj.data.get(name)
        if att is None:
            att = obj.data.attributes.new(
                name=name, type=atype.value.type_name, domain=domain.value.name
            )
        return Attribute(att)

    def from_array(self, array: np.ndarray) -> None:
        """
        Set the attribute data from a numpy array
        """
        if array.shape != self.shape:
            raise ValueError(
                f"Array shape {array.shape} does not match attribute shape {self.shape}"
            )

        self.attribute.data.foreach_set(self.value_name, array.reshape(-1))

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
        return "Attribute: {}, type: {}, size: {}".format(
            self.attribute.name, self.type_name, self.shape
        )


def store_named_attribute(
    obj: bpy.types.Object,
    data: np.ndarray,
    name: str,
    atype: str | AttributeType | None = None,
    domain: str | DomainType = Domains.POINT,
    overwrite: bool = True,
) -> bpy.types.Attribute:
    """
    Adds and sets the values of an attribute on the object.

    Parameters
    ----------
    obj : bpy.types.Object
        The Blender object.
    name : str
        The name of the attribute.
    data : np.ndarray
        The attribute data as a numpy array.
    atype : str, AttributeType, optional
        The attribute type to store the data as. One of the AttributeType enums or a string
        of the same name.
        'FLOAT_VECTOR', 'FLOAT_COLOR', 'FLOAT4X4', 'QUATERNION', 'FLOAT', 'INT', 'BOOLEAN'
    domain : str, optional
        The domain of the attribute. Defaults to 'POINT'. Currently, only 'POINT', 'EDGE',
        and 'FACE' have been tested.
    overwrite : bool
        Setting to false will create a new attribute if the given name is already an
        attribute on the mesh.

    Returns
    -------
    bpy.types.Attribute
        The added attribute.
    """

    if isinstance(atype, str):
        try:
            atype = AttributeTypes[atype]
        except KeyError:
            raise ValueError(
                f"Given data type {atype=} does not match any of the possible attribute types: {list(AttributeTypes)=}"
            )

    if atype is None:
        atype = guess_atype_from_array(data)

    attribute = obj.data.attributes.get(name)  # type: ignore
    if not attribute or not overwrite:
        attribute = obj.data.attributes.new(name, atype.value.type_name, str(domain))

    if len(data) != len(attribute.data):
        raise AttributeMismatchError(
            f"Data length {len(data)}, dimensions {data.shape} does not equal the size of the target domain {domain}, len={len(attribute.data)=}"
        )

    # the 'foreach_set' requires a 1D array, regardless of the shape of the attribute
    # so we have to flatten it first
    attribute.data.foreach_set(atype.value.value_name, data.reshape(-1))

    # attribute bug is fixed in 4.3+
    if bpy.app.version_string.startswith("4.2"):
        # TODO remove in later updates
        # The updating of data doesn't work 100% of the time (see:
        # https://projects.blender.org/blender/blender/issues/118507) so this resetting of a
        # single vertex is the current fix. Not great as I can see it breaking when we are
        # missing a vertex - but for now we shouldn't be dealing with any situations where this
        # is the case For now we will set a single vert to it's own position, which triggers a
        # proper refresh of the object data.
        try:
            obj.data.vertices[0].co = obj.data.vertices[0].co  # type: ignore
        except AttributeError:
            obj.data.update()  # type: ignore

    return attribute


def evaluate_object(obj: bpy.types.Object):
    "Return an object which has the modifiers evaluated."
    obj.update_tag()
    return obj.evaluated_get(bpy.context.evaluated_depsgraph_get())


def named_attribute(
    obj: bpy.types.Object, name="position", evaluate=False
) -> np.ndarray:
    """
    Get the named attribute data from the object, optionally evaluating modifiers first.

    Parameters:
        object (bpy.types.Object): The Blender object.
        name (str, optional): The name of the attribute. Defaults to 'position'.

    Returns:
        np.ndarray: The attribute data as a numpy array.
    """
    if evaluate:
        obj = evaluate_object(obj)
    verbose = False
    try:
        attr = Attribute(obj.data.attributes[name])
    except KeyError:
        message = f"The selected attribute '{name}' does not exist on the mesh."
        if verbose:
            message += f"Possible attributes are: {obj.data.attributes.keys()}"

        raise AttributeError(message)

    return attr.as_array()


def remove_named_attribute(
    obj: bpy.types.Object, name: str, domain: str | DomainType = Domains.POINT
):
    try:
        attr = obj.data.attributes[name]
        obj.data.attributes.remove(attr)
    except KeyError:
        raise AttributeError(
            f"The selected attribute '{name}' does not exist on the mesh."
        )
