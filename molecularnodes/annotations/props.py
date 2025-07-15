import typing
import bpy
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    FloatVectorProperty,
    IntProperty,
    StringProperty,
)
from ..blender.utils import viewport_tag_redraw
from .base import BaseAnnotation
from .utils import get_all_class_annotations, get_blender_supported_type


def _create_update_callback(func, prop_name):
    """Create an update property callback function"""
    if func is None:
        return None
    return lambda s, c: func(s, c, prop_name)


def create_annotation_type_inputs(
    annotation_class: BaseAnnotation,
    update_callback=None,
) -> bpy.types.PropertyGroup:
    """Create a dynamic PropertyGroup with annotation type inputs"""
    attributes = {"__annotations__": {}}
    py_annotations = get_all_class_annotations(annotation_class)
    for attr, atype in py_annotations.items():
        if attr == "name":
            continue
        stype = get_blender_supported_type(atype)
        if stype is str:
            prop = StringProperty(
                default=getattr(annotation_class, attr, ""),
                update=_create_update_callback(update_callback, attr),
            )
        elif stype is bool:
            prop = BoolProperty(
                default=getattr(annotation_class, attr, False),
                update=_create_update_callback(update_callback, attr),
            )
        elif stype is int:
            prop = IntProperty(
                default=getattr(annotation_class, attr, 0),
                update=_create_update_callback(update_callback, attr),
            )
        elif stype is float:
            prop = FloatProperty(
                default=getattr(annotation_class, attr, 0.0),
                update=_create_update_callback(update_callback, attr),
            )
        else:
            continue
        attributes["__annotations__"][attr] = prop
    # add the uuid to the annotation for lookup during update callback
    attributes["__annotations__"]["uuid"] = StringProperty()
    # add a boolean to indicate if validations succeeded
    attributes["__annotations__"]["valid_inputs"] = BoolProperty(default=True)
    AnnotationInputs = type("AnnotationInputs", (bpy.types.PropertyGroup,), attributes)
    return AnnotationInputs


def create_property_interface(
    prop: bpy.types.bpy_struct,
    attr: str,
    atype: typing.Any = None,
    instance: BaseAnnotation = None,
) -> property:
    """Create a property() interface for a blender property"""

    nbattr = f"_{attr}"  # non blender property

    def getter(self):
        # return the non blender property if set
        if hasattr(instance, nbattr):
            return getattr(instance, nbattr)
        return getattr(prop, attr)

    def setter(self, value):
        if atype is not None and not isinstance(value, atype):
            # not a supported blender type property
            # set a non blender property and validate
            setattr(instance, nbattr, value)
            instance.validate()
        else:
            # blender property
            setattr(prop, attr, value)
            # clear any corresponding non blender property
            if hasattr(instance, nbattr):
                delattr(instance, nbattr)
        viewport_tag_redraw()

    return property(getter, setter, doc=getattr(prop, "description", None))


class BaseAnnotationProperties(bpy.types.PropertyGroup):
    """Base annotation properties"""

    # annotation name (label)
    label: StringProperty()  # type: ignore

    # annotation type
    type: StringProperty()  # type: ignore

    # common properties
    visible: BoolProperty(
        name="Visible",
        description="Visibility of the annotation",
        default=True,
    )  # type: ignore

    text_color: FloatVectorProperty(
        name="Text Color",
        description="Font color",
        subtype="COLOR",
        size=4,
        default=(1, 1, 1, 1),
        min=0.0,
        max=1.0,
    )  # type: ignore

    text_size: IntProperty(
        name="Text Size",
        description="Font size",
        default=16,
        min=0,
        max=128,
    )  # type: ignore

    text_align: EnumProperty(
        name="Text Alignment",
        description="Alignment of text",
        items=(
            ("center", "Center", "Center"),
            ("left", "Left", "Left"),
            ("right", "Right", "Right"),
        ),
    )  # type: ignore

    text_rotation: FloatProperty(
        name="Text Rotation",
        description="Text rotation in degrees",
        default=0.0,
        min=0.0,
        max=360.0,
    )  # type: ignore

    offset_x: IntProperty(
        name="Offset X",
        description="Offset in X direction",
        default=0,
    )  # type: ignore

    offset_y: IntProperty(
        name="offset Y",
        description="Offset in Y direction",
        default=0,
    )  # type: ignore

    line_color: FloatVectorProperty(
        name="Line Color",
        description="Line color",
        subtype="COLOR",
        size=4,
        default=(1, 1, 1, 1),
        min=0.0,
        max=1.0,
    )  # type: ignore

    line_width: FloatProperty(
        name="Line Width",
        description="Line width",
        default=1.0,
        min=0.0,
    )  # type: ignore

    arrow_size: IntProperty(
        name="Arrow Size",
        description="Arrow size",
        default=16,
        min=0,
    )  # type: ignore

    pointer_length: IntProperty(
        name="Pointer Length",
        description="Pointer length",
        default=0,
        min=0,
    )  # type: ignore
