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
from .utils import get_all_class_annotations


def create_annotation_type_inputs(
    annotation_class: BaseAnnotation,
) -> bpy.types.PropertyGroup:
    """Create a dynamic PropertyGroup with annotation type inputs"""
    attributes = {"__annotations__": {}}
    py_annotations = get_all_class_annotations(annotation_class)
    for key, atype in py_annotations.items():
        if key == "name":
            continue
        if atype.__name__ == "str":
            prop = StringProperty(default=getattr(annotation_class, key, ""))
        elif atype.__name__ == "bool":
            prop = BoolProperty(default=getattr(annotation_class, key, False))
        elif atype.__name__ == "int":
            prop = IntProperty(default=getattr(annotation_class, key, 0))
        elif atype.__name__ == "float":
            prop = FloatProperty(default=getattr(annotation_class, key, 0.0))
        else:
            continue
        attributes["__annotations__"][key] = prop
    AnnotationInputs = type("AnnotationInputs", (bpy.types.PropertyGroup,), attributes)
    return AnnotationInputs


def create_property_interface(prop: bpy.types.bpy_struct, path: str) -> property:
    """Create a property() interface for a blender property"""

    def getter(self):
        return getattr(prop, path)

    def setter(self, value):
        setattr(prop, path, value)
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

    offset_x: FloatProperty(
        name="Offset X",
        description="Offset in X direction",
        default=0,
        min=0.0,
        max=1.0,
    )  # type: ignore

    offset_y: FloatProperty(
        name="offset Y",
        description="Offset in Y direction",
        default=0,
        min=0.0,
        max=1.0,
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
