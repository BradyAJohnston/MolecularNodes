from bpy.types import Operator
from bpy.props import StringProperty, BoolProperty, EnumProperty

from ...style import STYLE_ITEMS


class TrajectoryImportOperator(Operator):
    bl_label = "Import"
    bl_description = "Will import the given file and toplogy."
    bl_options = {"REGISTER"}

    topology: StringProperty(  # type: ignore
        name="Toplogy",
        description="File path for the topology file",
        subtype="FILE_PATH",
        maxlen=0,
    )
    trajectory: StringProperty(  # type: ignore
        name="Trajectory",
        description="File path for the trajectory file",
        subtype="FILE_PATH",
        maxlen=0,
    )
    name: StringProperty(  # type: ignore
        name="Name",
        description="Name for the object that will be created and linked to the trajectory",
        default="NewOrigami",
        maxlen=0,
    )
    style: EnumProperty(  # type: ignore
        name="Style",
        description="Starting style on import",
        default="spheres",
        items=STYLE_ITEMS,
    )
