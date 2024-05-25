import bpy
from . import parse

from .parse.ensemble import Ensemble
from typing import Union, Set
from pathlib import Path

bpy.types.Scene.MN_import_star_file_path = bpy.props.StringProperty(
    name="File",
    description="File path for the `.star` file to import.",
    subtype="FILE_PATH",
    maxlen=0,
)
bpy.types.Scene.MN_import_star_file_name = bpy.props.StringProperty(
    name="Name",
    description="Name of the created object.",
    default="NewStarInstances",
    maxlen=0,
)


def load(
    file_path: Union[str, Path], name: str = "NewStarInstances", node_setup: bool = True, world_scale: float = 0.01
) -> Ensemble:
    ensemble = parse.StarFile.from_starfile(file_path)
    ensemble.create_model(name=name, node_setup=node_setup, world_scale=world_scale)

    return ensemble


class MN_OT_Import_Star_File(bpy.types.Operator):  # type: ignore
    bl_idname = "mn.import_star_file"
    bl_label = "Load"
    bl_description = "Will import the given file, setting up the points to instance an object."
    bl_options = {"REGISTER"}

    @classmethod
    def poll(cls, context: bpy.types.Context) -> bool:
        return True

    def execute(self, context: bpy.types.Context) -> Set[str]:
        scene = context.scene
        load(
            file_path=scene.MN_import_star_file_path,
            name=scene.MN_import_star_file_name,
            node_setup=True,
        )
        return {"FINISHED"}


def panel(layout: bpy.types.UILayout, scene: bpy.types.Scene) -> None:
    layout.label(text="Load Star File", icon="FILE_TICK")
    layout.separator()
    row_import = layout.row()
    row_import.prop(scene, "MN_import_star_file_name")
    layout.prop(scene, "MN_import_star_file_path")
    row_import.operator("mn.import_star_file")
