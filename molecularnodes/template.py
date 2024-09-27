import os
import shutil

import bpy

from .utils import ADDON_DIR

SUBFOLDER = "Molecular Nodes"


def is_installed():
    base_path = bpy.utils.user_resource(
        "SCRIPTS", path=os.path.join("startup", "bl_app_templates_user"), create=False
    )
    molecular_nodes_path = os.path.join(base_path, "Molecular Nodes")
    return os.path.exists(molecular_nodes_path)


def install():
    # Construct the base path for app templates
    base_path = bpy.utils.user_resource(
        "SCRIPTS", path=os.path.join("startup", "bl_app_templates_user"), create=True
    )

    # Ensure the base path exists
    if not os.path.exists(base_path):
        os.makedirs(base_path)

    # Construct the full path to the "Molecular Nodes" subfolder
    path_app_templates = os.path.join(base_path, SUBFOLDER)

    # Create the "Molecular Nodes" folder if it doesn't exist
    if not os.path.exists(path_app_templates):
        os.makedirs(path_app_templates)

    # Define the path to the startup file
    startup_file = os.path.join(
        os.path.abspath(ADDON_DIR), "assets", "template", "startup.blend"
    )

    # Copy the startup file to the "Molecular Nodes" folder
    shutil.copy(startup_file, path_app_templates)

    # Refresh Blender's script paths to recognize the new template
    bpy.utils.refresh_script_paths()


def uninstall():
    for folder in bpy.utils.app_template_paths():
        path = os.path.join(os.path.abspath(folder), SUBFOLDER)
        if "Molecular Nodes" not in path:
            continue

        if os.path.exists(path):
            shutil.rmtree(path)
    bpy.utils.refresh_script_paths()
