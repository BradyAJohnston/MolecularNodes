import shutil
from pathlib import Path
import bpy

SUBFOLDER = "Molecular Nodes"


def is_installed():
    base_path = Path(
        bpy.utils.user_resource(
            "SCRIPTS", path=str(Path("startup") / "bl_app_templates_user"), create=False
        )
    )
    molecular_nodes_path = base_path / "Molecular Nodes"
    return molecular_nodes_path.exists()


def install():
    base_path = Path(
        bpy.utils.user_resource(
            "SCRIPTS", path=str(Path("startup") / "bl_app_templates_user"), create=True
        )
    )

    base_path.mkdir(parents=True, exist_ok=True)
    path_app_templates = base_path / SUBFOLDER
    path_app_templates.mkdir(parents=True, exist_ok=True)
    startup_file = Path(__file__).parent / "template/startup.blend"
    shutil.copy(startup_file, path_app_templates)
    bpy.utils.refresh_script_paths()


def uninstall():
    for folder in bpy.utils.app_template_paths():
        path = Path(folder).absolute() / SUBFOLDER
        if "Molecular Nodes" not in str(path):
            continue

        if path.exists():
            shutil.rmtree(path)
    bpy.utils.refresh_script_paths()
