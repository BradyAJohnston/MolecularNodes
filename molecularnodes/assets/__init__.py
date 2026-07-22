from pathlib import Path
import bpy
from . import data
from .template import install, uninstall

ASSET_DIR = Path(__file__).resolve().parent
ADDON_DIR = ASSET_DIR.parent
MN_DATA_FILE = ASSET_DIR / "node_data_file.blend"
LIBRARY_NAME = "Molecular Nodes"


__all__ = ["data", "install", "uninstall", "MN_DATA_FILE"]


def _libs(context: bpy.types.Context | None = None) -> bpy.types.AssetLibraryCollection:
    if not context:
        context = bpy.context
    return context.preferences.filepaths.asset_libraries


def _ensure_mn_asset_library() -> None:
    libs = _libs()
    if LIBRARY_NAME not in libs:
        libs.new(name=LIBRARY_NAME, directory=str(ASSET_DIR))


def _remove_mn_asset_library() -> None:
    libs = _libs()
    if LIBRARY_NAME in libs:
        libs.remove(libs[LIBRARY_NAME])
