from pathlib import Path
from . import data
from .template import install, uninstall

ASSET_DIR = Path(__file__).resolve().parent
ADDON_DIR = ASSET_DIR.parent
MN_DATA_FILE = ASSET_DIR / "node_data_file.blend"


__all__ = ["data", "install", "uninstall", "MN_DATA_FILE"]
