from pathlib import Path
from . import data
from .template import install, uninstall

ASSET_DIR = Path(__file__).resolve().parent
ADDON_DIR = ASSET_DIR.parent
MN_DATA_FILE = ASSET_DIR / "MN_data_file_4.4.blend"


__all__ = ["data", "install", "uninstall", "MN_DATA_FILE"]
