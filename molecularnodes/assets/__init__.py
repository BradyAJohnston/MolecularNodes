from pathlib import Path

from .template import install, uninstall
from . import data

ASSET_DIR = Path(__file__).resolve().parent
ADDON_DIR = ASSET_DIR.parent
MN_DATA_FILE = ASSET_DIR / "MN_data_file_4.2.blend"
