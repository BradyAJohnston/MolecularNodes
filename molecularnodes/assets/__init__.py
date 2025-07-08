from pathlib import Path
from . import data
from .template import install, uninstall

ASSET_DIR = Path(__file__).resolve().parent
ADDON_DIR = ASSET_DIR.parent
ASSET_FILE = ASSET_DIR / "asset_file.blend"


__all__ = ["data", "install", "uninstall", "ASSET_FILE"]
