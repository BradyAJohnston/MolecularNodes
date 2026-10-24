import os
from typing import Iterable
import bpy
from nodebpy.builder.asset import AssetNodeGroup


def link_asset_groups(groups: Iterable[type]) -> None:
    """
    Link the node groups of several asset classes, reading each asset library once.

    Each asset group otherwise reads its whole library when it is first created, a
    fixed cost of tens of milliseconds however few groups are taken from it. Linking
    the groups a tree is about to use up front means paying that once. Classes that
    aren't asset groups, and groups already in the file, are skipped; the asset
    classes then reuse the linked groups.

    Mirrors `nodebpy.builder.link_assets`, which can replace this once released.
    """
    names_by_library: dict[str, set[str]] = {}
    for cls in groups:
        if not (isinstance(cls, type) and issubclass(cls, AssetNodeGroup)):
            continue
        existing = bpy.data.node_groups.get(cls._asset_name)
        if existing is not None and existing.bl_idname == cls._tree_idname:
            continue
        path = cls._library.path()
        names_by_library.setdefault(path, set()).add(cls._asset_name)

    for path, names in names_by_library.items():
        # without a built library the asset classes build from source instead
        if not os.path.exists(path):
            continue
        with bpy.data.libraries.load(path, link=True, pack=True, assets_only=True) as (
            src,
            dst,
        ):
            dst.node_groups = sorted(name for name in names if name in src.node_groups)
