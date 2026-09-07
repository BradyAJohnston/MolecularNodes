"""
Run custom setup when MolecularNodes node group assets are added to a tree.

Since the node groups are shipped as regular assets in an asset library, no
MolecularNodes code runs when a user drags one into a node tree. The
``blend_import_post`` handler fires for every link/append operation — including
asset drag-and-drop and Add-menu asset adds, whether the asset library's import
method appends, links or packs, and for re-adds that reuse an already imported
group — which gives us a hook to restore the conveniences the old custom add
operator provided (assigning a default material) and to add new ones (building
assembly data for the entity the node was dropped onto).

Setup runs in two phases, because the handler fires while the add operator is
still importing the datablock, before the node instance exists in the tree:

- ``tree_setup`` runs immediately on the imported ``NodeTree`` datablock, and
  only when the tree was appended — linked (and packed) trees are read-only.
- ``node_setup`` runs on each new node instance via a one-shot timer scheduled
  by the handler, which fires right after the operator has created the node.
  Node instances are local even when their node group is linked, so their
  input socket values and custom properties are always writable.

Callbacks must be idempotent — they may run more than once for the same
datablock and are also invoked on nodes created programmatically, so they
should no-op when there is nothing to do (a value already set, no matching
entity, ...).
"""

import logging
import re
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Callable
import bpy
from bpy.app.handlers import persistent
from ..assets import MN_DATA_FILE

logger = logging.getLogger(__name__)

# stamped onto a node instance once its node_setup has been attempted, so a
# later add of the same asset doesn't re-run setup on pre-existing nodes
_MARKER = "mn_on_add_done"
_DUPLICATE_SUFFIX = re.compile(r"\.\d{3}$")


@dataclass(frozen=True)
class AssetOnAdd:
    """Setup callbacks to run when a node group asset is added to a tree."""

    tree_setup: Callable[[bpy.types.GeometryNodeTree], None] | None = None
    node_setup: (
        Callable[[bpy.types.GeometryNodeGroup, bpy.types.Object | None], None] | None
    ) = None


_registry: dict[str, AssetOnAdd] = {}
_predicates: list[tuple[Callable[[str], bool], AssetOnAdd]] = []
# (imported group name, drop-target tree name, context object name) awaiting
# node_setup once the add operator has created the node instance
_pending: list[tuple[str, str | None, str | None]] = []

# guards against re-entering the handler when a setup callback itself appends
# data (e.g. materials) from another .blend file; attached to
# bpy.types.WindowManager in ui.addon so the flag lives as a documented,
# runtime-only Blender property rather than in module state
importing_assets_property = bpy.props.BoolProperty(
    name="MN Importing Assets",
    description=(
        "True while MolecularNodes' own on-add setup is importing data from a "
        ".blend file (e.g. appending materials), so that its blend_import_post "
        "handler ignores imports it caused itself. Runtime-only, never saved"
    ),
    default=False,
    options={"HIDDEN", "SKIP_SAVE"},
)


def _in_setup() -> bool:
    return bpy.context.window_manager.mn_importing_assets


@contextmanager
def _own_imports_flagged():
    """Flag imports we cause ourselves so the handler ignores them."""
    wm = bpy.context.window_manager
    wm.mn_importing_assets = True
    try:
        yield
    finally:
        wm.mn_importing_assets = False


def register_on_add(
    match: str | Callable[[str], bool],
    *,
    tree_setup: Callable[[bpy.types.GeometryNodeTree], None] | None = None,
    node_setup: Callable[[bpy.types.GeometryNodeGroup, bpy.types.Object | None], None]
    | None = None,
) -> None:
    """
    Register setup callbacks for a node group asset.

    Parameters
    ----------
    match : str | Callable[[str], bool]
        Either the exact name of the asset's node group, or a predicate called
        with the (duplicate-suffix stripped) name of each imported group.
    tree_setup : Callable, optional
        Called with the imported ``NodeTree`` datablock during the import,
        before the dropped node instance exists. Only runs for appended trees,
        as linked ones are read-only.
    node_setup : Callable, optional
        Called with each newly added node instance and the object whose modifier
        tree it was added to (or ``None`` when that can't be determined).
    """
    entry = AssetOnAdd(tree_setup=tree_setup, node_setup=node_setup)
    if callable(match):
        _predicates.append((match, entry))
    else:
        _registry[match] = entry


def _lookup(name: str) -> AssetOnAdd | None:
    name = _DUPLICATE_SUFFIX.sub("", name)
    entry = _registry.get(name)
    if entry is not None:
        return entry
    for predicate, entry in _predicates:
        if predicate(name):
            return entry
    return None


def _is_mn_source(item) -> bool:
    """Whether an import item came from the MolecularNodes asset data file."""
    library = item.source_library
    if library is None:
        # source not recorded (e.g. append without keeping the library);
        # fall back to matching on name alone
        return True
    try:
        source = Path(bpy.path.abspath(library.filepath)).resolve()
    except (OSError, ValueError):
        return False
    return source == Path(MN_DATA_FILE).resolve()


@persistent
def node_asset_import_post(ctx: bpy.types.BlendImportContext) -> None:
    """
    Dispatch registered setup for MN node group assets after they are imported.

    Registered on ``bpy.app.handlers.blend_import_post``, which fires for asset
    drag-and-drop as well as any other link/append/pack, including adds that
    reuse an already imported group.
    """
    if _in_setup():
        return
    schedule = False
    for item in ctx.import_items:
        # only directly requested node trees, not their dependencies
        if item.id_type != "NODETREE" or item.import_info:
            continue
        tree = item.id
        if tree is None:
            continue
        if not _is_mn_source(item):
            continue
        entry = _lookup(item.name)
        if entry is None:
            continue
        # linked (rather than appended) trees are read-only, so only run
        # datablock-level setup on appended copies
        if entry.tree_setup is not None and tree.library is None:
            _run_tree_setup(entry.tree_setup, tree)
        if entry.node_setup is not None:
            # the node instance doesn't exist yet — capture where the drop is
            # happening while the operator context is available, and defer
            space = getattr(bpy.context, "space_data", None)
            target = getattr(space, "edit_tree", None)
            obj = getattr(bpy.context, "object", None)
            _pending.append(
                (
                    _DUPLICATE_SUFFIX.sub("", tree.name),
                    target.name if target is not None else None,
                    obj.name if obj is not None else None,
                )
            )
            schedule = True
    if schedule and not bpy.app.timers.is_registered(_process_pending):
        bpy.app.timers.register(_process_pending, first_interval=0.0)


def _run_tree_setup(
    setup: Callable[[bpy.types.GeometryNodeTree], None],
    tree: bpy.types.GeometryNodeTree,
) -> None:
    try:
        with _own_imports_flagged():
            setup(tree)
    except Exception:
        logger.exception(f"tree_setup failed for imported node group {tree.name!r}")


def _object_using_tree(tree: bpy.types.NodeTree) -> bpy.types.Object | None:
    """Find the object using this tree as a geometry nodes modifier, if any."""
    for obj in bpy.data.objects:
        for mod in obj.modifiers:
            if mod.type == "NODES" and mod.node_group == tree:
                return obj
    return None


def _process_pending() -> None:
    """Run deferred node_setup callbacks on newly added node instances."""
    entries, _pending[:] = list(_pending), []
    for group_name, target_name, object_name in entries:
        entry = _lookup(group_name)
        if entry is None or entry.node_setup is None:
            continue
        if target_name is not None and target_name in bpy.data.node_groups:
            hosts = [bpy.data.node_groups[target_name]]
        else:
            hosts = list(bpy.data.node_groups)
        fallback_obj = (
            bpy.data.objects.get(object_name) if object_name is not None else None
        )
        for host in hosts:
            if not isinstance(host, bpy.types.GeometryNodeTree):
                continue
            for node in host.nodes:
                if (
                    node.bl_idname != "GeometryNodeGroup"
                    or node.node_tree is None
                    # match by name so appended, linked and packed copies of
                    # the same asset group are all handled
                    or _DUPLICATE_SUFFIX.sub("", node.node_tree.name) != group_name
                    or node.get(_MARKER)
                ):
                    continue
                # a node in a nested group has no modifier of its own, so fall
                # back to the object that was active when the asset was dropped
                obj = _object_using_tree(host) or fallback_obj
                try:
                    with _own_imports_flagged():
                        entry.node_setup(node, obj)
                except Exception:
                    logger.exception(
                        f"node_setup failed for node {node.name!r} in {host.name!r}"
                    )
                node[_MARKER] = True
    return None


def unregister_pending() -> None:
    """Drop queued work and the timer; called when the add-on is unregistered."""
    _pending.clear()
    if bpy.app.timers.is_registered(_process_pending):
        bpy.app.timers.unregister(_process_pending)


def _style_node_material(
    node: bpy.types.GeometryNodeGroup, obj: bpy.types.Object | None
) -> None:
    """
    Give a newly added style node's Material input the default MN material.

    Set on the node instance rather than the tree interface so it works for
    linked and packed node groups too, matching what the old custom add
    operator did.
    """
    from .material import add_all_materials

    socket = next((s for s in node.inputs if s.bl_idname == "NodeSocketMaterial"), None)
    if socket is None or socket.default_value is not None:
        return
    socket.default_value = add_all_materials()["MN Default"]


register_on_add(lambda name: name.startswith("Style "), node_setup=_style_node_material)
