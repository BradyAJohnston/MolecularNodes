"""
Run custom setup when MolecularNodes node group assets are added to a tree.

Since the node groups are shipped as regular assets in an asset library, no
MolecularNodes code runs when a user drags one into a node tree. The
``blend_import_post`` handler fires for every link/append operation — including
asset drag-and-drop, and re-drops that reuse an already-appended group — which
gives us a hook to restore the conveniences the old custom add operator
provided (assigning a default material) and to add new ones (building assembly
data for the entity the node was dropped onto).

Setup runs in two phases, because the handler fires while the drop operator is
still appending the datablock, before the node instance exists in the tree:

- ``tree_setup`` runs immediately on the imported ``NodeTree`` datablock. Any
  interface defaults set here are inherited by the node instance the drop
  operator creates a moment later.
- ``node_setup`` runs on each new node instance via a one-shot timer scheduled
  by the handler, which fires right after the operator has created the node.

Callbacks must be idempotent — they may run more than once for the same
datablock and are also invoked on nodes created programmatically, so they
should no-op when there is nothing to do (a value already set, no matching
entity, ...).
"""

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable
import bpy
from bpy.app.handlers import persistent
from ..assets import MN_DATA_FILE

logger = logging.getLogger(__name__)

# stamped onto a node instance once its node_setup has been attempted, so a
# later drop of the same asset doesn't re-run setup on pre-existing nodes
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
# (imported tree name, drop-target tree name, context object name) awaiting
# node_setup once the drop operator has created the node instance
_pending: list[tuple[str, str | None, str | None]] = []
# guards against re-entering the handler when a callback itself appends data
# (e.g. materials) from another .blend file
_in_setup = False


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
        Called with the imported ``NodeTree`` datablock during the import, before
        the dropped node instance exists.
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
    drag-and-drop as well as any other link/append, including drops that reuse
    an already-appended group.
    """
    if _in_setup:
        return
    schedule = False
    for item in ctx.import_items:
        # only directly requested node trees, not their dependencies
        if item.id_type != "NODETREE" or item.import_info:
            continue
        tree = item.id
        # linked (rather than appended) data is read-only, so leave it alone
        if tree is None or tree.library is not None:
            continue
        if not _is_mn_source(item):
            continue
        entry = _lookup(item.name)
        if entry is None:
            continue
        if entry.tree_setup is not None:
            _run_tree_setup(entry.tree_setup, tree)
        if entry.node_setup is not None:
            # the node instance doesn't exist yet — capture where the drop is
            # happening while the operator context is available, and defer
            space = getattr(bpy.context, "space_data", None)
            target = getattr(space, "edit_tree", None)
            obj = getattr(bpy.context, "object", None)
            _pending.append(
                (
                    tree.name,
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
    global _in_setup
    _in_setup = True
    try:
        setup(tree)
    except Exception:
        logger.exception(f"tree_setup failed for imported node group {tree.name!r}")
    finally:
        _in_setup = False


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
    for tree_name, target_name, object_name in entries:
        tree = bpy.data.node_groups.get(tree_name)
        if tree is None:
            continue
        entry = _lookup(tree.name)
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
                    or node.node_tree != tree
                    or node.get(_MARKER)
                ):
                    continue
                # a node in a nested group has no modifier of its own, so fall
                # back to the object that was active when the asset was dropped
                obj = _object_using_tree(host) or fallback_obj
                try:
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


def _style_material_default(tree: bpy.types.GeometryNodeTree) -> None:
    """
    Give a style's Material input the default MN material.

    Setting the default on the tree interface means every node instance created
    from it — including the one about to be created by the asset drop — starts
    with the material assigned, matching what the old custom add operator did.
    """
    from .material import add_all_materials

    socket = next(
        (
            s
            for s in tree.interface.items_tree
            if s.item_type == "SOCKET"
            and s.in_out == "INPUT"
            and s.socket_type == "NodeSocketMaterial"
        ),
        None,
    )
    if socket is None or socket.default_value is not None:
        return
    socket.default_value = add_all_materials()["MN Default"]


register_on_add(
    lambda name: name.startswith("Style "), tree_setup=_style_material_default
)
