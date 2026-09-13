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

The handler fires while the add operator is still importing the datablock,
before the node instance exists in the tree, so it only records what was
imported and where. The setup itself runs on each new node instance via a
one-shot timer, which fires right after the operator has created the node.
Node instances are local even when their node group is linked, so their input
socket values and custom properties are always writable.

Callbacks must be idempotent — they may run more than once for the same
datablock and are also invoked on nodes created programmatically, so they
should no-op when there is nothing to do (a value already set, no matching
entity, ...).
"""

import logging
import re
from pathlib import Path
from typing import Callable
import bpy
from bpy.app.handlers import persistent
from ..assets import MN_DATA_FILE

logger = logging.getLogger(__name__)

NodeSetup = Callable[[bpy.types.GeometryNodeGroup, bpy.types.Object | None], None]

# stamped onto a node instance once its setup has been attempted, so a later
# add of the same asset doesn't re-run setup on pre-existing nodes
_MARKER = "mn_on_add_done"
_DUPLICATE_SUFFIX = re.compile(r"\.\d{3}$")

_registry: list[tuple[Callable[[str], bool], NodeSetup]] = []
# (imported group name, drop-target tree name, context object name) awaiting
# setup once the add operator has created the node instance
_pending: list[tuple[str, str | None, str | None]] = []


def register_on_add(match: str | Callable[[str], bool], setup: NodeSetup) -> None:
    """
    Register a setup callback for a node group asset.

    Parameters
    ----------
    match : str | Callable[[str], bool]
        Either the exact name of the asset's node group, or a predicate called
        with the (duplicate-suffix stripped) name of each imported group.
    setup : Callable
        Called with each newly added node instance and the object whose modifier
        tree it was added to (or ``None`` when that can't be determined).
    """
    if not callable(match):
        name = match
        match = lambda n: n == name  # noqa: E731
    _registry.append((match, setup))


def _lookup(name: str) -> NodeSetup | None:
    name = _DUPLICATE_SUFFIX.sub("", name)
    for predicate, setup in _registry:
        if predicate(name):
            return setup
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
    Queue registered setup for MN node group assets after they are imported.

    Registered on ``bpy.app.handlers.blend_import_post``, which fires for asset
    drag-and-drop as well as any other link/append/pack, including adds that
    reuse an already imported group.
    """
    schedule = False
    for item in ctx.import_items:
        # only directly requested node trees, not their dependencies
        if item.id_type != "NODETREE" or item.import_info:
            continue
        tree = item.id
        if tree is None or not _is_mn_source(item):
            continue
        if _lookup(item.name) is None:
            continue
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


def _object_using_tree(tree: bpy.types.NodeTree) -> bpy.types.Object | None:
    """Find the object using this tree as a geometry nodes modifier, if any."""
    for obj in bpy.data.objects:
        for mod in obj.modifiers:
            if mod.type == "NODES" and mod.node_group == tree:
                return obj
    return None


def _process_pending() -> None:
    """Run deferred setup callbacks on newly added node instances."""
    entries, _pending[:] = list(_pending), []
    for group_name, target_name, object_name in entries:
        setup = _lookup(group_name)
        if setup is None:
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
                    setup(node, obj)
                except Exception:
                    logger.exception(
                        f"setup failed for node {node.name!r} in {host.name!r}"
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
    from ..material import add_all_materials

    socket = next((s for s in node.inputs if s.bl_idname == "NodeSocketMaterial"), None)
    if socket is None or socket.default_value is not None:
        return
    socket.default_value = add_all_materials()["Default"]


register_on_add(lambda name: name.startswith("Style "), _style_node_material)
