"""
Update the inputs of MolecularNodes node groups from an entity on frame changes.

Some node groups take values that change along a trajectory as plain inputs,
such as the box lengths and angles of the periodic nodes. Node groups can't
read these from the universe themselves, so the entity sets them on every
frame change: an update function is registered for node groups with the
``updater`` decorator, and ``update_nodes`` calls it for every instance of
those groups in a tree.

Only nodes directly in the tree are updated, not ones inside other groups,
which get their values through the group's own inputs.

Updates must be cheap, as they run on every frame change, and should do
nothing when the entity has nothing to give (no box, an ``Update`` input
turned off, ...).
"""

import re
from typing import TYPE_CHECKING, Any, Callable, Sequence, cast
import bpy
import numpy as np
from nodebpy.builder import AssetGeometryGroup
from .geometry import (
    PeriodicArray,
    PeriodicBox,
    PeriodicCell,
    PeriodicOffset,
    PeriodicUnwrap,
)

if TYPE_CHECKING:
    from ..entities.molecule.base import Molecule

NodeUpdate = Callable[[bpy.types.GeometryNodeGroup, Any], None]

# a node group appended again gets a ".001" suffix, but is the same group
_DUPLICATE_SUFFIX = re.compile(r"\.\d{3}$")

_registry: dict[str, list[NodeUpdate]] = {}


def updater(
    nodes: Sequence[type[AssetGeometryGroup]],
) -> Callable[[NodeUpdate], NodeUpdate]:
    """
    Register the decorated function as an update for node groups, run on every
    frame change.

    The function is called with each instance of the node groups in an
    entity's tree, and the entity.

    Parameters
    ----------
    nodes : Sequence[type[AssetGeometryGroup]]
        The node group classes to update, such as ``[PeriodicBox]``.

    Examples
    --------
    ```python
    @updater([PeriodicBox, PeriodicArray])
    def _periodic_box(node, entity):
        ...
    ```
    """

    def register(update: NodeUpdate) -> NodeUpdate:
        for node in nodes:
            _registry.setdefault(node._name, []).append(update)
        return update

    return register


def update_nodes(tree: bpy.types.NodeTree, entity: Any) -> None:
    """Run the registered updates for every node group instance in a tree."""
    for node in tree.nodes:
        if not isinstance(node, bpy.types.GeometryNodeGroup) or node.node_tree is None:
            continue
        name = _DUPLICATE_SUFFIX.sub("", node.node_tree.name)
        for update in _registry.get(name, []):
            update(node, entity)


@updater([PeriodicBox, PeriodicArray, PeriodicOffset, PeriodicUnwrap, PeriodicCell])
def _periodic_box(node: bpy.types.GeometryNodeGroup, entity: "Molecule") -> None:
    """Set the box lengths and angles from the current frame while Update is on."""
    if not cast(bpy.types.NodeSocketBool, node.inputs["Update"]).default_value:
        return
    dimensions = entity.universe.trajectory.ts.dimensions
    if dimensions is None:
        return
    # lengths from Å to world units, angles stay in degrees
    box = np.asarray(dimensions, dtype=float)
    values = np.concatenate([box[:3] * entity.world_scale, box[3:]])
    for name, value in zip(["a", "b", "c", "alpha", "beta", "gamma"], values):
        cast(bpy.types.NodeSocketFloat, node.inputs[name]).default_value = value
