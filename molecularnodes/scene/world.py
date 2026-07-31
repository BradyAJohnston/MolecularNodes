from contextlib import contextmanager
from typing import Iterator, Tuple
import bpy
from bpy.types import ShaderNodeTree
from nodebpy import shader as s
from nodebpy.builder import ShaderSocket, TreeBuilder

mn_world_shader_node_name = "MN_world_shader"


class WorldTree(TreeBuilder[ShaderNodeTree]):
    """Builder for the scene world's shader node tree (lighting & background).

    Replaces reaching directly into the ``MN_world_shader`` node. Use the
    convenience properties for the common controls, or :meth:`reset` to build a
    world shader from scratch with ``nodebpy.shader`` nodes.

    >>> canvas.world.hdri_strength = 1.5
    >>> canvas.world.background = (0.0, 0.5, 0.5, 1.0)

    >>> from nodebpy import shader as s
    >>> with canvas.world.reset() as surface:
    ...     s.Background(color=(0.05, 0.05, 0.05, 1.0), strength=1.0) >> surface
    """

    def __init__(self, scene: bpy.types.Scene) -> None:
        self._scene = scene
        world = scene.world
        if world is None:
            world = bpy.data.worlds.new("World")
            scene.world = world
        if world.node_tree is None:
            world.use_nodes = True
        super().__init__(world.node_tree)

    @property
    def _shader(self) -> bpy.types.Node:
        """The ``MN_world_shader`` node backing the convenience properties."""
        node = self.tree.nodes.get(mn_world_shader_node_name)
        if node is None:
            raise ValueError(
                f"'{mn_world_shader_node_name}' node not found — it is removed by "
                "`reset()`. Rebuild the world shader or reload the template to use "
                "the `background`/`hdri_strength` shortcuts."
            )
        return node

    @property
    def background(self) -> Tuple[float, float, float, float]:
        """World background color as RGBA in the range [0, 1]."""
        return tuple(self._shader.inputs[3].default_value)

    @background.setter
    def background(self, value: Tuple[float, float, float, float]) -> None:
        self._shader.inputs[3].default_value = value

    @property
    def hdri_strength(self) -> float:
        """Strength multiplier for the world lighting."""
        return self._shader.inputs[1].default_value

    @hdri_strength.setter
    def hdri_strength(self, value: float) -> None:
        self._shader.inputs[1].default_value = value

    @contextmanager
    def reset(self) -> Iterator[ShaderSocket]:
        """Clear the world shader and build within it.

        Discards the existing shader — including ``MN_world_shader``, so the
        `background`/`hdri_strength` shortcuts stop working until it is rebuilt.
        Yields the ``Surface`` input of the world output to build towards.

        >>> with canvas.world.reset() as surface:
        ...     s.Background(strength=1.0) >> surface
        """
        with self:
            self.clear()
            world_output = s.WorldOutput(is_active_output=True)
            yield world_output.i.surface

    def clear(self) -> None:
        self.tree.nodes.clear()
