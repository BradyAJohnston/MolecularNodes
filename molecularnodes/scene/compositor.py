from contextlib import contextmanager
from typing import Iterator, NamedTuple
import bpy
from bpy.types import CompositorNodeTree
from nodebpy import compositor as c
from nodebpy.builder import ColorSocket, TreeBuilder

annotations_image = "mn_annotations"


class ResetCompositorSockets(NamedTuple):
    """The sockets to build between, as returned by `CompositorTree.reset()`."""

    image: ColorSocket
    output: ColorSocket


def _annotations_image() -> bpy.types.Image:
    """The placeholder image annotations are drawn into, created on demand."""
    image = bpy.data.images.get(annotations_image)
    if image is None:
        image = bpy.data.images.new(annotations_image, 1, 1)
    return image


class CompositorTree(TreeBuilder[CompositorNodeTree]):
    """Builder for the scene's compositor node tree.

    Mirrors :class:`~molecularnodes.entities.base.MolecularTree`: use it as a
    context manager to append nodes, or :meth:`reset` to start from a clean
    ``Render Layers -> output`` graph and build the post-processing chain
    yourself with ``nodebpy.compositor`` nodes.

    ```python
    from nodebpy import compositor as c

    with canvas.compositor.reset() as (image, output):
        image >> c.Glare.bloom() >> output
    canvas.compositor.add_annotations()
    ```

    Execution settings for the render-time compositor are exposed as
    properties: :attr:`device`, :attr:`precision`, :attr:`denoise_device`,
    :attr:`denoise_preview_quality` and :attr:`denoise_final_quality`.

    ```python
    # on a headless machine without a GPU
    canvas.compositor.device = "CPU"
    ```

    Notes
    -----
    The surface is intentionally minimal — raw ``nodebpy`` node construction. We
    may later add higher-level convenience wrappers (e.g. ``add_glare()``,
    ``add_vignette()``) that build common effects in one call, but for now
    effects are composed by the user directly.
    """

    def __init__(self, scene: bpy.types.Scene) -> None:
        self._scene = scene
        if scene.compositing_node_group is None:
            scene.compositing_node_group = bpy.data.node_groups.new(
                "Compositor Nodes", "CompositorNodeTree"
            )
        super().__init__(scene.compositing_node_group)

    def _wrap(self, socket: bpy.types.NodeSocket) -> ColorSocket:
        """Wrap an existing Blender socket as a socket bound to this tree."""
        wrapped = ColorSocket(socket)
        wrapped._tree = self
        return wrapped

    def _render_layers_node(self) -> bpy.types.Node:
        """The Render Layers node feeding the tree, adding it if it is missing."""
        for node in self.tree.nodes:
            if node.bl_idname == "CompositorNodeRLayers":
                return node
        return self.tree.nodes.new("CompositorNodeRLayers")

    @property
    def device(self) -> str:
        """Device the compositor executes on: ``"CPU"`` or ``"GPU"``.

        Blender defaults to ``"GPU"``. Set to ``"CPU"`` on headless machines
        without a GPU, where GPU compositing aborts the render.
        """
        return self._scene.render.compositor_device

    @device.setter
    def device(self, value: str) -> None:
        self._scene.render.compositor_device = value.upper()  # ty: ignore[invalid-assignment]

    @property
    def precision(self) -> str:
        """Precision the compositor executes at: ``"AUTO"`` or ``"FULL"``.

        ``"AUTO"`` uses reduced precision for final renders; ``"FULL"`` always
        uses full precision.
        """
        return self._scene.render.compositor_precision

    @precision.setter
    def precision(self, value: str) -> None:
        self._scene.render.compositor_precision = value.upper()  # ty: ignore[invalid-assignment]

    @property
    def denoise_device(self) -> str:
        """Device Denoise nodes execute on: ``"AUTO"``, ``"CPU"`` or ``"GPU"``."""
        return self._scene.render.compositor_denoise_device

    @denoise_device.setter
    def denoise_device(self, value: str) -> None:
        self._scene.render.compositor_denoise_device = value.upper()  # ty: ignore[invalid-assignment]

    @property
    def denoise_preview_quality(self) -> str:
        """Denoise node quality in preview renders: ``"HIGH"``, ``"BALANCED"``
        or ``"FAST"``."""
        return self._scene.render.compositor_denoise_preview_quality

    @denoise_preview_quality.setter
    def denoise_preview_quality(self, value: str) -> None:
        self._scene.render.compositor_denoise_preview_quality = value.upper()  # ty: ignore[invalid-assignment]

    @property
    def denoise_final_quality(self) -> str:
        """Denoise node quality in final renders: ``"HIGH"``, ``"BALANCED"``
        or ``"FAST"``."""
        return self._scene.render.compositor_denoise_final_quality

    @denoise_final_quality.setter
    def denoise_final_quality(self, value: str) -> None:
        self._scene.render.compositor_denoise_final_quality = value.upper()  # ty: ignore[invalid-assignment]

    @property
    def image(self) -> ColorSocket:
        """The rendered image to build the compositor chain from.

        ```python
        with canvas.compositor as tree:
            tree.image >> c.Glare() >> tree.output
        ```
        """
        return self._wrap(self._render_layers_node().outputs["Image"])

    @property
    def output(self) -> ColorSocket:
        """The final image output of the tree, adding it if it is missing.

        ```python
        with canvas.compositor as tree:
            tree.clear() # remove all nodes and interface items
            tree.image >> c.Glare() >> tree.output
        ```
        """
        for item in self.tree.interface.items_tree:
            if (
                item.item_type == "SOCKET"
                and item.in_out == "OUTPUT"
                and item.socket_type == "NodeSocketColor"
            ):
                return self._wrap(self._output_node().inputs[item.identifier])
        return self.outputs.color("Image")

    @contextmanager
    def reset(self) -> Iterator[ResetCompositorSockets]:
        """Clear the tree back to a default state and build within it.

        Discards the existing tree — including any annotation overlay, which can
        be restored with :meth:`add_annotations`. Use ``with canvas.compositor``
        instead to append to the existing tree.

        ```python
        with canvas.compositor.reset() as (image, output):
            image >> c.Glare.bloom(strength=2.0) >> output
        ```
        """
        with self:
            self.clear()
            render = c.RenderLayers()
            output = self.outputs.color("Image")
            # default passthrough so an empty reset still renders the raw image
            render.o.image >> output
            yield ResetCompositorSockets(render.o.image, output)

    def clear(self) -> None:
        self.tree.nodes.clear()
        if self.tree.interface:
            self.tree.interface.clear()

    def add_annotations(self) -> None:
        """Composite Molecular Nodes annotations on top of the current output.

        Alpha-composites the ``mn_annotations`` image over whatever currently
        feeds the output. Annotations are opt-in: :meth:`reset` clears them, so
        call this after building a custom compositor chain to draw them on top.
        """
        with self:
            output = self.output
            links = output.socket.links
            source = (
                links[0].from_socket
                if links
                else self._render_layers_node().outputs["Image"]
            )
            annotations = c.Image(image=_annotations_image())
            c.AlphaOver(self._wrap(source), annotations, 1.0) >> output


def setup_compositor(scene: bpy.types.Scene) -> CompositorTree:
    """Prepare the scene compositor with the default annotation overlay."""
    # lock interface when rendering
    scene.render.use_lock_interface = True
    tree = CompositorTree(scene)
    with tree.reset():
        pass
    tree.add_annotations()
    return tree
