import math
from collections.abc import Generator
from contextlib import contextmanager
from typing import Any, Literal, NamedTuple, cast
import bpy
from bpy.types import CompositorNodeTree
from nodebpy import compositor as c
from nodebpy.builder import ColorSocket, TreeBuilder
from ..nodes.compositor import ConeShadow, Illustrate
from ..session import get_session

annotations_image = "mn_annotations"

# world units per Angstrom that Molecular Nodes imports structures at
WORLD_SCALE = 0.1


def add_view_layer_aov(
    view_layer: bpy.types.ViewLayer, name: str, type: Literal["VALUE", "COLOR"]
) -> bpy.types.AOV:
    """Add a named AOV pass to the view layer, returning the existing one if
    the name is already taken."""
    aov = view_layer.aovs.get(name)
    if aov is None:
        aov = view_layer.aovs.add()
        aov.name = name
    aov.type = type
    return aov


def _links(socket: bpy.types.NodeSocket) -> bpy.types.NodeLinks:
    """The links of a socket, which the bpy stubs type as optional."""
    links = socket.links
    assert links is not None
    return links


def _entity_corner_depths(scene: bpy.types.Scene) -> list[float]:
    """Camera depths, in world units, of the bounding-box corners of every
    entity object in the scene (evaluated, so styles count)."""
    from mathutils import Vector

    camera = scene.camera
    if camera is None:
        return []
    depsgraph = bpy.context.evaluated_depsgraph_get()
    origin = camera.matrix_world.translation
    forward = camera.matrix_world.to_3x3() @ Vector((0.0, 0.0, -1.0))
    forward.normalize()
    depths = []
    for entity in get_session().entities.values():
        try:
            obj = entity.object
        except Exception:
            continue
        if obj is None:
            continue
        evaluated = obj.evaluated_get(depsgraph)
        for corner in evaluated.bound_box:
            point = evaluated.matrix_world @ Vector(corner)
            depths.append((point - origin).dot(forward))
    return depths


def _depth_range(scene: bpy.types.Scene) -> tuple[float, float]:
    """(near, far) camera depth of the entities in the scene, falling back to
    the camera clipping range."""
    depths = _entity_corner_depths(scene)
    if depths and max(depths) > min(depths):
        return max(min(depths), 0.0), max(depths)
    camera = scene.camera
    if camera is None:
        return 0.0, 100.0
    data = cast(bpy.types.Camera, camera.data)
    return data.clip_start, data.clip_end


def _pixel_size_angstrom(scene: bpy.types.Scene) -> float:
    """Width of one rendered pixel in Angstrom at the entities in the scene:
    from the orthographic scale, or the field of view at the mid-depth of the
    entities for a perspective camera."""
    camera = scene.camera
    render = scene.render
    pixels = max(render.resolution_x, render.resolution_y)
    pixels *= render.resolution_percentage / 100.0
    if camera is None or pixels <= 0:
        return 1.0
    data = cast(bpy.types.Camera, camera.data)
    if data.type == "ORTHO":
        width = data.ortho_scale
    else:
        near, far = _depth_range(scene)
        width = 2.0 * (near + far) / 2.0 * math.tan(data.angle / 2.0)
    return width / pixels / WORLD_SCALE


# shading source for `CompositorTree.illustrate()` -> (view layer pass
# attribute, Render Layers output accessor); "cone" is computed from the depth
# pass by a Cone Shadow node instead
_SHADING_SOURCES = {
    "cone": ("use_pass_z", "depth"),
    "ao": ("use_pass_ambient_occlusion", "ambient_occlusion"),
    "shadow": ("use_pass_shadow", "shadow"),
}


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
        assert scene.compositing_node_group is not None
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
        interface = self.tree.interface
        assert interface is not None
        for item in interface.items_tree:
            if (
                isinstance(item, bpy.types.NodeTreeInterfaceSocket)
                and item.in_out == "OUTPUT"
                and item.socket_type == "NodeSocketColor"
            ):
                return self._wrap(self._output_node().inputs[item.identifier])
        return self.outputs.color("Image")

    @contextmanager
    def reset(self) -> Generator[ResetCompositorSockets]:
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
            links = _links(output.socket)
            source = links[0].from_socket if links else None
            if source is None:
                source = self._render_layers_node().outputs["Image"]
            annotations = c.Image(image=_annotations_image())
            c.AlphaOver(self._wrap(source), annotations, 1.0) >> output

    def _remove_composites(self) -> None:
        """Remove any ``Illustrate`` node already in the tree, along with the
        ``Cone Shadow`` feeding it, handing its consumers back the rendered
        image, so :meth:`illustrate` replaces rather than stacks. Call inside
        ``with self``."""
        rl_node = self._render_layers_node()
        for node in list(self.tree.nodes):
            if (
                not isinstance(node, bpy.types.CompositorNodeGroup)
                or node.node_tree is None
            ):
                continue
            if node.node_tree.name != Illustrate._asset_name:
                continue
            for link in _links(node.inputs["Shadow"]):
                upstream = link.from_node
                if (
                    isinstance(upstream, bpy.types.CompositorNodeGroup)
                    and upstream.node_tree is not None
                    and upstream.node_tree.name == ConeShadow._asset_name
                ):
                    self.tree.nodes.remove(upstream)
            consumers = [
                to_socket
                for output in node.outputs
                for link in _links(output)
                if (to_socket := link.to_socket) is not None
            ]
            for to_socket in consumers:
                self.tree.links.new(rl_node.outputs["Image"], to_socket)
            self.tree.nodes.remove(node)

    def _insert_after_render(self, node) -> None:
        """Re-route every consumer of the rendered image to read from ``node``
        instead (the annotation overlay, or the output itself), or feed the
        output when nothing consumed it. Call inside ``with self``."""
        rl_node = self._render_layers_node()
        consumers = [
            to_socket
            for link in _links(rl_node.outputs["Image"])
            if link.to_node != node.node and (to_socket := link.to_socket) is not None
        ]
        if consumers:
            for to_socket in consumers:
                self.tree.links.new(node.node.outputs[0], to_socket)
        else:
            node >> self.output

    def illustrate(
        self,
        *,
        shading: Literal["cone", "ao", "shadow"] | None = "cone",
        outline: bool = True,
        outline_normals: bool = False,
        **inputs,
    ) -> Illustrate:
        """Insert an ``Illustrate`` node between Render Layers and whatever
        currently consumes the rendered image.

        Enables the render passes the node reads (depth, plus the normal pass
        for ``outline_normals`` and the chosen shading pass) on the view
        layer, adds the node group with its inputs linked to the Render Layers
        node, and re-routes every existing consumer of the rendered image (the
        annotation overlay, or the output itself) to read from it instead.
        Calling it again replaces the node.

        With ``shading="cone"`` the shading comes from a ``Cone Shadow`` node,
        David Goodsell's conical soft shadow computed from the depth pass
        alone, whose ``Pixel Size`` is set from the camera. It needs no
        lighting, so pair it with the ``Flat`` material
        (``mol.add_style("cartoon", material=mn.material.Flat())``) for
        Illustrate-style renders.

        ```python
        node = canvas.compositor.illustrate(shading="cone", outline=True)
        node.i.outline_size.default_value = 3
        node.i.outline_color.default_value = (0.1, 0.1, 0.2, 1.0)
        ```

        Parameters
        ----------
        shading : {"cone", "ao", "shadow"} or None, default "cone"
            What darkens the image: a ``Cone Shadow`` from the depth pass, the
            ambient occlusion pass, the shadow pass (EEVEE only; Cycles has no
            Shadow pass), or ``None`` for no shading at all.
        outline : bool, default True
            Draw lines where the depth changes sharply, from the ``Depth
            Outline`` node inside the group.
        outline_normals : bool, default False
            Also draw lines where the surface normal changes sharply; enables
            the normal pass.
        **inputs
            Any other input of the node by its Python name, for example
            ``outline_size=3`` or ``background=True``.

        Returns
        -------
        Illustrate
            The node handle; tweak its inputs afterwards through ``node.i``.
        """
        if shading is not None and shading not in _SHADING_SOURCES:
            raise ValueError(
                f"Unknown shading source {shading!r}; "
                f"expected one of {sorted(_SHADING_SOURCES)} or None."
            )
        view_layer = self._scene.view_layers[0]
        view_layer.use_pass_z = True
        if outline_normals:
            view_layer.use_pass_normal = True
        if shading is not None:
            pass_attr, accessor = _SHADING_SOURCES[shading]
            setattr(view_layer, pass_attr, True)

        with self:
            self._remove_composites()
            rl_node = self._render_layers_node()
            render = c.RenderLayers._from_node(rl_node)
            links: dict[str, Any] = {}
            if outline_normals:
                links["normal"] = render.o.normal
            if shading == "cone":
                links["shadow"] = ConeShadow(
                    depth=render.o.depth,
                    # the cone test compares against the depth pass, in world units
                    pixel_size=_pixel_size_angstrom(self._scene) * WORLD_SCALE,
                )
            elif shading is not None:
                if accessor not in {
                    s.name.lower().replace(" ", "_")
                    for s in rl_node.outputs
                    if s.enabled
                }:
                    raise ValueError(
                        f"The {shading!r} shading pass is not available with the "
                        f"{self._scene.render.engine} engine; the Shadow pass is EEVEE only."
                    )
                links["shadow"] = getattr(render.o, accessor)
            node = Illustrate(
                image=render.o.image,
                alpha=render.o.alpha,
                depth=render.o.depth,
                shading=shading is not None,
                outline=outline,
                outline_normals=outline_normals,
                **links,
                **inputs,
            )
            self._insert_after_render(node)
        return node


def setup_compositor(scene: bpy.types.Scene) -> CompositorTree:
    """Prepare the scene compositor with the default annotation overlay."""
    # lock interface when rendering
    scene.render.use_lock_interface = True
    tree = CompositorTree(scene)
    with tree.reset():
        pass
    tree.add_annotations()
    return tree
