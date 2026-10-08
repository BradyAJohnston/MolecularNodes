import math
from contextlib import contextmanager
from typing import Iterable, Iterator, Literal, NamedTuple
import bpy
from bpy.types import CompositorNodeTree
from nodebpy import compositor as c
from nodebpy.builder import ColorSocket, TreeBuilder
from .. import material as _material
from ..nodes.compositor import CompositeIllustrate, CompositeIllustrative
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
    aov.type = type  # ty: ignore[invalid-assignment]
    return aov


def _mn_materials() -> list[bpy.types.Material]:
    """The materials built from the Molecular Nodes presets ("Flat",
    "Flat.001", ...), which are the ones the styles use."""
    prefixes = tuple(_material.MATERIAL_NAMES)
    return [
        mat
        for mat in bpy.data.materials
        if mat.library is None
        and mat.node_tree is not None
        and (mat.name in prefixes or mat.name.rsplit(".", 1)[0] in prefixes)
    ]


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
    return camera.data.clip_start, camera.data.clip_end


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
    data = camera.data
    if data.type == "ORTHO":
        width = data.ortho_scale
    else:
        near, far = _depth_range(scene)
        width = 2.0 * (near + far) / 2.0 * math.tan(data.angle / 2.0)
    return width / pixels / WORLD_SCALE


# the composites `illustrative()` and `illustrate()` insert; a new call replaces
# whichever of them is already in the tree
_COMPOSITE_ASSETS = {CompositeIllustrative._asset_name, CompositeIllustrate._asset_name}

# shading source for `CompositorTree.illustrative()` -> (view layer pass
# attribute, Render Layers output accessor, group menu item)
_SHADING_SOURCES = {
    "ao": ("use_pass_ambient_occlusion", "ambient_occlusion", "Ambient Occlusion"),
    "shadow": ("use_pass_shadow", "shadow", "Shadow"),
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

    def illustrative(
        self,
        *,
        outline: bool = True,
        shading: Literal["ao", "shadow"] | None = "ao",
        flat: bool = False,
        **inputs,
    ) -> CompositeIllustrative:
        """Insert a ``Composite Illustrative`` node between Render Layers and
        whatever currently consumes the rendered image.

        Enables the render passes the node reads (depth, normal, diffuse
        colour and the chosen shading pass) on the view layer, adds the node
        group with its inputs linked to the Render Layers node, and re-routes
        every existing consumer of the rendered image (the annotation overlay,
        or the output itself) to read from it instead. Calling it again, or
        calling :meth:`illustrate`, replaces the node.

        ```python
        node = canvas.compositor.illustrative(outline=True, shading="ao")
        node.i.outline_size.default_value = 3
        node.i.outline_color.default_value = (0.1, 0.1, 0.2, 1.0)
        ```

        Parameters
        ----------
        outline : bool, default True
            Draw lines where the depth (and optionally the normal) changes
            sharply, from the ``Composite Outline Mask`` node inside the group.
        shading : {"ao", "shadow"} or None, default "ao"
            Which pass darkens the base colour: the ambient occlusion pass, the
            shadow pass (EEVEE only; Cycles has no Shadow pass), or ``None``
            for no shading at all.
        flat : bool, default False
            Shade the flat ``Diffuse Color`` pass, which ignores lighting and
            materials, instead of the rendered image.
        **inputs
            Any other input of the node by its Python name, for example
            ``outline_size=3`` or ``depth_threshold=4.0`` (Angstrom).

        Returns
        -------
        CompositeIllustrative
            The node handle; tweak its inputs afterwards through ``node.i``.
        """
        if shading is not None and shading not in _SHADING_SOURCES:
            raise ValueError(
                f"Unknown shading source {shading!r}; "
                f"expected one of {sorted(_SHADING_SOURCES)} or None."
            )
        view_layer = self._scene.view_layers[0]
        view_layer.use_pass_z = True
        view_layer.use_pass_normal = True
        view_layer.use_pass_diffuse_color = True
        links: dict[str, object] = {}
        if shading is not None:
            pass_attr, accessor, menu_item = _SHADING_SOURCES[shading]
            setattr(view_layer, pass_attr, True)
            inputs.setdefault("shading_source", menu_item)

        with self:
            self._remove_composites()
            rl_node = self._render_layers_node()
            render = c.RenderLayers._from_node(rl_node)
            if shading is not None:
                if accessor not in {
                    s.name.lower().replace(" ", "_")
                    for s in rl_node.outputs
                    if s.enabled
                }:
                    raise ValueError(
                        f"The {shading!r} shading pass is not available with the "
                        f"{self._scene.render.engine} engine; the Shadow pass is EEVEE only."
                    )
                links[accessor] = getattr(render.o, accessor)
            node = CompositeIllustrative(
                image=render.o.image,
                alpha=render.o.alpha,
                depth=render.o.depth,
                normal=render.o.normal,
                diffuse_color=render.o.diffuse_color,
                base_color="Diffuse Color" if flat else "Image",
                shading=shading is not None,
                outline=outline,
                **links,
                **inputs,
            )
            self._insert_after_render(node)
        return node

    def _remove_composites(self) -> None:
        """Remove any ``Composite Illustrative`` / ``Composite Illustrate``
        node already in the tree, handing its consumers back the rendered
        image, so the helpers replace rather than stack. Call inside
        ``with self``."""
        rl_node = self._render_layers_node()
        for node in list(self.tree.nodes):
            if node.bl_idname != "CompositorNodeGroup" or node.node_tree is None:
                continue
            if node.node_tree.name not in _COMPOSITE_ASSETS:
                continue
            consumers = [
                link.to_socket for output in node.outputs for link in output.links
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
            link.to_socket
            for link in rl_node.outputs["Image"].links
            if link.to_node != node.node
        ]
        if consumers:
            for to_socket in consumers:
                self.tree.links.new(node.node.outputs[0], to_socket)
        else:
            node >> self.output

    def illustrate(
        self,
        *,
        shadow: bool = True,
        fog: bool = False,
        contour: bool = True,
        chain_outline: bool = False,
        residue_outline: bool = False,
        flat: bool = False,
        materials: Iterable[bpy.types.Material] | None = None,
        **inputs,
    ) -> CompositeIllustrate:
        """Insert a ``Composite Illustrate`` node, David Goodsell's Illustrate
        composition, between Render Layers and whatever consumes the image.

        Illustrate has no lighting model: form comes from a conical soft
        shadow computed from the depth pass, depth fog and outlines, on flat
        per-atom colour. Use it with the ``Flat`` material
        (``mol.add_style("cartoon", material=mn.material.Flat())``), which
        renders the atom colours unlit; the node reads the ``Diffuse Color``
        pass by default, so any material gives flat colour, but the Flat
        material also keeps the rendered image consistent with it.

        Enables the depth and diffuse colour passes, and for the chain and
        residue outlines declares ``chain_id`` / ``res_id`` AOV passes on the
        view layer and adds an ``AOV Output`` writing that attribute (offset
        by one so the first ID differs from the empty background) to every
        Molecular Nodes material in the file that does not write it yet. The
        cone shadow's ``Pixel Size`` is set from the camera and the fog's
        ``Near``/``Far`` from the camera depth of the entities' bounds, unless
        given in ``inputs``.

        ```python
        node = canvas.compositor.illustrate(shadow=True, contour=True, chain_outline=True)
        node.i.max_darkening.default_value = 0.3
        ```

        Parameters
        ----------
        shadow : bool, default True
            Conical soft shadow (``Composite Cone Shadow``).
        fog : bool, default False
            Depth fog towards ``fog_color`` (``Composite Depth Fog``); set
            ``back_fog`` below 1 to see it.
        contour : bool, default True
            Contour outlines on depth steps (``Composite Contour Outline``).
        chain_outline : bool, default False
            Outlines between chains from a ``chain_id`` AOV.
        residue_outline : bool, default False
            Outlines between residues from a ``res_id`` AOV.
        flat : bool, default False
            Use the ``Diffuse Color`` pass as the base colour instead of the
            rendered image: flat, unlit colour from any material. With the
            Flat material the rendered image is already flat, and emission
            shaders leave the diffuse pass empty, so the image is the default.
        materials : iterable of bpy.types.Material, optional
            Materials to add the AOV outputs to; defaults to every material
            built from a Molecular Nodes preset.
        **inputs
            Any other input of the node by its Python name, for example
            ``radius=30.0`` or ``back_fog=0.4``.

        Returns
        -------
        CompositeIllustrate
            The node handle; tweak its inputs afterwards through ``node.i``.
        """
        view_layer = self._scene.view_layers[0]
        view_layer.use_pass_z = True
        view_layer.use_pass_diffuse_color = True
        aovs = {
            "chain_id": (chain_outline, "chain_id"),
            "res_id": (residue_outline, "residue_id"),
        }
        targets = list(materials) if materials is not None else _mn_materials()
        for aov, (enabled, _) in aovs.items():
            if not enabled:
                continue
            add_view_layer_aov(view_layer, aov, "VALUE")
            for mat in targets:
                if not _material.has_aov(mat, aov):
                    _material.add_aov(mat, aov, offset=1.0)
        inputs.setdefault("pixel_size", _pixel_size_angstrom(self._scene))
        if fog and not ({"near", "far"} <= set(inputs)):
            near, far = _depth_range(self._scene)
            inputs.setdefault("near", near)
            inputs.setdefault("far", far)

        with self:
            self._remove_composites()
            render = c.RenderLayers._from_node(self._render_layers_node())
            links = {
                key: getattr(render.o, aov)
                for aov, (enabled, key) in aovs.items()
                if enabled
            }
            node = CompositeIllustrate(
                image=render.o.image,
                alpha=render.o.alpha,
                depth=render.o.depth,
                diffuse_color=render.o.diffuse_color,
                base_color="Diffuse Color" if flat else "Image",
                shadow=shadow,
                fog=fog,
                contour_outline=contour,
                chain_outline=chain_outline,
                residue_outline=residue_outline,
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
