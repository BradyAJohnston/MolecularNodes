"""Render layers: render groups of styles separately and composite them in order.

Geometry is put on a layer with the ``Set Render Layer`` node, which stores
``mn_render_layer`` (layer + 1, so 0 is "every layer") on points and instances
for the selection. Each layer is a Blender view layer tagged through
``ViewLayer.mn_layer``; MN creates the view layers and the compositor stack, so
users only deal with a list of layers and ``layer=`` on styles.

Entity objects themselves are hidden in renders. One ``_ViewLayer`` object
sources their evaluated geometry with Object Info and runs it through
``Filter Render Layer``, which keeps the geometry of the layer being rendered.
The layer comes from one driver on the shared ``.MN_render_layer_current``
group, reading the plain ``["mn_layer_index"]`` property of the active view
layer: each view layer's depsgraph resolves it for itself, the expression is a
simple one that needs no Python auto-run, and a plain custom property resolves
without the add-on installed. The entity trees never depend on the view layer,
so they evaluate once per render and only the small filter tree differs.
"""

from __future__ import annotations
from collections.abc import Iterator
from typing import TYPE_CHECKING, cast
from uuid import uuid4
import bpy
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    FloatVectorProperty,
    IntProperty,
    StringProperty,
)

if TYPE_CHECKING:
    from .compositor import CompositorTree

CURRENT_GROUP = ".MN_render_layer_current"
INDEX_PROP = "mn_layer_index"
PROXY_NAME = "_ViewLayer"
PROXY_TREE = ".MN_ViewLayer"
SOURCE_PROP = "mn_render_layer_source"
BUILT_PROP = "mn_layers_built"
BACKGROUND_NODE = "MN Layers Background"
NO_LAYER = -1


def _current_nodes() -> list[bpy.types.FunctionNodeInputInt]:
    """The Integer node of every copy of the shared current-layer group."""
    nodes = []
    for group in bpy.data.node_groups:
        if not group.name.startswith(CURRENT_GROUP):
            continue
        for node in group.nodes:
            if node.bl_idname == "FunctionNodeInputInt":
                nodes.append(cast(bpy.types.FunctionNodeInputInt, node))
    return nodes


def _ensure_drivers() -> None:
    """Drive the current layer from the active view layer's index property."""
    for node in _current_nodes():
        tree = cast(bpy.types.NodeTree, node.id_data)
        path = f'nodes["{node.name}"].integer'
        anim = tree.animation_data
        if anim is not None and any(fc.data_path == path for fc in anim.drivers):
            continue
        driver = cast(bpy.types.FCurve, tree.driver_add(path)).driver
        assert driver is not None
        driver.type = "SCRIPTED"
        var = driver.variables.new()
        var.name = "idx"
        var.type = "CONTEXT_PROP"
        target = var.targets[0]
        target.context_property = "ACTIVE_VIEW_LAYER"
        target.data_path = f'["{INDEX_PROP}"]'
        driver.expression = "idx"


class MN_RenderLayerSettings(bpy.types.PropertyGroup):
    """Molecular Nodes render layer settings, stored on each view layer.

    ``blend``, ``outline`` and ``outline_normals`` decide which compositor
    nodes are built. The other values are what a rebuild starts from: once
    built, they live on the compositor nodes, which the panel edits directly.
    """

    uid: StringProperty(  # type: ignore
        name="ID",
        description="Names this layer's compositor nodes, so they are found "
        "after the view layer is renamed or moved",
    )
    is_layer: BoolProperty(  # type: ignore
        name="Render Layer",
        description="Composite this view layer as a Molecular Nodes render layer",
        default=False,
    )
    index: IntProperty(  # type: ignore
        name="Index",
        description="Layer index that Set Render Layer nodes assign. Higher layers "
        "are composited over lower ones",
        default=0,
        min=0,
    )
    opacity: FloatProperty(  # type: ignore
        name="Opacity",
        description="Opacity of this layer over the layers beneath it",
        default=1.0,
        min=0.0,
        max=1.0,
        subtype="FACTOR",
    )
    blend: EnumProperty(  # type: ignore
        name="Blend",
        description="How this layer is combined with the layers beneath it",
        items=(
            ("OVER", "Over", "Draw over the layers beneath, ignoring depth"),
            (
                "DEPTH",
                "Depth",
                "Combine by depth, so lower layers in front still show",
            ),
        ),
        default="OVER",
    )
    outline: BoolProperty(  # type: ignore
        name="Outline",
        description="Draw outlines from this layer's depth pass",
        default=True,
    )
    outline_color: FloatVectorProperty(  # type: ignore
        name="Outline Color",
        subtype="COLOR",
        size=4,
        min=0.0,
        max=1.0,
        default=(0.0, 0.0, 0.0, 1.0),
    )
    outline_size: IntProperty(  # type: ignore
        name="Outline Size",
        description="Line thickness in pixels",
        default=2,
        min=1,
        max=20,
    )
    outline_threshold: FloatProperty(  # type: ignore
        name="Depth Threshold",
        description="Depth gradient, in world units like the depth pass, above which a line is drawn",
        default=6.0,
        min=0.0,
    )
    outline_normals: BoolProperty(  # type: ignore
        name="Normal Lines",
        description="Also draw lines where the surface normal changes sharply",
        default=False,
    )


def _layer_index(view_layer: bpy.types.ViewLayer) -> int:
    settings = view_layer.mn_layer  # ty: ignore[unresolved-attribute]
    return settings.index if settings.is_layer else NO_LAYER


# --- manager ----------------------------------------------------------------


# settings that live on a compositor node once built: node role, input socket
LIVE_SETTINGS = {
    "opacity": ("Opacity", None),
    "outline_color": ("Outline", "Color"),
    "outline_size": ("Lines", "Size"),
    "outline_threshold": ("Lines", "Threshold"),
}


def _node_name(uid: str, role: str) -> str:
    return f"MN Layer {uid} {role}"


def _opacity_socket(node: bpy.types.Node) -> bpy.types.NodeSocket:
    # Alpha Over's factor for OVER, the alpha multiplier for DEPTH
    if node.bl_idname == "CompositorNodeAlphaOver":
        return node.inputs["Factor"]
    return node.inputs[1]


class RenderLayer:
    """One render layer: a view layer with Molecular Nodes layer settings.

    ``opacity`` and the outline colour, size and threshold read and write the
    compositor node inputs once the compositor is built, so the panel, the
    Python API and the compositor itself all edit the same value. Before that,
    and for the settings that change the compositor's structure, they read and
    write the ``mn_layer`` properties of the view layer.
    """

    def __init__(self, view_layer: bpy.types.ViewLayer) -> None:
        self.view_layer = view_layer

    @property
    def settings(self) -> MN_RenderLayerSettings:
        return self.view_layer.mn_layer  # ty: ignore[unresolved-attribute]

    @property
    def name(self) -> str:
        return self.view_layer.name

    @property
    def index(self) -> int:
        return self.settings.index

    def socket(self, name: str) -> bpy.types.NodeSocket | None:
        """The compositor input that holds the setting ``name``, if built."""
        if name not in LIVE_SETTINGS or not self.settings.uid:
            return None
        role, input = LIVE_SETTINGS[name]
        scene = cast(bpy.types.Scene, self.view_layer.id_data)
        # read only: the panel calls this while drawing
        tree = scene.compositing_node_group
        node = tree.nodes.get(_node_name(self.settings.uid, role)) if tree else None
        if node is None:
            return None
        return _opacity_socket(node) if input is None else node.inputs.get(input)

    def pull(self) -> None:
        """Copy the values from the compositor nodes back into the settings."""
        for name in LIVE_SETTINGS:
            socket = self.socket(name)
            if socket is not None:
                value = socket.default_value  # ty: ignore[unresolved-attribute]
                setattr(self.settings, name, value)

    def __getattr__(self, name: str):
        if name in MN_RenderLayerSettings.__annotations__:
            socket = self.socket(name)
            if socket is not None:
                return socket.default_value  # ty: ignore[unresolved-attribute]
            return getattr(self.settings, name)
        raise AttributeError(name)

    def __setattr__(self, name: str, value) -> None:
        if name in MN_RenderLayerSettings.__annotations__:
            setattr(self.settings, name, value)
            socket = self.socket(name)
            if socket is not None:
                socket.default_value = value  # ty: ignore[unresolved-attribute]
            elif name not in LIVE_SETTINGS:
                # changes which nodes exist: rebuild from a script straight away
                scene = cast(bpy.types.Scene, self.view_layer.id_data)
                RenderLayers(scene).sync()
        else:
            super().__setattr__(name, value)

    def __repr__(self) -> str:
        return f"<RenderLayer {self.index}: {self.name!r}>"


class RenderLayers:
    """The render layers of a scene, as ``canvas.layers``.

    ```python
    canvas = mn.Canvas()
    mol = mn.Molecule.fetch("6n2y")
    mol.add_style("cartoon")                     # layer 0, the base layer
    mol.add_style("surface", layer=1)
    canvas.layers.add("Surface", opacity=0.4, outline_size=3)
    canvas.snapshot("layers.png")
    ```

    Layer 0 is the scene's first view layer. Adding layers creates further
    view layers and rebuilds the compositor to stack them, bottom to top.
    """

    def __init__(self, scene: bpy.types.Scene) -> None:
        self.scene = scene

    def _view_layers(self) -> list[bpy.types.ViewLayer]:
        return sorted(
            (vl for vl in self.scene.view_layers if vl.mn_layer.is_layer),  # ty: ignore[unresolved-attribute]
            key=_layer_index,
        )

    def __iter__(self) -> Iterator[RenderLayer]:
        return (RenderLayer(vl) for vl in self._view_layers())

    def __len__(self) -> int:
        return len(self._view_layers())

    def __getitem__(self, key: int | str) -> RenderLayer:
        for layer in self:
            if key == layer.index or key == layer.name:
                return layer
        raise KeyError(key)

    def _ensure_base(self) -> None:
        base = self.scene.view_layers[0]
        settings = base.mn_layer  # ty: ignore[unresolved-attribute]
        if not settings.is_layer:
            settings.is_layer = True
            settings.index = 0
        if not settings.uid:
            settings.uid = uuid4().hex[:8]

    def add(self, name: str | None = None, **settings) -> RenderLayer:
        """Add a render layer on top of the existing ones.

        Parameters
        ----------
        name : str, optional
            Name of the view layer, ``"Layer <index>"`` by default.
        **settings
            Any layer setting: ``opacity``, ``blend``, ``outline``,
            ``outline_color``, ``outline_size``, ``outline_threshold``,
            ``outline_normals``.

        Returns
        -------
        RenderLayer
            The new layer. Pass its ``index`` (or name) to
            ``add_style(layer=...)``.
        """
        self._ensure_base()
        base = self.scene.view_layers[0]
        index = max(layer.index for layer in self) + 1
        view_layer = self.scene.view_layers.new(name or f"Layer {index}")
        # passes the base layer renders are needed by every layer
        for attr in ("use_pass_z", "use_pass_normal"):
            setattr(view_layer, attr, getattr(base, attr))
        layer = RenderLayer(view_layer)
        layer.settings.is_layer = True
        layer.settings.index = index
        layer.settings.uid = uuid4().hex[:8]
        for key, value in settings.items():
            if key not in MN_RenderLayerSettings.__annotations__:
                raise TypeError(f"Unknown render layer setting {key!r}")
            setattr(layer.settings, key, value)
        self.sync()
        return layer

    def remove(self, key: int | str) -> None:
        layer = self[key]
        if layer.view_layer == self.scene.view_layers[0]:
            raise ValueError("The base layer can't be removed")
        self.scene.view_layers.remove(layer.view_layer)
        self.sync()

    def swap_style_layers(self, a: int, b: int) -> None:
        """Swap the layer of every ``Set Render Layer`` node on layer ``a``
        with those on layer ``b``, in every node tree."""
        for tree in bpy.data.node_groups:
            if tree.bl_idname != "GeometryNodeTree":
                continue
            for node in tree.nodes:
                if not _is_set_node(node):
                    continue
                socket = cast(bpy.types.NodeSocketInt, node.inputs["Layer"])
                if socket.default_value == a:
                    socket.default_value = b
                elif socket.default_value == b:
                    socket.default_value = a

    def resolve(self, layer: int | str | RenderLayer) -> int:
        """The layer index for a layer given by index, name or handle."""
        if isinstance(layer, RenderLayer):
            return layer.index
        if isinstance(layer, int):
            return layer
        return self[layer].index

    def _signature(self) -> str:
        """What the compositor's structure depends on. Values that only set a
        node input (opacity, outline colour, ...) are left out."""
        layers = list(self)
        if len(layers) < 2:
            return ""
        parts = [
            (
                layer.settings.uid,
                layer.name,
                layer.index,
                layer.settings.blend,
                layer.settings.outline,
                layer.settings.outline_normals,
            )
            for layer in layers
        ]
        has_background = self.scene.mn_layers_background[3] > 0  # ty: ignore[unresolved-attribute]
        return repr((parts, has_background))

    @property
    def is_stale(self) -> bool:
        """Whether the layers changed in a way the compositor doesn't show
        until it is rebuilt."""
        return self.scene.get(BUILT_PROP, "") != self._signature()

    def sync(self, rebuild: bool = False) -> None:
        """Bring the hidden machinery in line with the layers: the index
        property on every view layer, the ``_ViewLayer`` object and its
        sources, and the driver. The compositor is rebuilt only when it is
        stale or ``rebuild`` is set, so edits made to its nodes survive."""
        for view_layer in self.scene.view_layers:
            view_layer[INDEX_PROP] = _layer_index(view_layer)
        _sync_proxy(self.scene, enabled=len(self) > 1)
        if rebuild or self.is_stale:
            self.composite()

    def composite(self) -> None:
        """Rebuild the compositor to stack the layers, bottom to top, followed
        by the annotation overlay. Values edited on the previous compositor's
        nodes carry over."""
        from .compositor import CompositorTree

        layers = list(self)
        tree = CompositorTree(self.scene)
        for layer in layers:
            layer.pull()
        if tree.tree is not None and BACKGROUND_NODE in tree.tree.nodes:
            node = tree.tree.nodes[BACKGROUND_NODE]
            self.scene.mn_layers_background = node.inputs["Background"].default_value  # ty: ignore[unresolved-attribute]
        self.scene[BUILT_PROP] = self._signature()
        if len(layers) < 2:
            with tree.reset():
                pass
            tree.add_annotations()
            return
        # upper layers must be transparent to stack, and film transparency is
        # per scene, so the background is added back by the compositor
        self.scene.render.film_transparent = True
        for layer in layers:
            vl = layer.view_layer
            vl.use_pass_z = True
            if layer.outline_normals:
                vl.use_pass_normal = True
        with tree.reset() as (_, output):
            result = _stack(tree, layers)
            background = tuple(self.scene.mn_layers_background)  # ty: ignore[unresolved-attribute]
            if background[3] > 0:
                from nodebpy import compositor as c

                result = c.AlphaOver(background=background, foreground=result)
                result.node.name = BACKGROUND_NODE
            result >> output
        # reset() adds its own Render Layers node; the stack makes one per layer
        for node in list(tree.tree.nodes):
            if node.bl_idname == "CompositorNodeRLayers" and not any(
                socket.is_linked for socket in node.outputs
            ):
                tree.tree.nodes.remove(node)
        tree.add_annotations()

    def background_socket(self) -> bpy.types.NodeSocket | None:
        """The compositor input that holds the background colour, if built."""
        tree = self.scene.compositing_node_group
        node = tree.nodes.get(BACKGROUND_NODE) if tree else None
        return node.inputs["Background"] if node is not None else None

    @property
    def background(self) -> tuple[float, float, float, float]:
        """Colour the stacked layers are composited over. Alpha 0 keeps the
        render transparent; switching between transparent and not needs a
        rebuild."""
        socket = self.background_socket()
        if socket is not None:
            return tuple(socket.default_value)  # ty: ignore[unresolved-attribute]
        return tuple(self.scene.mn_layers_background)  # ty: ignore[unresolved-attribute]

    @background.setter
    def background(self, value) -> None:
        self.scene.mn_layers_background = value  # ty: ignore[unresolved-attribute]
        socket = self.background_socket()
        if socket is not None:
            socket.default_value = value  # ty: ignore[unresolved-attribute]
        self.sync()


def _layer_image(layer: RenderLayer):
    """The layer's rendered image with its outlines, and its depth pass."""
    from nodebpy import compositor as c
    from ..nodes.compositor import DepthOutline

    settings = layer.settings
    render = c.RenderLayers(layer=layer.name)
    image = render.o.image
    mask = None
    if settings.outline:
        mask = DepthOutline(
            depth=render.o.depth,
            size=settings.outline_size,
            threshold=settings.outline_threshold,
            normals=settings.outline_normals,
            normals_pass=render.o.normal if settings.outline_normals else None,
        )
        mask.node.name = _node_name(settings.uid, "Lines")
    return render, image, mask


def _stack(tree: CompositorTree, layers: list[RenderLayer]):
    from nodebpy import compositor as c
    from ..nodes.compositor import Outline

    def outline(layer: RenderLayer, image, mask):
        settings = layer.settings
        node = Outline(image=image, mask=mask, color=tuple(settings.outline_color))
        node.node.name = _node_name(settings.uid, "Outline")
        return node

    base = layers[0]
    render, result, mask = _layer_image(base)
    if mask is not None:
        result = outline(base, result, mask)
    depth = render.o.depth
    for layer in layers[1:]:
        settings = layer.settings
        render, image, mask = _layer_image(layer)
        if settings.blend == "DEPTH":
            alpha = render.o.alpha * settings.opacity
            opacity = alpha.node
            faded = c.SetAlpha(image, alpha=alpha)
            combined = c.DepthCombine(
                a=result, depth_a=depth, b=faded, depth_b=render.o.depth, use_alpha=True
            )
            result = combined.o.result
            depth = combined.o.depth
        else:
            result = opacity = c.AlphaOver(
                background=result, foreground=image, factor=settings.opacity
            )
            opacity = opacity.node
        opacity.name = _node_name(settings.uid, "Opacity")
        if mask is not None:
            result = outline(layer, result, mask)
    return result


# --- registration -----------------------------------------------------------


def register() -> None:
    bpy.utils.register_class(MN_RenderLayerSettings)
    bpy.types.ViewLayer.mn_layer = bpy.props.PointerProperty(  # ty: ignore[unresolved-attribute]
        type=MN_RenderLayerSettings
    )
    bpy.types.Scene.mn_layers_active = IntProperty(default=0)  # ty: ignore[unresolved-attribute]
    bpy.types.Scene.mn_layers_background = FloatVectorProperty(  # ty: ignore[unresolved-attribute]
        name="Background",
        description="Colour the layers are composited over. Alpha 0 keeps the render transparent",
        subtype="COLOR",
        size=4,
        min=0.0,
        max=1.0,
        default=(1.0, 1.0, 1.0, 0.0),
    )


def unregister() -> None:
    del bpy.types.Scene.mn_layers_active  # ty: ignore[unresolved-attribute]
    del bpy.types.Scene.mn_layers_background  # ty: ignore[unresolved-attribute]
    del bpy.types.ViewLayer.mn_layer  # ty: ignore[unresolved-attribute]
    bpy.utils.unregister_class(MN_RenderLayerSettings)


# --- proxy -------------------------------------------------------------------


def _entity_objects(scene: bpy.types.Scene) -> list[bpy.types.Object]:
    """Objects in the scene with a Molecular Nodes modifier."""
    return [
        obj
        for obj in scene.objects
        if obj.name != PROXY_NAME and obj.modifiers.get("Molecular Nodes") is not None
    ]


def _sync_proxy(scene: bpy.types.Scene, enabled: bool) -> None:
    """Create, update or remove the ``_ViewLayer`` object that renders every
    entity through ``Filter Render Layer``."""
    from nodebpy import TreeBuilder
    from nodebpy import geometry as g
    from ..nodes._assets import link_asset_groups
    from ..nodes.geometry import FilterRenderLayer

    proxy = bpy.data.objects.get(PROXY_NAME)
    if not enabled:
        for obj in bpy.data.objects:
            if obj.get(SOURCE_PROP):
                obj.hide_render = False
                del obj[SOURCE_PROP]
        if proxy is not None:
            bpy.data.objects.remove(proxy)
        return

    if proxy is None:
        proxy = bpy.data.objects.new(PROXY_NAME, bpy.data.meshes.new(PROXY_NAME))
        proxy.hide_viewport = True
        proxy.hide_select = True
    if proxy.name not in scene.collection.objects:
        scene.collection.objects.link(proxy)
    modifier = proxy.modifiers.get("Filter") or proxy.modifiers.new("Filter", "NODES")
    tree = cast(
        bpy.types.GeometryNodeTree,
        bpy.data.node_groups.get(PROXY_TREE)
        or bpy.data.node_groups.new(PROXY_TREE, "GeometryNodeTree"),
    )
    cast(bpy.types.NodesModifier, modifier).node_group = tree

    sources = _entity_objects(scene)
    link_asset_groups([FilterRenderLayer])
    tree.nodes.clear()
    if tree.interface is not None:
        tree.interface.clear()
    with TreeBuilder(tree) as builder:
        output = builder.outputs.geometry("Geometry")
        joined = g.JoinGeometry(
            [
                g.ObjectInfo(object=obj, transform_space="RELATIVE").o.geometry
                for obj in sources
            ]
        )
        joined >> FilterRenderLayer() >> output
    for obj in sources:
        obj.hide_render = True
        obj[SOURCE_PROP] = True
    _ensure_drivers()


# --- style assignment (used by the UI) --------------------------------------

SET_NODE = "Set Render Layer"


def _is_set_node(node: bpy.types.Node | None) -> bool:
    return (
        isinstance(node, bpy.types.GeometryNodeGroup)
        and node.node_tree is not None
        and node.node_tree.name.startswith(SET_NODE)
    )


def style_layer_node(style: bpy.types.Node) -> bpy.types.GeometryNodeGroup | None:
    """The ``Set Render Layer`` node directly after a style node, if any."""
    for link in style.outputs[0].links or ():
        if _is_set_node(link.to_node):
            return cast(bpy.types.GeometryNodeGroup, link.to_node)
    return None


def set_style_layer(style: bpy.types.Node, index: int) -> None:
    """Put a whole style on a render layer, or on every layer with
    ``index=-1``, inserting or removing the ``Set Render Layer`` node after
    it."""
    from ..nodes._assets import link_asset_groups
    from ..nodes.geometry import SetRenderLayer

    tree = cast(bpy.types.GeometryNodeTree, style.id_data)
    node = style_layer_node(style)
    if index < 0:
        if node is not None:
            targets = [
                socket
                for link in node.outputs[0].links or ()
                if (socket := link.to_socket) is not None
            ]
            tree.nodes.remove(node)
            for socket in targets:
                tree.links.new(style.outputs[0], socket)
        return
    if node is None:
        links = list(style.outputs[0].links or ())
        targets = [socket for link in links if (socket := link.to_socket) is not None]
        # Join Geometry is multi-input: new links add rather than replace
        for link in links:
            tree.links.remove(link)
        link_asset_groups([SetRenderLayer])
        node = cast(bpy.types.GeometryNodeGroup, tree.nodes.new("GeometryNodeGroup"))
        node.node_tree = cast(
            bpy.types.GeometryNodeTree, bpy.data.node_groups[SetRenderLayer._asset_name]
        )
        node.name = node.label = SET_NODE
        node.location = (style.location.x + 220, style.location.y)
        tree.links.new(style.outputs[0], node.inputs[0])
        for socket in targets:
            tree.links.new(node.outputs[0], socket)
    cast(bpy.types.NodeSocketInt, node.inputs["Layer"]).default_value = index
