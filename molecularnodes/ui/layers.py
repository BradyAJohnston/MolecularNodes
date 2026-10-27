"""Render layers panel: a list of layers that hides the view layers behind it.

The list reads the scene's view layers directly, filtered to those tagged as
Molecular Nodes layers and ordered by layer index (top layer first, like an
image editor's layer stack). Operators are thin wrappers over
``molecularnodes.scene.layers``.
"""

from typing import cast
import bpy
from bpy.props import EnumProperty, StringProperty
from bpy.types import Operator, UILayout
from ..scene.layers import RenderLayers, set_style_layer, style_layer_node
from .utils import operator_props


def _settings(view_layer: bpy.types.ViewLayer):
    return view_layer.mn_layer  # ty: ignore[unresolved-attribute]


def _active_view_layer(scene: bpy.types.Scene) -> bpy.types.ViewLayer | None:
    index = scene.mn_layers_active  # ty: ignore[unresolved-attribute]
    if 0 <= index < len(scene.view_layers):
        view_layer = scene.view_layers[index]
        if _settings(view_layer).is_layer:
            return view_layer
    return None


class MN_UL_RenderLayers(bpy.types.UIList):
    """The Molecular Nodes render layers of the scene, top layer first."""

    def draw_item(
        self,
        context,
        layout: UILayout,
        data,
        item,
        icon,
        active_data,
        active_property,
        index=0,
        flt_flag=0,
    ):
        settings = _settings(item)
        row = layout.row(align=True)
        row.label(text=str(settings.index))
        row.prop(item, "name", text="", emboss=False)
        sub = row.row(align=True)
        sub.scale_x = 0.8
        sub.prop(settings, "opacity", text="", slider=True)
        row.prop(
            settings,
            "outline",
            text="",
            icon="MOD_OUTLINE" if settings.outline else "MESH_CIRCLE",
            emboss=False,
        )
        row.prop(item, "use", text="", emboss=False)

    def filter_items(self, context, data, property):
        items = getattr(data, property)
        flags = [
            self.bitflag_filter_item if _settings(vl).is_layer else 0 for vl in items
        ]
        # highest layer at the top of the list
        indices = [_settings(vl).index for vl in items]
        ranked = sorted(range(len(items)), key=lambda i: -indices[i])
        order = [0] * len(items)
        for position, i in enumerate(ranked):
            order[i] = position
        return flags, order


class MN_OT_RenderLayerAdd(Operator):
    bl_idname = "mn.render_layer_add"
    bl_label = "Add Render Layer"
    bl_description = "Add a render layer on top of the stack"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        scene = context.scene
        assert scene is not None
        layer = RenderLayers(scene).add()
        scene.mn_layers_active = list(scene.view_layers).index(layer.view_layer)
        return {"FINISHED"}


class MN_OT_RenderLayerRemove(Operator):
    bl_idname = "mn.render_layer_remove"
    bl_label = "Remove Render Layer"
    bl_description = "Remove the selected render layer and its view layer"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        scene = context.scene
        view_layer = _active_view_layer(scene) if scene else None
        return view_layer is not None and view_layer != scene.view_layers[0]

    def execute(self, context):
        scene = context.scene
        assert scene is not None
        view_layer = _active_view_layer(scene)
        assert view_layer is not None
        RenderLayers(scene).remove(view_layer.name)
        scene.mn_layers_active = 0
        return {"FINISHED"}


class MN_OT_RenderLayerMove(Operator):
    bl_idname = "mn.render_layer_move"
    bl_label = "Move Render Layer"
    bl_description = "Move the selected layer up or down the stack"
    bl_options = {"REGISTER", "UNDO"}

    direction: EnumProperty(items=(("UP", "Up", ""), ("DOWN", "Down", "")))  # type: ignore

    def execute(self, context):
        scene = context.scene
        assert scene is not None
        view_layer = _active_view_layer(scene)
        if view_layer is None:
            return {"CANCELLED"}
        layers = RenderLayers(scene)
        index = _settings(view_layer).index
        target = index + (1 if self.direction == "UP" else -1)
        try:
            other = layers[target]
        except KeyError:
            return {"CANCELLED"}
        # swap the indices; styles follow the index, so they swap layers too
        with layers._paused():
            _settings(other.view_layer).index = index
            _settings(view_layer).index = target
        layers.swap_style_layers(index, target)
        layers.sync()
        return {"FINISHED"}


class MN_OT_RenderLayersComposite(Operator):
    bl_idname = "mn.render_layers_composite"
    bl_label = "Sync Render Layers"
    bl_description = "Update the view layers, the render object and the compositor from the layer list"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        assert context.scene is not None
        RenderLayers(context.scene).sync()
        return {"FINISHED"}


def _layer_items(self, context):
    items = [("-1", "All Layers", "Draw the style on every layer", "RENDERLAYERS", 0)]
    if context is None or context.scene is None:
        return items
    for i, layer in enumerate(RenderLayers(context.scene)):
        items.append((str(layer.index), layer.name, "", "RENDER_RESULT", i + 1))
    return items


class MN_OT_StyleSetLayer(Operator):
    bl_idname = "mn.style_set_layer"
    bl_label = "Set Render Layer"
    bl_description = "Choose the render layer this style is drawn on"
    bl_options = {"REGISTER", "UNDO"}

    name_tree: StringProperty()  # type: ignore
    name_node: StringProperty()  # type: ignore
    layer: EnumProperty(items=_layer_items, name="Layer")  # type: ignore

    def execute(self, context):
        tree = bpy.data.node_groups.get(self.name_tree)
        if tree is None or self.name_node not in tree.nodes:
            return {"CANCELLED"}
        set_style_layer(tree.nodes[self.name_node], int(self.layer))
        assert context.scene is not None
        RenderLayers(context.scene).sync()
        return {"FINISHED"}


def panel_layers(layout: UILayout, context: bpy.types.Context) -> None:
    scene = context.scene
    assert scene is not None
    row = layout.row()
    row.template_list(
        "MN_UL_RenderLayers",
        "mn_render_layers",
        scene,
        "view_layers",
        scene,
        "mn_layers_active",
        rows=3,
    )
    col = row.column(align=True)
    col.operator("mn.render_layer_add", icon="ADD", text="")
    col.operator("mn.render_layer_remove", icon="REMOVE", text="")
    col.separator()
    operator_props(
        col.operator("mn.render_layer_move", icon="TRIA_UP", text="")
    ).direction = "UP"
    operator_props(
        col.operator("mn.render_layer_move", icon="TRIA_DOWN", text="")
    ).direction = "DOWN"

    layout.use_property_split = True
    layout.use_property_decorate = False
    view_layer = _active_view_layer(scene)
    if view_layer is not None:
        settings = _settings(view_layer)
        col = layout.column()
        col.prop(settings, "opacity")
        col.prop(settings, "blend")
        header, body = layout.panel("mn_render_layer_outline")
        header.use_property_split = False
        header.prop(settings, "outline")
        if body is not None:
            body.active = settings.outline
            body.prop(settings, "outline_color", text="Color")
            body.prop(settings, "outline_size", text="Size")
            body.prop(settings, "outline_threshold", text="Threshold")
            body.prop(settings, "outline_normals")
    header, body = layout.panel("mn_render_layers_options", default_closed=True)
    header.label(text="Options")
    if body is not None:
        body.prop(scene, "mn_layers_background")
        body.prop(scene, "mn_layers_auto_composite")
        body.operator("mn.render_layers_composite")


def draw_style_layer(layout: UILayout, context, style: bpy.types.Node) -> None:
    """The render layer picker for a style, drawn in the Styles panel."""
    assert context.scene is not None
    layers = RenderLayers(context.scene)
    if len(layers) < 2:
        return
    node = style_layer_node(style)
    current = "All Layers"
    if node is not None:
        index = cast(bpy.types.NodeSocketInt, node.inputs["Layer"]).default_value
        try:
            current = layers[index].name
        except KeyError:
            current = f"Layer {index} (missing)"
    row = layout.row(align=True)
    row.label(text="Render Layer:")
    op = operator_props(
        row.operator_menu_enum("mn.style_set_layer", "layer", text=current)
    )
    assert style.id_data is not None
    op.name_tree = style.id_data.name
    op.name_node = style.name


class MN_PT_RenderLayers(bpy.types.Panel):
    bl_idname = "MN_PT_render_layers"
    bl_label = "Layers"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Molecular Nodes"
    bl_parent_id = "MN_PT_render"

    def draw(self, context):
        assert self.layout is not None
        panel_layers(self.layout, context)


CLASSES = [
    MN_UL_RenderLayers,
    MN_OT_RenderLayerAdd,
    MN_OT_RenderLayerRemove,
    MN_OT_RenderLayerMove,
    MN_OT_RenderLayersComposite,
    MN_OT_StyleSetLayer,
    MN_PT_RenderLayers,
]
