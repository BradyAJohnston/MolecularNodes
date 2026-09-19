"""Compositor node assets and the Canvas plumbing that wires them up."""

import bpy
import numpy as np
import pytest
import molecularnodes as mn
from molecularnodes.nodes import compositor as mc
from .constants import data_dir

GROUPS = [
    mc.CompositeOutlineMask,
    mc.CompositeOutline,
    mc.CompositeIllustrative,
    mc.ValueToMask,
    mc.MaskPixels,
]


@pytest.fixture
def canvas():
    return mn.Canvas()


def _nodes(tree, bl_idname):
    return [n for n in tree.nodes if n.bl_idname == bl_idname]


@pytest.mark.parametrize("group", GROUPS)
def test_compositor_group_instantiates(canvas, group):
    "Each asset appends from nodes.blend into the scene compositor without error."
    with canvas.compositor.reset():
        node = group()
    tree = canvas.compositor.tree
    group_nodes = [n for n in _nodes(tree, "CompositorNodeGroup") if n.node_tree]
    assert group._asset_name in [n.node_tree.name for n in group_nodes]
    assert node.node.node_tree.bl_idname == "CompositorNodeTree"
    # every documented input and output resolves to a socket on the node
    for name in group._Inputs.__annotations__:
        getattr(node.i, name)
    for name in group._Outputs.__annotations__:
        getattr(node.o, name)


def test_illustrative_enables_passes_and_inserts_group(canvas):
    node = canvas.compositor.illustrative(outline=True, shading="ao")
    assert {"z", "normal", "diffuse_color", "ambient_occlusion"} <= set(canvas.passes)
    assert "shadow" not in canvas.passes

    tree = canvas.compositor.tree
    (render_layers,) = _nodes(tree, "CompositorNodeRLayers")
    group = node.node
    assert group.node_tree.name == "Composite Illustrative"
    for name in ["Image", "Alpha", "Depth", "Normal", "Diffuse Color"]:
        (link,) = group.inputs[name].links
        assert link.from_node == render_layers
    (link,) = group.inputs["Ambient Occlusion"].links
    assert link.from_socket.name == "Ambient Occlusion"
    assert not group.inputs["Shadow"].links
    assert node.i.outline.default_value is True
    assert node.i.shading.default_value is True
    assert node.i.shading_source.default_value == "Ambient Occlusion"
    assert node.i.base_color.default_value == "Image"

    # the annotation overlay still feeds the output, now reading the composite
    (alpha_over,) = _nodes(tree, "CompositorNodeAlphaOver")
    assert alpha_over.inputs["Background"].links[0].from_node == group
    assert alpha_over.outputs["Image"].links[0].to_node == tree.nodes["Group Output"]
    # the render image itself now only feeds the group
    assert [link.to_node for link in render_layers.outputs["Image"].links] == [group]


def test_illustrative_options(canvas):
    node = canvas.compositor.illustrative(
        outline=False, shading="shadow", flat=True, outline_size=4
    )
    assert "shadow" in canvas.passes
    assert "ambient_occlusion" not in canvas.passes
    assert node.i.shading_source.default_value == "Shadow"
    assert node.i.base_color.default_value == "Diffuse Color"
    assert node.i.outline.default_value is False
    assert node.i.outline_size.default_value == 4
    (link,) = node.node.inputs["Shadow"].links
    assert link.from_socket.name == "Shadow"

    with pytest.raises(ValueError):
        canvas.compositor.illustrative(shading="mist")


def test_illustrative_without_shading_after_reset(canvas):
    with canvas.compositor.reset():
        pass
    node = canvas.compositor.illustrative(shading=None)
    assert node.i.shading.default_value is False
    assert not node.node.inputs["Ambient Occlusion"].links
    tree = canvas.compositor.tree
    assert node.node.outputs[0].links[0].to_node == tree.nodes["Group Output"]


def test_shadow_pass(canvas):
    canvas.passes = ["combined", "shadow"]
    assert canvas.scene.view_layers[0].use_pass_shadow is True
    assert set(canvas.passes) == {"combined", "shadow"}


def test_add_aov(canvas):
    aov = canvas.add_aov("chain_id")
    assert aov.name == "chain_id"
    assert aov.type == "VALUE"
    # adding the same name again returns the existing pass
    canvas.add_aov("chain_id")
    assert canvas.aovs == ["chain_id"]
    canvas.add_aov("tint", type="COLOR")
    assert canvas.aovs == ["chain_id", "tint"]
    assert canvas.scene.view_layers[0].aovs["tint"].type == "COLOR"
    # the pass is an output of the Render Layers node for Value to Mask to read
    with canvas.compositor as tree:
        (render_layers,) = _nodes(tree.tree, "CompositorNodeRLayers")
        assert "chain_id" in [s.name for s in render_layers.outputs]


def test_material_add_aov():
    mat = mn.material.Flat()
    n_before = len(mat.material.node_tree.nodes)
    node = mat.add_aov("chain_id")
    tree = mat.material.node_tree
    assert len(tree.nodes) == n_before + 2
    assert node.node.bl_idname == "ShaderNodeOutputAOV"
    assert node.node.aov_name == "chain_id"
    (link,) = node.node.inputs["Value"].links
    assert link.from_node.bl_idname == "ShaderNodeAttribute"
    assert link.from_node.attribute_name == "chain_id"
    assert link.from_socket.name == "Factor"

    colour = mn.material.add_aov(mat.material, "tint", attribute="Color", type="COLOR")
    (link,) = colour.node.inputs["Color"].links
    assert link.from_node.attribute_name == "Color"
    assert link.from_socket.name == "Color"


def _render_pixels(canvas, path) -> np.ndarray:
    canvas.snapshot(path)
    image = bpy.data.images.load(str(path))
    try:
        return np.array(image.pixels[:]).reshape(image.size[1], image.size[0], 4)
    finally:
        bpy.data.images.remove(image)


def test_render_outline_changes_pixels(tmp_path):
    "A low-res render with the outline on differs from the same render with it off."
    canvas = mn.Canvas(resolution=(64, 64), transparent=True)
    canvas.engine = mn.scene.Cycles(samples=4, device="CPU", denoise=False)
    canvas.compositor.device = "CPU"
    mol = mn.Molecule.fetch("4ozs", cache=data_dir, format="bcif")
    mol.add_style("spheres")
    canvas.look_at(mol, viewpoint="front")

    node = canvas.compositor.illustrative(
        outline=False, shading=None, outline_color=(1.0, 0.0, 0.0, 1.0)
    )
    plain = _render_pixels(canvas, tmp_path / "plain.png")
    node.i.outline.default_value = True
    outlined = _render_pixels(canvas, tmp_path / "outlined.png")

    assert plain.shape == outlined.shape == (64, 64, 4)
    assert plain[..., 3].max() > 0.5  # something was rendered
    assert not np.allclose(plain, outlined)
    # the lines are opaque red pixels that were not there before
    red = (outlined[..., 0] > 0.5) & (outlined[..., 1] < 0.2) & (outlined[..., 3] > 0.5)
    assert red.sum() > 0
    assert red.sum() > ((plain[..., 0] > 0.5) & (plain[..., 1] < 0.2)).sum()
