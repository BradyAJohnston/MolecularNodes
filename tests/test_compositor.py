"""Compositor node assets and the Canvas plumbing that wires them up."""

import bpy
import numpy as np
import pytest
import molecularnodes as mn
from molecularnodes.nodes import compositor as mc
from .constants import data_dir

GROUPS = [
    mc.ConeShadow,
    mc.CompositeDepthFog,
    mc.DepthOutline,
    mc.IDOutline,
    mc.Illustrate,
    mc.Outline,
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


def test_illustrate_enables_passes_and_inserts_group(canvas):
    node = canvas.compositor.illustrate(shading="ao")
    assert {"z", "ambient_occlusion"} <= set(canvas.passes)
    assert "shadow" not in canvas.passes
    assert "normal" not in canvas.passes

    tree = canvas.compositor.tree
    (render_layers,) = _nodes(tree, "CompositorNodeRLayers")
    group = node.node
    assert group.node_tree.name == "Illustrate"
    for name in ["Image", "Alpha", "Depth"]:
        (link,) = group.inputs[name].links
        assert link.from_node == render_layers
    (link,) = group.inputs["Shadow"].links
    assert link.from_socket.name == "Ambient Occlusion"
    assert not group.inputs["Normal"].links
    assert node.i.outline.default_value is True
    assert node.i.shading.default_value is True
    assert node.i.outline_normals.default_value is False

    # the annotation overlay still feeds the output, now reading the composite
    (alpha_over,) = _nodes(tree, "CompositorNodeAlphaOver")
    assert alpha_over.inputs["Background"].links[0].from_node == group
    assert alpha_over.outputs["Image"].links[0].to_node == tree.nodes["Group Output"]
    # the render image itself now only feeds the group
    assert [link.to_node for link in render_layers.outputs["Image"].links] == [group]


def test_illustrate_options(canvas):
    canvas.engine = mn.scene.EEVEE()
    node = canvas.compositor.illustrate(
        shading="shadow", outline=False, outline_normals=True, outline_size=4
    )
    assert {"shadow", "normal"} <= set(canvas.passes)
    assert "ambient_occlusion" not in canvas.passes
    assert node.i.outline.default_value is False
    assert node.i.outline_normals.default_value is True
    assert node.i.outline_size.default_value == 4
    (link,) = node.node.inputs["Shadow"].links
    assert link.from_socket.name == "Shadow"
    (link,) = node.node.inputs["Normal"].links
    assert link.from_socket.name == "Normal"

    with pytest.raises(ValueError):
        canvas.compositor.illustrate(shading="mist")


def test_illustrate_cone_shadow(canvas):
    mol = mn.Molecule.fetch("4ozs", cache=data_dir, format="bcif")
    mol.add_style("cartoon", material=mn.material.Flat())
    canvas.look_at(mol, viewpoint="front")
    node = canvas.compositor.illustrate()
    (link,) = node.node.inputs["Shadow"].links
    cone = link.from_node
    assert cone.node_tree.name == "Cone Shadow"
    (depth,) = cone.inputs["Depth"].links
    assert depth.from_socket.name == "Depth"
    # pixel size from the camera framing, not the node default
    assert cone.inputs["Pixel Size"].default_value > 0
    assert cone.inputs["Pixel Size"].default_value != pytest.approx(0.1)

    # a second call replaces the node and its cone shadow rather than stacking
    canvas.compositor.illustrate(shading="ao")
    groups = [
        n.node_tree.name
        for n in _nodes(canvas.compositor.tree, "CompositorNodeGroup")
        if n.node_tree
    ]
    assert groups.count("Illustrate") == 1
    assert "Cone Shadow" not in groups
    (illustrate,) = [
        n
        for n in _nodes(canvas.compositor.tree, "CompositorNodeGroup")
        if n.node_tree and n.node_tree.name == "Illustrate"
    ]
    assert illustrate.inputs["Image"].links[0].from_node.bl_idname == (
        "CompositorNodeRLayers"
    )


def test_illustrate_shadow_needs_eevee(canvas):
    canvas.engine = mn.scene.Cycles(samples=1, device="CPU")
    with pytest.raises(ValueError, match="EEVEE"):
        canvas.compositor.illustrate(shading="shadow")


def test_illustrate_without_shading_after_reset(canvas):
    with canvas.compositor.reset():
        pass
    node = canvas.compositor.illustrate(shading=None)
    assert node.i.shading.default_value is False
    assert not node.node.inputs["Shadow"].links
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

    node = canvas.compositor.illustrate(
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


def test_cone_shadow_sample_count(canvas):
    "32 samples: 8 directions on 4 rings, as the docs state."
    with canvas.compositor.reset():
        mc.ConeShadow()
    tree = bpy.data.node_groups["Cone Shadow"]
    samples = [
        n
        for n in tree.nodes
        if n.bl_idname == "CompositorNodeGroup"
        and n.node_tree
        and n.node_tree.name == "Cone Shadow Sample"
    ]
    assert len(samples) == 32


def _dilate(mask: np.ndarray, radius: int) -> np.ndarray:
    "Binary dilation of a 2D mask by a square of the given radius."
    h, w = mask.shape
    padded = np.pad(mask, radius)
    out = np.zeros_like(mask)
    for dy in range(-radius, radius + 1):
        for dx in range(-radius, radius + 1):
            out |= padded[radius + dy : radius + dy + h, radius + dx : radius + dx + w]
    return out


def test_render_illustrate(tmp_path):
    "Illustrate changes the image and only extends alpha next to the structure."
    canvas = mn.Canvas(resolution=(64, 64), transparent=True)
    canvas.engine = mn.scene.Cycles(samples=4, device="CPU", denoise=False)
    canvas.compositor.device = "CPU"
    mol = mn.Molecule.fetch("4ozs", cache=data_dir, format="bcif")
    mol.add_style("cartoon", material=mn.material.Flat())
    canvas.look_at(mol, viewpoint="front")

    plain = _render_pixels(canvas, tmp_path / "plain.png")
    canvas.compositor.illustrate(shading="cone", outline=True)
    illustrated = _render_pixels(canvas, tmp_path / "illustrate.png")

    assert plain.shape == illustrated.shape == (64, 64, 4)
    assert not np.allclose(plain, illustrated)
    covered = plain[..., 3] > 0
    assert covered.any()
    # the shadow and outlines darken the structure overall
    assert illustrated[covered, :3].mean() < plain[covered, :3].mean()
    # outlines reach at most a few pixels past the structure; beyond that alpha stays 0
    # (8-bit PNG dithering can leave single-step noise, hence the tolerance)
    near_structure = _dilate(covered, 3)
    assert (~near_structure).any()
    assert illustrated[~near_structure, 3].max() <= 1.5 / 255
