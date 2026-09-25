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

    # a second call replaces the node rather than stacking another
    canvas.compositor.illustrative(outline=True, shading="ao")
    composites = [
        n
        for n in _nodes(canvas.compositor.tree, "CompositorNodeGroup")
        if n.node_tree and n.node_tree.name.startswith("Composite Illustrat")
    ]
    assert len(composites) == 1
    assert composites[0].inputs["Image"].links[0].from_node.bl_idname == (
        "CompositorNodeRLayers"
    )


def test_illustrative_shadow_needs_eevee(canvas):
    canvas.engine = mn.scene.Cycles(samples=1, device="CPU")
    with pytest.raises(ValueError, match="EEVEE"):
        canvas.compositor.illustrative(shading="shadow")


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


ILLUSTRATE_GROUPS = [
    mc.CompositeConeShadow,
    mc.CompositeContourOutline,
    mc.CompositeIDOutline,
    mc.CompositeDepthFog,
    mc.CompositeIllustrate,
]


@pytest.mark.parametrize("group", ILLUSTRATE_GROUPS)
def test_illustrate_group_instantiates(canvas, group):
    with canvas.compositor.reset():
        node = group()
    tree = canvas.compositor.tree
    group_nodes = [n for n in _nodes(tree, "CompositorNodeGroup") if n.node_tree]
    assert group._asset_name in [n.node_tree.name for n in group_nodes]
    for name in group._Inputs.__annotations__:
        getattr(node.i, name)
    for name in group._Outputs.__annotations__:
        getattr(node.o, name)


def test_cone_shadow_sample_count(canvas):
    "32 samples: 8 directions on 4 rings, as the docs state."
    with canvas.compositor.reset():
        mc.CompositeConeShadow()
    tree = bpy.data.node_groups["Composite Cone Shadow"]
    samples = [
        n
        for n in tree.nodes
        if n.bl_idname == "CompositorNodeGroup"
        and n.node_tree
        and n.node_tree.name == "Cone Shadow Sample"
    ]
    assert len(samples) == 32


def test_illustrate_wires_passes_and_aovs(canvas):
    mol = mn.Molecule.fetch("4ozs", cache=data_dir, format="bcif")
    mol.add_style("cartoon", material=mn.material.Flat())
    canvas.look_at(mol, viewpoint="front")
    node = canvas.compositor.illustrate(
        shadow=True, fog=True, contour=True, chain_outline=True, residue_outline=True
    )
    assert {"z", "diffuse_color"} <= set(canvas.passes)
    assert canvas.aovs == ["chain_id", "res_id"]
    tree = canvas.compositor.tree
    (render_layers,) = _nodes(tree, "CompositorNodeRLayers")
    group = node.node
    for name, source in [
        ("Image", "Image"),
        ("Alpha", "Alpha"),
        ("Depth", "Depth"),
        ("Diffuse Color", "Diffuse Color"),
        ("Chain ID", "chain_id"),
        ("Residue ID", "res_id"),
    ]:
        (link,) = group.inputs[name].links
        assert link.from_node == render_layers
        assert link.from_socket.name == source
    for name in [
        "shadow",
        "fog",
        "contour_outline",
        "chain_outline",
        "residue_outline",
    ]:
        assert getattr(node.i, name).default_value is True
    assert node.i.base_color.default_value == "Image"
    # pixel size from the camera and the fog range from the molecule's bounds
    assert node.i.pixel_size.default_value > 0
    assert 0 < node.i.near.default_value < node.i.far.default_value

    # every Molecular Nodes material now writes both AOVs, exactly once, with
    # the attribute offset by one so the first chain differs from the background
    flat = next(m for m in bpy.data.materials if m.name.startswith("Flat"))
    for aov in ["chain_id", "res_id"]:
        assert mn.material.has_aov(flat, aov)
        aov_nodes = [
            n
            for n in flat.node_tree.nodes
            if n.bl_idname == "ShaderNodeOutputAOV" and n.aov_name == aov
        ]
        assert len(aov_nodes) == 1
        (link,) = aov_nodes[0].inputs["Value"].links
        assert link.from_node.bl_idname == "ShaderNodeMath"
    # a second call does not duplicate the AOV outputs
    canvas.compositor.illustrate(chain_outline=True)
    chain_nodes = [
        n
        for n in flat.node_tree.nodes
        if n.bl_idname == "ShaderNodeOutputAOV" and n.aov_name == "chain_id"
    ]
    assert len(chain_nodes) == 1
    (alpha_over,) = _nodes(tree, "CompositorNodeAlphaOver")
    assert alpha_over.outputs["Image"].links[0].to_node == tree.nodes["Group Output"]


def test_illustrate_options(canvas):
    node = canvas.compositor.illustrate(
        shadow=False, contour=False, flat=True, radius=20.0, near=1.0, far=2.0
    )
    assert node.i.shadow.default_value is False
    assert node.i.contour_outline.default_value is False
    assert node.i.base_color.default_value == "Diffuse Color"
    assert node.i.radius.default_value == pytest.approx(20.0)
    assert node.i.near.default_value == pytest.approx(1.0)
    assert canvas.aovs == []
    assert not node.node.inputs["Chain ID"].links


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
    canvas.compositor.illustrate(shadow=True, contour=True, chain_outline=True)
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
