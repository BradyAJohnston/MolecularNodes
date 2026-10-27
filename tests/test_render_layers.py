"""Render layers: Set / Filter Render Layer, the `_ViewLayer` object and the
`canvas.layers` bookkeeping."""

import bpy
import pytest
import molecularnodes as mn
from molecularnodes.nodes import geometry as mg
from molecularnodes.scene import layers as ml


@pytest.fixture
def canvas():
    return mn.Canvas()


def _verts(obj: bpy.types.Object, view_layer: bpy.types.ViewLayer) -> int:
    with bpy.context.temp_override(view_layer=view_layer):
        depsgraph = bpy.context.evaluated_depsgraph_get()
        depsgraph.update()
        mesh = obj.evaluated_get(depsgraph).data
        assert isinstance(mesh, bpy.types.Mesh)
        return len(mesh.vertices)


def _proxy_verts(scene: bpy.types.Scene) -> dict[str, int]:
    """Vertices the `_ViewLayer` object renders on each view layer. It is
    hidden in the viewport, so show it while evaluating."""
    proxy = bpy.data.objects[ml.PROXY_NAME]
    proxy.hide_viewport = False
    try:
        return {vl.name: _verts(proxy, vl) for vl in scene.view_layers}
    finally:
        proxy.hide_viewport = True


def _style_verts(fetch, style: str) -> int:
    mol = fetch("1cbs")
    mol.add_style(style)
    scene = bpy.context.scene
    assert scene is not None
    count = _verts(mol.object, scene.view_layers[0])
    mol.object.modifiers.clear()
    return count


def test_whole_styles_split_by_layer(canvas, fetch):
    cartoon = _style_verts(fetch, "cartoon")
    surface = _style_verts(fetch, "surface")
    mol = fetch("1cbs")
    top = canvas.layers.add("Top")
    mol.add_style("cartoon", layer=0)
    mol.add_style("surface", layer=top)

    assert mol.object.hide_render
    counts = _proxy_verts(canvas.scene)
    assert counts == {"ViewLayer": cartoon, "Top": surface}
    # the viewport still shows everything, from the entity object itself
    assert _verts(mol.object, canvas.scene.view_layers[0]) == cartoon + surface


def test_unassigned_geometry_renders_on_every_layer(canvas, fetch):
    cartoon = _style_verts(fetch, "cartoon")
    surface = _style_verts(fetch, "surface")
    mol = fetch("1cbs")
    canvas.layers.add("Top")
    mol.add_style("cartoon")
    mol.add_style("surface", layer=1)

    assert _proxy_verts(canvas.scene) == {
        "ViewLayer": cartoon,
        "Top": cartoon + surface,
    }


def test_set_render_layer_selection(canvas, fetch):
    "A selection puts part of the geometry on a layer, the rest keeps its own."
    mol = fetch("1cbs")
    canvas.layers.add("Top")
    with mol.tree.reset() as (atoms, join):
        selection = mg.ResidueID().o.res_id < 50
        (
            atoms
            >> mg.SetRenderLayer(layer=0)
            >> mg.SetRenderLayer(selection=selection, layer=1)
            >> join
        )
    canvas.layers.sync()

    total = len(mol.object.data.vertices)
    counts = _proxy_verts(canvas.scene)
    assert counts["ViewLayer"] + counts["Top"] == total
    assert 0 < counts["Top"] < total


def test_sync_index_properties_and_driver(canvas, fetch):
    mol = fetch("1cbs")
    canvas.layers.add("Top")
    mol.add_style("cartoon", layer=1)
    assert [vl[ml.INDEX_PROP] for vl in canvas.scene.view_layers] == [0, 1]
    (node,) = ml._current_nodes()
    tree = node.id_data
    assert isinstance(tree, bpy.types.NodeTree) and tree.animation_data is not None
    drivers = tree.animation_data.drivers
    assert len(drivers) == 1
    driver = drivers[0].driver
    assert driver is not None
    assert driver.is_valid
    # a simple expression needs no Python auto-run
    assert driver.is_simple_expression
    # a plain custom property resolves without the add-on registered
    assert driver.variables[0].targets[0].data_path == f'["{ml.INDEX_PROP}"]'


def test_removing_last_layer_restores_objects(canvas, fetch):
    mol = fetch("1cbs")
    canvas.layers.add("Top")
    mol.add_style("cartoon", layer=1)
    assert ml.PROXY_NAME in bpy.data.objects
    canvas.layers.remove("Top")
    assert ml.PROXY_NAME not in bpy.data.objects
    assert not mol.object.hide_render
    assert ml.SOURCE_PROP not in mol.object


def test_compositor_stacks_layers(canvas, fetch):
    fetch("1cbs").add_style("cartoon", layer=0)
    canvas.layers.add("Top", opacity=0.5)
    tree = canvas.compositor.tree
    render_layers = [n for n in tree.nodes if n.bl_idname == "CompositorNodeRLayers"]
    assert sorted(n.layer for n in render_layers) == ["Top", "ViewLayer"]
    assert canvas.scene.render.film_transparent
    (over,) = [
        n
        for n in tree.nodes
        if n.bl_idname == "CompositorNodeAlphaOver"
        and n.inputs["Factor"].default_value == pytest.approx(0.5)
    ]
    assert over is not None


def test_move_swaps_style_layers(canvas, fetch):
    mol = fetch("1cbs")
    canvas.layers.add("Top")
    mol.add_style("cartoon", layer=0)
    mol.add_style("surface", layer=1)
    canvas.layers.swap_style_layers(0, 1)
    values = sorted(
        (n.inputs["Layer"].default_value, n.inputs[0].links[0].from_node.name)
        for n in mol.tree.tree.nodes
        if ml._is_set_node(n)
    )
    assert values == [(0, "Style Surface"), (1, "Style Cartoon")]
