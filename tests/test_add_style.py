import bpy
import pytest
import molecularnodes as mn
from .constants import data_dir


def _sphere(mol) -> str | None:
    "The `Sphere Geometry` value on whichever style node has one."
    for node in mol.tree.nodes:
        for socket in node.inputs:
            if socket.name == "Sphere Geometry":
                return socket.default_value
    return None


def test_add_style_with_selection():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    mol.store_named_attribute(mol.named_attribute("is_side_chain"), "show_side_chains")
    mol.add_style("ball_and_stick", selection="show_side_chains")

    node_style = [
        node
        for node in mol.modifier_node_tree.nodes
        if (
            isinstance(node, bpy.types.GeometryNodeGroup)
            and node.node_tree.name == "Style Ball and Stick"
        )
    ][0]

    assert (
        node_style.inputs["Selection"].links[0].from_node.inputs["Name"].default_value
        == "show_side_chains"
    )

    with pytest.warns(UserWarning):
        mol.add_style("cartoon", selection="non_existing_selection")


def _nodes_using_tree(node_group, tree_name):
    return [
        node
        for node in node_group.nodes
        if isinstance(node, bpy.types.GeometryNodeGroup)
        and node.node_tree.name == tree_name
    ]


def test_add_style_color_scheme_and_name():
    mol = mn.Molecule.load(data_dir / "1cd3.cif").add_style(
        "cartoon", color="common", name="My Style"
    )
    ng = mol.modifier_node_tree

    # a Set Color node is inserted upstream of the style, driven by the
    # element-colors chain with per-chain random carbons
    (set_color,) = _nodes_using_tree(ng, "Set Color")
    color_source = set_color.inputs["Color"].links[0].from_node
    assert color_source.node_tree.name == "Color Element"

    (style,) = _nodes_using_tree(ng, "Style Cartoon")
    assert style.label == "My Style"
    assert style.inputs["Atoms"].links[0].from_node == set_color


def test_add_style_uniform_color():
    color = (0.1, 0.2, 0.3, 1.0)
    mol = mn.Molecule.load(data_dir / "1cd3.cif").add_style("spheres", color=color)
    (set_color,) = _nodes_using_tree(mol.modifier_node_tree, "Set Color")
    assert not set_color.inputs["Color"].links
    assert tuple(set_color.inputs["Color"].default_value) == pytest.approx(color)


def test_add_style_assembly():
    mol = mn.Molecule.load(data_dir / "1cd3.cif").add_style("cartoon", assembly=True)
    ng = mol.modifier_node_tree

    (instance,) = _nodes_using_tree(ng, "Assembly Instance")
    (style,) = _nodes_using_tree(ng, "Style Cartoon")
    assert instance.inputs["Geometry"].links[0].from_node == style
    assert instance.inputs["Data Object"].default_value is not None


def test_add_style_operator_color_scheme():
    """mn.add_style applies its color scheme and label through Molecule.add_style."""
    mol = mn.Molecule.load(data_dir / "1cd3.cif")
    res = bpy.ops.mn.add_style(
        uuid=mol.uuid, style="spheres", color_scheme="common", name="Op Style"
    )
    assert res == {"FINISHED"}

    ng = mol.modifier_node_tree
    assert _nodes_using_tree(ng, "Set Color")
    (style,) = _nodes_using_tree(ng, "Style Spheres")
    assert style.label == "Op Style"


def test_styles_panel_style_node_selection():
    """The Styles panel identifies the active style node from the list index."""
    from molecularnodes.ui.panel import is_style_node

    mol = mn.Molecule.load(data_dir / "1cd3.cif").add_style("cartoon")
    node_group = mol.modifier_node_tree
    style_names = [n.name for n in node_group.nodes if is_style_node(n)]
    assert style_names, "adding a style should create a style node"

    # the panel treats styles_active_index as an index into node_group.nodes and
    # only shows details when it points at a style node
    index = node_group.nodes.find(style_names[0])
    mol.props.styles_active_index = index

    active = node_group.nodes[mol.props.styles_active_index]
    assert is_style_node(active)
    assert active.name == style_names[0]


def test_styles_panel_node_group_lookup_without_session():
    """Style nodes are found from the object alone, without a session link."""
    from molecularnodes.ui.panel import get_entity_node_group, is_style_node

    session = mn.session.get_session()
    mol = mn.Molecule.load(data_dir / "1cd3.cif").add_style("cartoon")
    obj = mol.object

    # drop the session link (as after re-opening a .blend without the pickle)
    session.remove_entity(obj.uuid)
    assert session.get(obj.uuid) is None

    # the node group and its style nodes are still reachable from the object
    node_group = get_entity_node_group(obj)
    assert node_group is not None
    style_names = [n.name for n in node_group.nodes if is_style_node(n)]
    assert style_names


def test_swap_style_operator():
    """mn.swap_style swaps the style node's tree in place, keeping connections."""
    from molecularnodes.ui.panel import is_style_node

    mol = mn.Molecule.load(data_dir / "1cd3.cif").add_style("cartoon")
    ng = mol.modifier_node_tree
    style = [n for n in ng.nodes if is_style_node(n)][0]
    assert style.node_tree.name == "Style Cartoon"

    res = bpy.ops.mn.swap_style(
        name_tree=ng.name, name_node=style.name, style="surface"
    )
    assert res == {"FINISHED"}

    # still exactly one style node, now referencing the Surface style tree
    style_nodes = [n for n in ng.nodes if is_style_node(n)]
    assert len(style_nodes) == 1
    assert style_nodes[0].node_tree.name == "Style Surface"


def test_add_selection_to_style_operator():
    """mn.add_selection_to_style wires a named-attribute selection into the style node."""
    from molecularnodes.ui.panel import is_style_node

    mol = mn.Molecule.load(data_dir / "1cd3.cif").add_style("cartoon")
    ng = mol.modifier_node_tree
    style = [n for n in ng.nodes if is_style_node(n)][0]
    assert not style.inputs["Selection"].links

    bpy.context.view_layer.objects.active = mol.object
    res = bpy.ops.mn.add_selection_to_style(node_tree=ng.name, node_name=style.name)
    assert res == {"FINISHED"}

    # the style's Selection input is now driven by a named attribute node whose
    # name matches a selection on the molecule, as the panel expects
    links = style.inputs["Selection"].links
    assert links
    attr_node = links[0].from_node
    assert isinstance(attr_node, bpy.types.GeometryNodeInputNamedAttribute)
    selection = mol.selections.get(attr_node.inputs["Name"].default_value)
    assert selection is not None
    assert selection.string == "all"


def test_add_selection_to_style_operator_missing_targets():
    """mn.add_selection_to_style errors cleanly when the tree or node is gone."""
    from molecularnodes.ui.panel import is_style_node

    mol = mn.Molecule.load(data_dir / "1cd3.cif").add_style("cartoon")
    ng = mol.modifier_node_tree
    style = [n for n in ng.nodes if is_style_node(n)][0]
    bpy.context.view_layer.objects.active = mol.object

    with pytest.raises(RuntimeError, match="not found"):
        bpy.ops.mn.add_selection_to_style(
            node_tree="does not exist", node_name=style.name
        )

    with pytest.raises(RuntimeError, match="not found"):
        bpy.ops.mn.add_selection_to_style(node_tree=ng.name, node_name="does not exist")

    # the style node is untouched by the failed attempts
    assert not style.inputs["Selection"].links


def test_add_style_invalid_style_raises():
    "A style string that isn't dispatchable should raise, not KeyError deeper down."
    mol = mn.Molecule.fetch("4ozs")
    # 'vdw' used to pass validation against `styles_mapping` and then KeyError
    for style in ("vdw", "atoms", "sphere", "bogus"):
        with pytest.raises(ValueError, match="Invalid style"):
            mol.add_style(style)


@pytest.mark.parametrize("color", ["chain", "lipophobicity", "b_factor", "nonsense"])
def test_add_style_unknown_color_warns(color):
    """Unknown colors, and non-color attributes, must warn rather than render black.

    `lipophobicity` and `b_factor` exist but are FLOAT attributes - reading them as a
    color silently produced black geometry.
    """
    mol = mn.Molecule.fetch("4ozs")
    with pytest.warns(UserWarning, match="is neither a known color keyword"):
        mol.add_style("cartoon", color=color)
    # no Set Color node should have been added
    assert not _nodes_using_tree(mol.modifier_node_tree, "Set Color")


@pytest.mark.parametrize("color", ["common", "default", "plddt", "Color"])
def test_add_style_valid_color_adds_node(color):
    mol = mn.Molecule.fetch("4ozs")
    mol.add_style("cartoon", color=color)
    assert _nodes_using_tree(mol.modifier_node_tree, "Set Color")


def test_add_style_callable_color():
    "A callable is evaluated inside the tree context, reaching any Color* node."
    from molecularnodes.nodes import geometry as mg

    mol = mn.Molecule.fetch("4ozs")
    mol.add_style("cartoon", color=lambda: mg.ColorRainbow())
    assert _nodes_using_tree(mol.modifier_node_tree, "Set Color")
    assert _nodes_using_tree(mol.modifier_node_tree, "Color Rainbow")


def test_add_style_callable_selection():
    from molecularnodes.nodes import geometry as mg

    mol = mn.Molecule.fetch("4ozs")
    mol.add_style("sticks", selection=lambda: mg.IsPeptide() & mg.IsSideChain())
    style = _nodes_using_tree(mol.modifier_node_tree, "Style Sticks")[0]
    # the selection input is driven by a link, not a named attribute lookup
    assert style.inputs["Selection"].links


def test_add_style_callable_style():
    from molecularnodes.nodes import geometry as mg

    mol = mn.Molecule.fetch("4ozs")
    mol.add_style(lambda: mg.StyleCartoon(quality=5))
    style = _nodes_using_tree(mol.modifier_node_tree, "Style Cartoon")[0]
    assert style.inputs["Quality"].default_value == 5


@pytest.mark.parametrize(
    "kwargs",
    [
        {"selection": "protein"},
        {"material": "MN Default"},
        {"quality": 5},
    ],
)
def test_add_style_callable_style_rejects_owned_args(kwargs):
    "Args belonging to the style node must not be silently dropped."
    from molecularnodes.nodes import geometry as mg

    mol = mn.Molecule.fetch("4ozs")
    with pytest.raises(TypeError, match="cannot also be passed"):
        mol.add_style(lambda: mg.StyleCartoon(), **kwargs)


def test_add_style_callable_style_allows_color_and_name():
    "color/assembly/name are separate nodes, so they compose with a callable style."
    from molecularnodes.nodes import geometry as mg

    mol = mn.Molecule.fetch("4ozs")
    mol.add_style(lambda: mg.StyleCartoon(), color=lambda: mg.ColorRainbow(), name="x")
    assert _nodes_using_tree(mol.modifier_node_tree, "Color Rainbow")
    assert _nodes_using_tree(mol.modifier_node_tree, "Style Cartoon")[0].label == "x"


@pytest.mark.parametrize(
    "engine,expected",
    [("EEVEE", "Instance"), ("CYCLES", "Point")],
)
def test_spheres_geometry_follows_the_engine(engine, expected):
    """Point clouds are only ray-traced by Cycles.

    Left on the `"Point"` default under EEVEE, spheres render as octahedra -
    invisible at whole-protein zoom, glaring on a close-up.
    """
    canvas = mn.Canvas(engine=engine)
    mol = mn.Molecule.fetch("4ozs").add_style("spheres")
    assert _sphere(mol) == expected
    assert canvas.scene.render.engine == (
        "BLENDER_EEVEE" if engine == "EEVEE" else "CYCLES"
    )


def test_explicit_sphere_is_not_overridden():
    mn.Canvas(engine="EEVEE")
    mol = mn.Molecule.fetch("4ozs").add_style("spheres", sphere="Point")
    assert _sphere(mol) == "Point"


def test_ball_and_stick_already_instances_and_is_left_alone():
    "It defaults to `Instance` already, so nothing should change under either engine."
    for engine in ("EEVEE", "CYCLES"):
        mn.Canvas(engine=engine)
        mol = mn.Molecule.fetch("4ozs").add_style("ball_and_stick")
        assert _sphere(mol) == "Instance"


def test_callable_style_keeps_its_own_geometry():
    "A callable builds the node itself, so nothing is imposed on it."
    from molecularnodes.nodes.geometry import StyleSpheres

    mn.Canvas(engine="EEVEE")
    mol = mn.Molecule.fetch("4ozs").add_style(lambda: StyleSpheres())
    assert _sphere(mol) == "Point"
