import bpy
import numpy as np
import pytest
import molecularnodes as mn


def test_selection_string_changes(fetch):
    mol = fetch("1BNA")
    sel = mol.selections.from_string("resid 1:10")

    with mol.tree.reset() as (atoms, join):
        atoms >> mn.nodes.geometry.SeparateAtoms(selection=sel.node()) >> join

    current = mol.named_attribute("position", evaluate=True)
    sel.string = "resid 1:5"
    # narrowing the selection reduces the number of separated atoms
    assert len(current) > len(mol.named_attribute("position", evaluate=True))


def test_selection_updating_toggle_rebuilds_atomgroup(fetch):
    """Toggling `updating` swaps between an UpdatingAtomGroup and a static one."""
    mol = fetch("1BNA")
    sel = mol.selections.from_string("resid 1:10")
    manager = mol.selections

    ag = manager.atomgroups[sel.name]
    assert manager.ag_is_updating(ag)

    sel.updating = False
    ag_static = manager.atomgroups[sel.name]
    assert ag_static is not ag
    assert not manager.ag_is_updating(ag_static)
    assert sel.message == ""

    sel.updating = True
    ag_updating = manager.atomgroups[sel.name]
    assert ag_updating is not ag_static
    assert manager.ag_is_updating(ag_updating)


def test_selection_periodic_toggle_rebuilds_atomgroup(fetch):
    """Toggling `periodic` recreates the AtomGroup and re-stores the attribute."""
    mol = fetch("1BNA")
    sel = mol.selections.from_string("around 5.0 resid 1")
    manager = mol.selections

    ag = manager.atomgroups[sel.name]
    sel.periodic = False
    ag_after = manager.atomgroups[sel.name]
    assert ag_after is not ag
    assert sel.message == ""

    # the stored boolean attribute matches the rebuilt AtomGroup
    mask = mol.named_attribute(sel.name)
    assert mask.sum() == ag_after.n_atoms


def test_selection_bad_string_sets_message_and_recovers(fetch):
    """An invalid selection string flags an error and a valid one clears it."""
    mol = fetch("1BNA")
    sel = mol.selections.from_string("resid 1:10")
    mask_before = mol.named_attribute(sel.name)

    sel.string = "not a valid selection"
    assert sel.message != ""
    # the attribute keeps its last good value
    assert np.array_equal(mol.named_attribute(sel.name), mask_before)

    sel.string = "resid 1:5"
    assert sel.message == ""
    assert mol.named_attribute(sel.name).sum() < mask_before.sum()


def test_selections_find_by_string(fetch):
    "`find` searches by selection string; `get` searches by selection name."
    mol = fetch("4ozs")
    assert mol.selections.find("protein") is None
    item = mol.selections.from_string("protein")
    assert mol.selections.find("protein").name == item.name
    # `get` keys on the name, not the string
    assert mol.selections.get("protein") is None
    assert mol.selections.get(item.name) is not None
    # the flags are part of the identity
    assert mol.selections.find("protein", updating=False) is None


def test_selections_node_reuses_existing(fetch):
    "Repeated calls must not pile up duplicate selections and mesh attributes."
    mol = fetch("4ozs")
    with mol.tree:
        for _ in range(5):
            mol.selections.node("protein")
    assert len(mol.selections) == 1


def test_selections_node_flags_are_distinct(fetch):
    mol = fetch("4ozs")
    with mol.tree:
        mol.selections.node("protein")
        mol.selections.node("protein", updating=False)
    assert len(mol.selections) == 2


def test_selections_node_reuses_by_name(fetch):
    mol = fetch("4ozs")
    with mol.tree:
        mol.selections.node("protein", name="my_sel")
        mol.selections.node("a different phrase", name="my_sel")
    assert len([i for i in mol.selections.ui_items if i.name == "my_sel"]) == 1


def test_selections_node_from_atomgroup(fetch):
    mol = fetch("4ozs")
    with mol.tree:
        node = mol.selections.node(mol.universe.select_atoms("resid 1:20"))
    assert node.node.inputs["Name"].default_value == mol.selections.ui_items[0].name


def test_selections_node_requires_tree_context(fetch):
    mol = fetch("4ozs")
    with pytest.raises(RuntimeError, match="TreeBuilder context"):
        mol.selections.node("protein")


def test_selections_node_in_add_style_callable(fetch):
    "The motivating case: an MDA selection inside a callable style."
    from molecularnodes.nodes import geometry as g

    mol = fetch("4ozs")
    mol.add_style(lambda: g.StyleSpheres(selection=mol.selections.node("resid 1:20")))
    style = [
        n
        for n in mol.modifier_node_tree.nodes
        if isinstance(n, bpy.types.GeometryNodeGroup)
        and n.node_tree.name == "Style Spheres"
    ][0]
    assert style.inputs["Selection"].links


def test_select_mask(fetch):
    mol = fetch("1BNA")
    ag = mol.universe.select_atoms("chainID A and resid 1:5")
    n_atoms = mol.universe.atoms.n_atoms

    mask = mol.select("chainID A and resid 1:5")
    assert mask.dtype == bool and mask.shape == (n_atoms,)
    assert np.array_equal(np.flatnonzero(mask), ag.ix)
    assert np.array_equal(mol.select(ag), mask)
    assert np.array_equal(mol.select(ag.ix), mask)
    assert np.array_equal(mol.select(mask), mask)
    assert mol.select().all()
    # an existing boolean attribute is read rather than parsed as a phrase
    assert np.array_equal(mol.select("is_nucleic"), mol.named_attribute("is_nucleic"))
    with pytest.raises(TypeError):
        mol.select("Color")


def test_select_indexing(fetch):
    mol = fetch("1BNA")
    mask = mol.select("resid 1:5")
    before = mol.named_attribute("b_factor")

    # indexing an existing attribute writes back only the selected atoms
    mol["b_factor"][mask] = 0.0
    after = mol.named_attribute("b_factor")
    assert np.all(after[mask] == 0.0)
    assert np.array_equal(after[~mask], before[~mask])

    # a new attribute is stored for every atom
    mol["highlight"] = np.where(mask, 1.0, 0.0)
    assert np.array_equal(mol.named_attribute("highlight"), mask)
    mol["is_site"] = mask
    assert mol.object.data.attributes["is_site"].data_type == "BOOLEAN"
    assert np.array_equal(mol.select("is_site"), mask)


def test_set_color(fetch):
    mol = fetch("1BNA")
    mask = mol.select("chainID A and resid 1:5")
    before = mol.named_attribute("Color")

    # RGB keeps the current alpha
    mol.set_color((0.1, 0.7, 0.7), selection="chainID A and resid 1:5")
    after = mol.named_attribute("Color")
    assert np.allclose(after[mask, :3], (0.1, 0.7, 0.7))
    assert np.allclose(after[:, 3], before[:, 3])
    assert np.allclose(after[~mask], before[~mask])

    mol.set_color((1.0, 0.0, 0.0, 0.5), selection=mask)
    assert np.allclose(mol.named_attribute("Color")[mask], (1.0, 0.0, 0.0, 0.5))

    with pytest.raises(ValueError):
        mol.set_color((1.0, 0.0))


def test_set_color_hex(fetch):
    mol = fetch("1BNA")
    mask = mol.select("chainID A")
    before = mol.named_attribute("Color")

    mol.set_color("#00FFCC", selection="chainID A")
    after = mol.named_attribute("Color")
    assert np.allclose(after[mask, :3], mn.color.parse_color("#00FFCC"))
    assert np.allclose(after[:, 3], before[:, 3])

    mol.set_color("#FF000080", selection=mask)
    assert np.allclose(mol.named_attribute("Color")[mask, 3], 128 / 255)
