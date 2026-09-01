import bpy
import numpy as np
import pytest
import molecularnodes as mn


def test_selection_string_changes():
    mol = mn.Molecule.fetch("1bna")
    sel = mol.selections.from_string("resid 1:10")

    with mol.tree.reset() as (atoms, join):
        atoms >> mn.nodes.geometry.SeparateAtoms(selection=sel.node()) >> join

    current = mol.named_attribute("position", evaluate=True)
    sel.string = "resid 1:5"
    # narrowing the selection reduces the number of separated atoms
    assert len(current) > len(mol.named_attribute("position", evaluate=True))


def test_selection_updating_toggle_rebuilds_atomgroup():
    """Toggling `updating` swaps between an UpdatingAtomGroup and a static one."""
    mol = mn.Molecule.fetch("1bna")
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


def test_selection_periodic_toggle_rebuilds_atomgroup():
    """Toggling `periodic` recreates the AtomGroup and re-stores the attribute."""
    mol = mn.Molecule.fetch("1bna")
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


def test_selection_bad_string_sets_message_and_recovers():
    """An invalid selection string flags an error and a valid one clears it."""
    mol = mn.Molecule.fetch("1bna")
    sel = mol.selections.from_string("resid 1:10")
    mask_before = mol.named_attribute(sel.name)

    sel.string = "not a valid selection"
    assert sel.message != ""
    # the attribute keeps its last good value
    assert np.array_equal(mol.named_attribute(sel.name), mask_before)

    sel.string = "resid 1:5"
    assert sel.message == ""
    assert mol.named_attribute(sel.name).sum() < mask_before.sum()


def test_selections_find_by_string():
    "`find` searches by selection string; `get` searches by selection name."
    mol = mn.Molecule.fetch("4ozs")
    assert mol.selections.find("protein") is None
    item = mol.selections.from_string("protein")
    assert mol.selections.find("protein").name == item.name
    # `get` keys on the name, not the string
    assert mol.selections.get("protein") is None
    assert mol.selections.get(item.name) is not None
    # the flags are part of the identity
    assert mol.selections.find("protein", updating=False) is None


def test_selections_node_reuses_existing():
    "Repeated calls must not pile up duplicate selections and mesh attributes."
    mol = mn.Molecule.fetch("4ozs")
    with mol.tree:
        for _ in range(5):
            mol.selections.node("protein")
    assert len(mol.selections) == 1


def test_selections_node_flags_are_distinct():
    mol = mn.Molecule.fetch("4ozs")
    with mol.tree:
        mol.selections.node("protein")
        mol.selections.node("protein", updating=False)
    assert len(mol.selections) == 2


def test_selections_node_reuses_by_name():
    mol = mn.Molecule.fetch("4ozs")
    with mol.tree:
        mol.selections.node("protein", name="my_sel")
        mol.selections.node("a different phrase", name="my_sel")
    assert len([i for i in mol.selections.ui_items if i.name == "my_sel"]) == 1


def test_selections_node_from_atomgroup():
    mol = mn.Molecule.fetch("4ozs")
    with mol.tree:
        node = mol.selections.node(mol.universe.select_atoms("resid 1:20"))
    assert node.node.inputs["Name"].default_value == mol.selections.ui_items[0].name


def test_selections_node_requires_tree_context():
    mol = mn.Molecule.fetch("4ozs")
    with pytest.raises(RuntimeError, match="TreeBuilder context"):
        mol.selections.node("protein")


def test_selections_node_in_add_style_callable():
    "The motivating case: an MDA selection inside a callable style."
    from molecularnodes.nodes import geometry as g

    mol = mn.Molecule.fetch("4ozs")
    mol.add_style(lambda: g.StyleSpheres(selection=mol.selections.node("resid 1:20")))
    style = [
        n
        for n in mol.modifier_node_tree.nodes
        if isinstance(n, bpy.types.GeometryNodeGroup)
        and n.node_tree.name == "Style Spheres"
    ][0]
    assert style.inputs["Selection"].links
