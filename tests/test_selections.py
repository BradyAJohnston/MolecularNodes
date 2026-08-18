import numpy as np
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
