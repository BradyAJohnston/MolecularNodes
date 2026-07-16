import molecularnodes as mn
from molecularnodes.nodes.geometry import StyleCartoon, StyleRibbon, StyleSpheres


def join_nodes(tree):
    return [n for n in tree.tree.nodes if n.bl_idname == "GeometryNodeJoinGeometry"]


def branches(tree):
    "The nodes feeding the tree's final join"
    return [link.from_node for link in join_nodes(tree)[0].inputs[0].links]


def test_atoms_and_join_created_on_demand():
    mol = mn.Molecule.fetch("4ozs")
    with mol.tree as tree:
        tree.atoms >> StyleCartoon() >> tree.join

    assert [i.name for i in mol.tree.tree.interface.items_tree] == ["Geometry", "Atoms"]
    assert len(join_nodes(mol.tree)) == 1
    # the default input -> output passthrough is dropped rather than joined in, so
    # unstyled geometry isn't rendered alongside the style
    assert len(branches(mol.tree)) == 1


def test_atoms_and_join_are_reused_across_branches():
    mol = mn.Molecule.fetch("4ozs")
    with mol.tree as tree:
        tree.atoms >> StyleCartoon() >> tree.join
    with mol.tree as tree:
        tree.atoms >> StyleSpheres() >> tree.join
        tree.atoms >> StyleRibbon() >> tree.join

    # every branch lands on the same join, from the same input, without duplicating
    # the tree interface
    assert len(join_nodes(mol.tree)) == 1
    assert len(branches(mol.tree)) == 3
    assert [i.name for i in mol.tree.tree.interface.items_tree] == ["Geometry", "Atoms"]
    assert len(mol.named_attribute("position", evaluate=True)) > 0


def test_join_keeps_existing_work_feeding_the_output():
    mol = mn.Molecule.fetch("1BNA")
    # a style wired straight to the output, with no join in between
    with mol.tree as tree:
        tree.atoms >> StyleCartoon() >> tree.geometry
    assert join_nodes(mol.tree) == []

    # asking for the join now has to insert one, without orphaning the existing style
    with mol.tree as tree:
        tree.atoms >> StyleSpheres() >> tree.join

    assert len(join_nodes(mol.tree)) == 1
    assert len(branches(mol.tree)) == 2


def test_reset_returns_the_same_sockets():
    mol = mn.Molecule.fetch("4ozs")
    with mol.tree.reset() as (atoms, join):
        atoms >> StyleCartoon() >> join

    assert mol.tree.atoms.socket == atoms.socket
    assert mol.tree.join.socket == join.socket
    assert len(join_nodes(mol.tree)) == 1


def test_reset_with_custom_input_name():
    mol = mn.Molecule.fetch("4ozs")
    with mol.tree.reset(input="Volume") as (volume, join):
        volume >> StyleCartoon() >> join

    assert [i.name for i in mol.tree.tree.interface.items_tree] == [
        "Geometry",
        "Volume",
    ]
    # atoms finds the geometry input whatever it is named, rather than adding another
    assert mol.tree.atoms.socket == volume.socket
    assert [i.name for i in mol.tree.tree.interface.items_tree] == [
        "Geometry",
        "Volume",
    ]
