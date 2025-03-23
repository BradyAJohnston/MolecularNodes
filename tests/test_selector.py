import molecularnodes as mn


def test_select_ca():
    mol = mn.Molecule.fetch("8H1B")
    mol.select.atom_name("CA").chain_id(["A"]).store_selection("custom_selection")
    assert sum(mol.named_attribute("custom_selection")) == 190
    mol.select.reset().atom_name("CA").store_selection("custom_selection")
    assert sum(mol.named_attribute("custom_selection")) == 384
    (
        mol.select.reset()
        .is_peptide_backbone()
        .res_id(range(50))
        .chain_id(["B"])
        .res_name(["LYS", "ARG", "HIS"])
        .store_selection("custom_selection")
    )
    assert sum(mol.named_attribute("custom_selection")) == 18


def test_evalute_on_array():
    mol = mn.Molecule.fetch("8H1B")
    sel = mn.entities.MoleculeSelector()
    sel.atom_name("CA").chain_id(["A"]).res_id(range(50))

    assert sum(sel.evaluate_on_array(mol.array)) == 49


def test_adjust_seletion():
    mol = mn.Molecule.fetch("8H1B")
    sel = mol.select

    (
        sel.res_id(range(200))
        .chain_id(["A", "B"])
        .res_name(["LYS", "ARG", "HIS", "GLY", "PRO"])
        .is_backbone()
    )

    assert sum(sel.evaluate_on_array(mol.array)) == 400
    sel.res_id(range(100, 200))
    assert sum(sel.evaluate_on_array(mol.array)) == 192
