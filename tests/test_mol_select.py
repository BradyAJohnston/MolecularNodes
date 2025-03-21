import molecularnodes as mn


def test_select_ca():
    mol = mn.Molecule.fetch("8H1B")
    mol.selector.atom_name("CA").chain_id(["A"]).store_named_attribute(
        "custom_selection"
    )
    assert sum(mol.named_attribute("custom_selection")) == 190
    mol.selector.clear_selections().atom_name("CA").store_named_attribute(
        "custom_selection"
    )
    assert sum(mol.named_attribute("custom_selection")) == 384
    (
        mol.selector.clear_selections()
        .peptide_backbone()
        .res_id(range(50))
        .chain_id(["B"])
        .res_name(["LYS", "ARG", "HIS"])
        .store_named_attribute("custom_selection")
    )
    assert sum(mol.named_attribute("custom_selection")) == 51
