import molecularnodes as mn

mol = mn.Molecule.fetch("4ozs")
mol.add_style(
    mn.StyleCartoon(peptide_loop_radius=0.4), material=mn.material.AmbientOcclusion()
)
