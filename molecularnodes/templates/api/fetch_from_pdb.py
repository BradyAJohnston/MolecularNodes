import bl_ext.blender_org.molecularnodes as mn  # type: ignore

mol = mn.Molecule.fetch("4ozs")
mol.add_style("cartoon", material="MN Ambient Occlusion")
