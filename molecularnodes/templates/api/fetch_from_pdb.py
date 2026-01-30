import addon_utils

addon_id = "bl_ext.blender_org.molecularnodes"
is_enabled, is_loaded = addon_utils.check(addon_id)
if is_enabled and is_loaded:
    import bl_ext.blender_org.molecularnodes as mn  # type: ignore
else:
    import molecularnodes as mn

mol = mn.Molecule.fetch("4ozs")
mol.add_style(
    mn.StyleCartoon(peptide_loop_radius=0.4), material=mn.material.AmbientOcclusion()
)
