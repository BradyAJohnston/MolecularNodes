import molecularnodes as mn
import bpy
import shutil


file_to_zip = mn.pkg.ADDON_DIR / "assets/template/MolecularNodes"
zip_file_name = mn.pkg.ADDON_DIR / "assets/template/MolecularNodes"

shutil.make_archive(zip_file_name, 'zip', file_to_zip)

def test_template():
    mn.util.utils.template_install()
    bpy.ops.wm.read_homefile(app_template = "MolecularNodes")
    assert not bpy.data.objects.get('Cube')