import molecularnodes as mn
import bpy


def test_template():
    mn.utils.template_install()
    bpy.ops.wm.read_homefile(app_template="Molecular Nodes")
    assert not bpy.data.objects.get('Cube')
