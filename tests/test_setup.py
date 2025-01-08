import molecularnodes as mn
import bpy


def test_template():
    mn.template.install()
    bpy.ops.wm.read_homefile(app_template="Molecular Nodes")
    assert not bpy.data.objects.get("Cube")

    mn.template.uninstall()
    try:
        bpy.ops.wm.read_homefile(app_template="Molecular Nodes")
        assert False
    except Exception:
        assert True

    bpy.ops.wm.read_homefile(app_template="")
