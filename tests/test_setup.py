import bpy
import molecularnodes as mn


def test_template():
    mn.assets.template.install()
    bpy.ops.wm.read_homefile(app_template="Molecular Nodes")
    assert not bpy.data.objects.get("Cube")

    mn.assets.template.uninstall()
    try:
        bpy.ops.wm.read_homefile(app_template="Molecular Nodes")
        assert False
    except Exception:
        assert True

    bpy.ops.wm.read_homefile(app_template="")
