import os
import subprocess
import sys
import bpy

# Run in a subprocess with its own user scripts folder: uninstalling the template
# from the shared one races with other pytest-xdist workers loading it, which then
# render Blender's factory startup scene instead
_TEMPLATE_SCRIPT = """
import bpy
import pytest
import molecularnodes as mn

mn.assets.template.install()
assert mn.assets.template.find_startup("Molecular Nodes") is not None
bpy.ops.wm.read_homefile(app_template="Molecular Nodes")
assert not bpy.data.objects.get("Cube")

mn.assets.template.uninstall()
assert mn.assets.template.find_startup("Molecular Nodes") is None
# by name alone Blender only warns and falls back to the factory startup file
bpy.ops.wm.read_homefile(app_template="Molecular Nodes")
assert bpy.data.objects.get("Cube")
# a Canvas installs the template itself, but a missing one fails to load loudly
canvas = mn.Canvas()
mn.assets.template.uninstall()
with pytest.raises(ValueError, match="Molecular Nodes"):
    canvas.load_preset("Molecular Nodes")
"""


def test_template(tmp_path):
    env = os.environ.copy()
    env["BLENDER_USER_SCRIPTS"] = str(tmp_path)
    if bpy.app.binary_path:  # inside a full Blender
        cmd = [
            bpy.app.binary_path,
            "-b",
            "--factory-startup",
            "--python-exit-code",
            "1",
            "--python-expr",
        ]
    else:
        cmd = [sys.executable, "-c"]
    result = subprocess.run(
        [*cmd, _TEMPLATE_SCRIPT], capture_output=True, text=True, env=env
    )
    assert result.returncode == 0, result.stdout + result.stderr
