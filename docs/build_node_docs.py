import bpy
import sys

sys.path.insert(0, "./builder")

from documenter import TreeDocumenter

try:
    from bl_ext.user_default import molecularnodes as mn
except ImportError:
    try:
        from bl_ext.vscode_development import molecularnodes as mn
    except ImportError:
        import molecularnodes as mn

import os
import sys
import pathlib


folder = pathlib.Path(__file__).resolve().parent
file_output_qmd = os.path.join(folder, "nodes/index.qmd")
bpy.ops.wm.open_mainfile(filepath=mn.blender.nodes.MN_DATA_FILE)


header = """---
toc: true
toc-depth: 3
fig-align: center
---
"""

for submenu in mn.ui.node_menu.menu_items.menus:
    with open(os.path.join(folder, f"nodes/{submenu.name}.qmd"), "w") as file:
        file.write(header)
        file.write(f"# {submenu.title}\n\n")
        for item in submenu.items:
            if item.is_break:
                continue
            if item.backup is not None:
                name = item.backup
            else:
                name = item.name

            try:
                file.write(TreeDocumenter(bpy.data.node_groups[name]).printable())
                file.write("\n\n")
            except Exception as e:
                print(e)
                pass
