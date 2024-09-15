import bpy
import sys

sys.path.insert(0, "./")

from noodlenotes.documenter import documenter

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
toc-depth: 2
fig-align: center
---
"""

for submenu in mn.ui.node_menu.menu_items.submenus:
    with open(os.path.join(folder, f"nodes/{submenu.name}.qmd"), "w") as file:
        file.write(header)
        file.write(f"# {submenu.title}\n\n")
        if submenu.description:
            file.write(submenu.description)
            file.write("\n\n")
        for item in submenu.items:
            if item.is_break:
                continue
            if item.backup is not None:
                name = item.backup
            else:
                name = item.name
            doc = documenter(item)

            file.write(doc.as_markdown())
            file.write("\n\n")
