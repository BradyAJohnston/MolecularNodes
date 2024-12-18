import bpy
import pathlib

import molecularnodes as mn
import nodepad

DOCS_FOLDER = pathlib.Path(__file__).resolve().parent

# load the data file which contains all of the nodes to build docs for
bpy.ops.wm.open_mainfile(filepath=mn.utils.MN_DATA_FILE)


header = """---
toc: true
toc-depth: 2
fig-align: center
---
"""

for submenu in mn.ui.node_menu.menu_items.submenus:
    with open(DOCS_FOLDER / f"nodes/{submenu.name}.qmd", "w") as file:
        file.write(header)
        file.write(f"# {submenu.title}\n\n")
        if submenu.description:
            file.write(submenu.description)
            file.write("\n\n")
        for menu_item in submenu.items:
            if menu_item.is_break:
                continue
            if menu_item.backup is not None:
                name = menu_item.backup
            else:
                name = menu_item.name
            doc = nodepad.Documenter(menu_item.tree)
            try:
                doc.lookup_info(menu_item.to_dict())
            except AttributeError as e:
                print(e)

            if menu_item.description != "":
                doc.description += "\n\n" + menu_item.description

            file.write(doc.as_markdown())
            file.write("\n\n")
