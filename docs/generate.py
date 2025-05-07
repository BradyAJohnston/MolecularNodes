import pathlib
import bpy
import nodepad
import molecularnodes as mn

DOCS_FOLDER = pathlib.Path(__file__).resolve().parent

# load the data file which contains all of the nodes to build docs for
bpy.ops.wm.open_mainfile(filepath=str(mn.assets.MN_DATA_FILE))


header = """---
toc: true
toc-depth: 2
fig-align: center
---
"""

# generate all of the node related docs
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


# write the data table page
with open(DOCS_FOLDER / "data_table.qmd", "w") as file:
    file.write(header)
    file.write("# Data Tables\n\n")
    file.write(
        "The different lookup tables that are used to conver strings to integers in Molecular Nodes.\n\n"
        "Code for this can be found on the [GitHub Page](https://github.com/BradyAJohnston/MolecularNodes/blob/main/molecularnodes/assets/data.py)\n\n"
    )
    file.write("### Residue Names\n\n| Name | Integer |\n|----------:|:------------|\n")
    for name, res in mn.assets.data.residues.items():
        file.write(f"| {name} | `{res['res_name_num']}::Int` |\n")
    file.write("\n")
    file.write("\n")

    file.write("### Atom Names\n\n| Name | Integer |\n|----------:|:------------|\n")
    for name, value in mn.assets.data.atom_names.items():
        file.write(f"| {name} | `{value}::Int` |\n")
    file.write("\n")
