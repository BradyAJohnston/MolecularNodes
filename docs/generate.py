"""Generate the node reference structure for the quartodoc API docs.

The .blend asset file is the source of truth: every node group marked as an
asset in `molecularnodes/assets/node_data_file.blend` is documented, grouped
by its asset catalog (the same grouping the GUI shows). Each top-level catalog
becomes one quartodoc page listing the generated node classes for that
category (see generate_node_classes.py), written into the marked block of
`docs/_quarto.yml` before `quartodoc build` runs.

Extra information that cannot live on the nodes themselves (long-form prose,
demo videos) lives in `docs/nodes.yml`, keyed by node group name, and is
injected into the rendered class pages by the custom renderer in
`docs/_renderer.py`. Entries marked `custom: true` describe the node groups
generated per imported structure; they have no class in the asset file, so
they are documented on `docs/api/generated_nodes.qmd` instead.
"""

import pathlib
import textwrap
from collections import defaultdict
import bpy
import yaml
from nodebpy.assets import _codegen
import molecularnodes as mn

DOCS_FOLDER = pathlib.Path(__file__).resolve().parent
QUARTO_YML = DOCS_FOLDER / "_quarto.yml"
CATS_FILE = pathlib.Path(mn.assets.MN_DATA_FILE).parent / "blender_assets.cats.txt"

BEGIN_MARKER = "  # -- begin generated node sections (docs/generate.py) --\n"
END_MARKER = "  # -- end generated node sections --\n"

NODES_DESC = (
    "The node groups included with Molecular Nodes, grouped by the same"
    " categories as the add menu inside of Geometry Nodes. Each node group is"
    " also available as a typed Python class for scripting."
)

CATEGORY_DESCRIPTIONS = {
    "Style": "Generate geometry for the different molecular styles",
    "Select": "Create boolean selections based on the atomic attributes",
    "Color": "Set and manipulate the colors of atoms",
    "Animate": "Animate values and geometry over frames",
    "Topology": "Work with the bond and residue topology of structures",
    "Attributes": "Read attributes from structures",
    "Density": "Work with volumetric density data",
    "DNA": "Nodes for working with oxDNA and other DNA models",
    "Ensemble": "Instance and manipulate ensembles and assemblies",
    "Simulation": "Nodes for simulations inside of Geometry Nodes",
    "Curves": "Create and manipulate curves",
    "Geometry": "General geometry processing utilities",
    "Utilities": "Small helper node groups used throughout Molecular Nodes",
    "Materials": "Shader node groups for the pre-built materials",
}

header = """---
toc: true
toc-depth: 3
fig-align: center
---
"""


def catalog_paths() -> dict[str, str]:
    "Map catalog UUID -> catalog path (e.g. 'Molecular Nodes/Style'), in file order."
    paths = {}
    for line in CATS_FILE.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line == "VERSION 1":
            continue
        uuid, path, _ = line.split(":", 2)
        paths[uuid] = path
    return paths


def categorise(path: str) -> tuple[str, str]:
    "Split a catalog path into (page category, subcategory within the page)."
    parts = path.split("/")
    if parts[0] == "Molecular Nodes":
        parts = parts[1:]
    if not parts or parts == [""]:
        return "General", ""
    return parts[0], "/".join(parts[1:])


def video_markdown(urls: list[str]) -> str:
    lines = []
    for url in urls:
        if not url.endswith(".mp4"):
            url += ".mp4"
        lines.append(f"![]({url})")
    return "\n\n".join(lines)


def class_path(group) -> str | None:
    "Qualified path of the generated class for a node group, or None if absent."
    module = _codegen._TREE_MODULES.get(group.bl_idname)
    if module is None or not hasattr(mn.nodes, module):
        return None
    cls_name = _codegen._class_name(group.name)
    if not hasattr(getattr(mn.nodes, module), cls_name):
        return None
    return f"nodes.{module}.{cls_name}"


def generate_node_sections() -> None:
    extras: dict[str, dict] = yaml.safe_load((DOCS_FOLDER / "nodes.yml").read_text())
    bpy.ops.wm.open_mainfile(filepath=str(mn.assets.MN_DATA_FILE))
    paths = catalog_paths()

    # (category, subcategory) -> list of class paths; subcategories get their
    # own page (e.g. 'Utilities / Color'). Page order follows cats.txt.
    pages: dict[tuple[str, str], list[str]] = defaultdict(list)
    for path in paths.values():
        pages.setdefault(categorise(path), [])

    documented = set()
    for group in bpy.data.node_groups:
        # skip non-assets and groups linked in from Blender's bundled libraries
        if group.asset_data is None or group.library is not None:
            continue
        path = class_path(group)
        if path is None:
            print(f"no generated class for asset node group: {group.name}")
            continue
        pages[categorise(paths.get(group.asset_data.catalog_id, ""))].append(path)
        documented.add(group.name)

    for name, extra in extras.items():
        if not extra.get("custom") and name not in documented:
            print(f"nodes.yml entry matches no asset in the .blend file: {name}")

    contents = []
    for (category, sub), items in pages.items():
        if not items:
            continue
        items.sort()
        name = f"{category} / {sub}" if sub else category
        slug = f"{category}_{sub}" if sub else category
        contents.append(
            {
                "kind": "page",
                "path": "nodes." + slug.lower().replace(" ", "_").replace("/", "_"),
                "summary": {
                    "name": name,
                    "desc": "" if sub else CATEGORY_DESCRIPTIONS.get(category, ""),
                },
                "contents": [{"name": path, "children": "embedded"} for path in items],
            }
        )
    section = {"title": "Nodes", "desc": NODES_DESC, "contents": contents}

    text = QUARTO_YML.read_text()
    head, found, rest = text.partition(BEGIN_MARKER)
    _, found_end, tail = rest.partition(END_MARKER)
    if not (found and found_end):
        raise RuntimeError(f"generated node section markers not found in {QUARTO_YML}")
    block = yaml.safe_dump([section], sort_keys=False, width=88)
    QUARTO_YML.write_text(
        head + BEGIN_MARKER + textwrap.indent(block, "  ") + END_MARKER + tail
    )


def generate_custom_nodes_page() -> None:
    extras: dict[str, dict] = yaml.safe_load((DOCS_FOLDER / "nodes.yml").read_text())
    with open(DOCS_FOLDER / "api/generated_nodes.qmd", "w") as file:
        file.write(header)
        file.write("# Generated Node Groups\n\n")
        file.write(
            "These node groups are generated for each imported structure, rather"
            " than being included in the asset file, so they have no corresponding"
            " Python class. The name of each generated node group starts with the"
            " prefix shown and ends with the name of the structure it was"
            " generated for.\n"
        )
        for name, extra in extras.items():
            if not extra.get("custom"):
                continue
            file.write(f"\n## {extra.get('label', name)} (`{name}`)\n")
            if extra.get("description"):
                file.write(f"\n{extra['description']}\n")
            if extra.get("videos"):
                file.write(f"\n{video_markdown(extra['videos'])}\n")


def generate_data_table() -> None:
    with open(DOCS_FOLDER / "data_table.qmd", "w") as file:
        file.write(header)
        file.write("# Data Tables\n\n")
        file.write(
            "The different lookup tables that are used to conver strings to integers in Molecular Nodes.\n\n"
            "Code for this can be found on the [GitHub Page](https://github.com/BradyAJohnston/MolecularNodes/blob/main/molecularnodes/assets/data.py)\n\n"
        )
        file.write(
            "### Residue Names\n\n| Name | Integer |\n|----------:|:------------|\n"
        )
        for name, res in mn.assets.data.residues.items():
            file.write(f"| {name} | `{res['res_name_num']}::Int` |\n")
        file.write("\n")
        file.write("\n")

        file.write(
            "### Atom Names\n\n| Name | Integer |\n|----------:|:------------|\n"
        )
        for name, value in mn.assets.data.atom_names.items():
            file.write(f"| {name} | `{value}::Int` |\n")
        file.write("\n")


if __name__ == "__main__":
    generate_node_sections()
    generate_custom_nodes_page()
    generate_data_table()
