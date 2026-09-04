"""Generate the node documentation: GUI pages and the quartodoc API structure.

The .blend asset file is the source of truth: every node group marked as an
asset in `molecularnodes/assets/node_data_file.blend` is documented, grouped
by its asset catalog (the same grouping the GUI shows). Each node group is
documented twice, with the two versions linking to each other:

- a GUI page per category under `docs/nodes/` — node name, demo video,
  description and the input/output socket tables read straight from the node
  group interface — listed in the "Nodes" sidebar (written into the marked
  sidebar block of `docs/_quarto.yml`)
- an API reference page per category under `docs/api/reference/`, rendered by
  quartodoc from the generated node classes (the quartodoc section is written
  into the marked block of `docs/_quarto.yml` before `quartodoc build` runs);
  `docs/_node_links.yml` maps each class to its GUI entry so the renderer in
  `docs/_renderer.py` can link back

Extra information that cannot live on the nodes themselves (long-form prose,
demo videos) lives in `docs/nodes.yml`, keyed by node group name. Entries
marked `custom: true` describe the node groups generated per imported
structure; they have no class in the asset file, so they are documented on
`docs/nodes/generated_nodes.qmd` instead.
"""

import pathlib
import re
import textwrap
from collections import defaultdict
import bpy
import yaml
from nodebpy.assets import _codegen
import molecularnodes as mn

DOCS_FOLDER = pathlib.Path(__file__).resolve().parent
QUARTO_YML = DOCS_FOLDER / "_quarto.yml"
NODES_FOLDER = DOCS_FOLDER / "nodes"
LINKS_FILE = DOCS_FOLDER / "_node_links.yml"
CATS_FILE = pathlib.Path(mn.assets.MN_DATA_FILE).parent / "blender_assets.cats.txt"

BEGIN_MARKER = "  # -- begin generated node sections (docs/generate.py) --\n"
END_MARKER = "  # -- end generated node sections --\n"
SIDEBAR_BEGIN_MARKER = "  # -- begin generated nodes sidebar (docs/generate.py) --\n"
SIDEBAR_END_MARKER = "  # -- end generated nodes sidebar --\n"

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


def _anchor(name: str) -> str:
    "Stable heading anchor for a node's entry on its GUI page."
    slug = "".join(c if c.isalnum() else "-" for c in name.lower())
    return re.sub("-+", "-", slug).strip("-")


def _table_cell(text: str) -> str:
    return text.replace("|", "\\|").replace("\n", " ").strip()


# socket types docs/filters.lua knows how to color (as `Type::Type` tags)
SOCKET_KEYWORDS = {
    "Float", "Int", "Vector", "Geometry", "Bool", "Matrix", "Rotation",
    "Material", "Color", "Collection", "String", "Object", "Menu", "Shader",
    "Image", "Bundle", "Closure",
}  # fmt: skip


def _socket_type_code(socket_type: str) -> str:
    "Inline code for a socket type, tagged for coloring by docs/filters.lua."
    name = socket_type.removeprefix("NodeSocket")
    if name in SOCKET_KEYWORDS:
        return f"`{name}::{name}`"
    print(f"socket type not colorable by docs/filters.lua: {socket_type}")
    return f"`{name}`"


def _socket_table(group, in_out: str) -> str:
    "Markdown table of the group's input or output sockets, as named in the GUI."
    rows = []
    for item in group.interface.items_tree:
        if item.item_type != "SOCKET" or item.in_out != in_out:
            continue
        socket_type = _socket_type_code(item.socket_type)
        rows.append(
            f"| {item.name} | {socket_type} | {_table_cell(item.description)} |"
        )
    if not rows:
        return ""
    label = "Inputs" if in_out == "INPUT" else "Outputs"
    return "\n".join(
        [
            f"**{label}**",
            "",
            "| Name | Type | Description |",
            "|:-----|:-----|:------------|",
            *rows,
        ]
    )


def _gui_entry(group, cls_path: str, page_slug: str, extra: dict) -> str:
    "One node's entry on a GUI page: name, video, description, socket tables."
    parts = [f"## {group.name} {{#{_anchor(group.name)}}}"]
    if extra.get("videos"):
        parts.append(video_markdown(extra["videos"]))
    if group.description:
        parts.append(group.description)
    if extra.get("description"):
        parts.append(extra["description"])
    for in_out in ("INPUT", "OUTPUT"):
        table = _socket_table(group, in_out)
        if table:
            parts.append(table)
    api_href = f"/api/reference/nodes.{page_slug}.qmd#molecularnodes.{cls_path}"
    parts.append(f"API Reference: [`{cls_path}`]({api_href})")
    return "\n\n".join(parts)


def _write_marked_block(begin: str, end: str, block: str) -> None:
    "Replace the block between the given markers in _quarto.yml."
    text = QUARTO_YML.read_text()
    head, found, rest = text.partition(begin)
    _, found_end, tail = rest.partition(end)
    if not (found and found_end):
        raise RuntimeError(f"marker not found in {QUARTO_YML}: {begin.strip()!r}")
    QUARTO_YML.write_text(head + begin + block + end + tail)


def generate_node_sections() -> None:
    extras: dict[str, dict] = yaml.safe_load((DOCS_FOLDER / "nodes.yml").read_text())
    bpy.ops.wm.open_mainfile(filepath=str(mn.assets.MN_DATA_FILE))
    paths = catalog_paths()

    # (category, subcategory) -> list of (node group, class path); subcategories
    # get their own page (e.g. 'Utilities / Color'). Page order follows cats.txt.
    pages: dict[tuple[str, str], list] = defaultdict(list)
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
        key = categorise(paths.get(group.asset_data.catalog_id, ""))
        pages[key].append((group, path))
        documented.add(group.name)

    for name, extra in extras.items():
        if not extra.get("custom") and name not in documented:
            print(f"nodes.yml entry matches no asset in the .blend file: {name}")

    NODES_FOLDER.mkdir(exist_ok=True)
    contents = []
    sidebar_contents = [{"href": "nodes/index.qmd", "text": "Overview"}]
    links: dict[str, str] = {}
    for (category, sub), items in pages.items():
        if not items:
            continue
        items.sort(key=lambda item: item[1])
        name = f"{category} / {sub}" if sub else category
        slug = f"{category}_{sub}" if sub else category
        slug = slug.lower().replace(" ", "_").replace("/", "_")
        desc = "" if sub else CATEGORY_DESCRIPTIONS.get(category, "")

        # quartodoc page of the generated node classes (the API version)
        contents.append(
            {
                "kind": "page",
                "path": "nodes." + slug,
                "summary": {"name": name, "desc": desc},
                "contents": [
                    {"name": path, "children": "embedded"} for _, path in items
                ],
            }
        )

        # the GUI version of the same nodes, linking to the API version
        page = (
            f'---\ntitle: "{name}"\ntoc: true\ntoc-depth: 2\nfig-align: center\n---\n'
        )
        if desc:
            page += f"\n{desc}\n"
        for group, path in items:
            page += f"\n{_gui_entry(group, path, slug, extras.get(group.name, {}))}\n"
            links[f"molecularnodes.{path}"] = f"/nodes/{slug}.qmd#{_anchor(group.name)}"
        (NODES_FOLDER / f"{slug}.qmd").write_text(page)
        sidebar_contents.append({"href": f"nodes/{slug}.qmd", "text": name})

    sidebar_contents.append(
        {"href": "nodes/generated_nodes.qmd", "text": "Generated Node Groups"}
    )

    section = {"title": "Nodes", "desc": NODES_DESC, "contents": contents}
    _write_marked_block(
        BEGIN_MARKER,
        END_MARKER,
        textwrap.indent(yaml.safe_dump([section], sort_keys=False, width=88), "  "),
    )

    sidebar = [{"id": "nodes", "collapse-level": 2, "contents": sidebar_contents}]
    _write_marked_block(
        SIDEBAR_BEGIN_MARKER,
        SIDEBAR_END_MARKER,
        textwrap.indent(yaml.safe_dump(sidebar, sort_keys=False, width=88), "  "),
    )

    LINKS_FILE.write_text(
        "# Generated by docs/generate.py: node class -> its entry on the GUI node\n"
        "# pages, used by docs/_renderer.py to link the API pages back.\n"
        + yaml.safe_dump(links, sort_keys=True, width=88)
    )


def generate_custom_nodes_page() -> None:
    extras: dict[str, dict] = yaml.safe_load((DOCS_FOLDER / "nodes.yml").read_text())
    with open(NODES_FOLDER / "generated_nodes.qmd", "w") as file:
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
