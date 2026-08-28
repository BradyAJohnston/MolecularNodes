"""Generate the per-category node documentation pages in docs/nodes/ and the
lookup-table page in docs/user_guide/.

The .blend asset file is the source of truth: every node group marked as an
asset in `molecularnodes/assets/node_data_file.blend` is documented, grouped
by its asset catalog (the same grouping the GUI shows). Socket names,
defaults, tooltips and descriptions all come from the .blend interface —
the same data that generates the typed classes in `molecularnodes/nodes/`
(see generate_node_classes.py).

Extra information that cannot live on the nodes themselves (long-form prose,
demo videos) is merged in from `docs/nodes.yml`, keyed by node group name.
"""

import pathlib
from collections import defaultdict
import bpy
import yaml
from nodebpy.assets import _codegen
import molecularnodes as mn

DOCS_FOLDER = pathlib.Path(__file__).resolve().parent
NODES_FOLDER = DOCS_FOLDER / "nodes"
DATA_TABLE_FILE = DOCS_FOLDER / "user_guide" / "21-data-table.qmd"
CATS_FILE = pathlib.Path(mn.assets.MN_DATA_FILE).parent / "blender_assets.cats.txt"

def header(title: str, **extra: str) -> str:
    lines = [f"title: {title}"] + [f"{k}: {v}" for k, v in extra.items()]
    return "---\n" + "\n".join(lines) + "\ntoc: true\ntoc-depth: 3\nfig-align: center\n---\n"


def catalog_paths() -> dict[str, str]:
    "Map catalog UUID -> catalog path (e.g. 'Molecular Nodes/Style')."
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


def format_default(socket: _codegen._Socket) -> str:
    if socket.default in ("None", ""):
        return ""
    return f"`{socket.default}`"


def socket_tables(cls: _codegen._AssetClass) -> str:
    text = ""
    if cls.inputs:
        text += "\n\n| Input | Type | Default | Description |\n|---|---|---|---|\n"
        for s in cls.inputs:
            doc = s.description
            if s.menu_items:
                options = ", ".join(f"`{item}`" for item in s.menu_items)
                doc = f"{doc} Options: {options}".strip()
            text += f"| {s.name} | `{s.socket_class.removesuffix('Socket')}` | {format_default(s)} | {doc} |\n"
    if cls.outputs:
        text += "\n\n| Output | Type | Description |\n|---|---|---|\n"
        for s in cls.outputs:
            text += f"| {s.name} | `{s.socket_class.removesuffix('Socket')}` | {s.description} |\n"
    return text


def node_markdown(group, extra: dict) -> str:
    cls = _codegen._introspect_group(group, group.name, "")
    text = f"## {group.name}\n"
    description = group.description.strip()
    prose = (extra.get("description") or "").strip()
    if description:
        text += f"\n{description}\n"
    if prose and prose != description:
        text += f"\n{prose}\n"
    if extra.get("videos"):
        text += f"\n{video_markdown(extra['videos'])}\n"
    text += socket_tables(cls)
    return text


def custom_markdown(name: str, extra: dict) -> str:
    text = f"## {extra.get('label', name)}\n"
    text += (
        f"\n*Node groups with the prefix `{name}` are generated for each"
        " imported structure, rather than being included in the asset file.*\n"
    )
    if extra.get("description"):
        text += f"\n{extra['description']}\n"
    if extra.get("videos"):
        text += f"\n{video_markdown(extra['videos'])}\n"
    return text


def generate_node_pages() -> None:
    extras: dict[str, dict] = yaml.safe_load((DOCS_FOLDER / "nodes.yml").read_text())
    bpy.ops.wm.open_mainfile(filepath=str(mn.assets.MN_DATA_FILE))
    paths = catalog_paths()

    pages: dict[str, list] = defaultdict(list)
    documented = set()
    for group in bpy.data.node_groups:
        # skip non-assets and groups linked in from Blender's bundled libraries
        if group.asset_data is None or group.library is not None:
            continue
        category, sub = categorise(paths.get(group.asset_data.catalog_id, ""))
        pages[category].append((sub, group))
        documented.add(group.name)

    for name, extra in extras.items():
        if extra.get("custom"):
            pages[extra["category"]].append(("Generated Node Groups", name))
        elif name not in documented:
            print(f"nodes.yml entry matches no asset in the .blend file: {name}")

    for category, items in pages.items():
        items.sort(key=lambda item: (item[0], getattr(item[1], "name", item[1])))
        with open(NODES_FOLDER / f"{category.lower()}.qmd", "w") as file:
            file.write(header(category))
            current_sub = ""
            for sub, item in items:
                if sub != current_sub:
                    current_sub = sub
                    file.write(f"\n# {sub}\n")
                if isinstance(item, str):
                    file.write("\n" + custom_markdown(item, extras[item]) + "\n")
                else:
                    file.write(
                        "\n" + node_markdown(item, extras.get(item.name, {})) + "\n"
                    )


def generate_data_table() -> None:
    with open(DATA_TABLE_FILE, "w") as file:
        file.write(header("Data Tables", **{"guide-section": "Reference"}))
        file.write("\n")
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
    generate_node_pages()
    generate_data_table()
