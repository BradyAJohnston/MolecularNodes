"""Custom quartodoc renderer for the node class pages.

The node groups are documented twice: GUI pages under `docs/nodes/` (generated
by `docs/generate.py`, holding the long-form prose and demo videos) and the
quartodoc class pages this renderer produces. For each node class this renders
the Inputs / Outputs docstring sections as tables, titles the entry with the
node group's name, and injects a link to the node's entry on the GUI pages —
using the mapping `docs/generate.py` writes to `docs/_node_links.yml`.

Configured in `_quarto.yml` via `quartodoc: renderer: style: _renderer.py`.
"""

import ast
import pathlib
import re
import yaml
from plum import dispatch
from quartodoc import MdRenderer, layout
from quartodoc._griffe_compat import docstrings as ds
from quartodoc.renderers.md_renderer import ParamRow

DOCS_FOLDER = pathlib.Path(__file__).resolve().parent
LINKS_FILE = DOCS_FOLDER / "_node_links.yml"
MARKER = "<!-- link to the GUI node documentation (docs/_renderer.py) -->"


def _gui_links() -> dict[str, str]:
    if not LINKS_FILE.exists():
        return {}
    return yaml.safe_load(LINKS_FILE.read_text()) or {}


# socket types docs/filters.lua knows how to color (as `Type::Type` tags);
# the nodebpy types map onto them by dropping the "Socket" suffix of the
# accessor types (`GeometrySocket`) or the "Input" prefix of the parameter
# annotations (`InputGeometry`)
_SOCKET_KEYWORDS = {
    "Float", "Int", "Vector", "Geometry", "Bool", "Matrix", "Rotation",
    "Material", "Color", "Collection", "String", "Object", "Menu", "Shader",
    "Image", "Bundle", "Closure",
}  # fmt: skip
_SOCKET_ALIASES = {"Boolean": "Bool", "Integer": "Int"}

# the interlink-style annotation links quartodoc renders for the Parameters
# table; nodebpy has no interlinks inventory so these would render as dead
# links stripped to plain code
_ANNOTATION_LINK_RE = re.compile(r"\[(Input(\w+))\]\(`nodebpy\.types\.\1`\)")


def _socket_code(text: str, type_name: str) -> str:
    "Inline code for a socket type, tagged for coloring by docs/filters.lua."
    keyword = _SOCKET_ALIASES.get(type_name, type_name)
    if keyword in _SOCKET_KEYWORDS:
        return f"`{text}::{keyword}`"
    return f"`{text}`"


def _annotation_code(annotation: str) -> str:
    return _socket_code(annotation, annotation.removesuffix("Socket"))


def _tag_annotation_links(text: str) -> str:
    return _ANNOTATION_LINK_RE.sub(lambda m: _socket_code(m.group(1), m.group(2)), text)


def _socket_rows(text: str) -> list[ParamRow]:
    """Parse the numpydoc-style entries of an Inputs / Outputs docstring
    section (``i.atoms : GeometrySocket`` with indented description lines)."""
    rows = []
    for line in text.splitlines():
        if line.startswith(" ") and rows:
            rows[-1].description = f"{rows[-1].description} {line.strip()}".strip()
        elif " : " in line:
            name, _, annotation = line.partition(" : ")
            rows.append(
                ParamRow(
                    name.strip(), "", annotation=_annotation_code(annotation.strip())
                )
            )
    return rows


def _node_title(el: layout.Doc) -> str | None:
    """The node group's name (the generated class's `_name` attribute) for
    generated node classes, used as the page heading instead of the class path."""
    if not el.obj.path.startswith("molecularnodes.nodes."):
        return None
    member = el.obj.members.get("_name")
    if member is None or member.value is None:
        return None
    try:
        title = ast.literal_eval(str(member.value))
    except (ValueError, SyntaxError):
        return None
    return title if isinstance(title, str) else None


class Renderer(MdRenderer):
    style = "markdown_node_extras"

    _links = _gui_links()

    @dispatch
    def render(self, el: ds.DocstringSectionAdmonition):
        # the Inputs / Outputs sections of the generated node class docstrings
        # are parsed by griffe as admonitions; render them as tables like the
        # Parameters section instead of as raw text
        if el.title in ("Inputs", "Outputs"):
            rows = _socket_rows(el.value.description)
            if rows:
                return self._render_table(
                    rows, ["Name", "Type", "Description"], "attributes"
                )
        return super().render(el)

    @dispatch
    def render_header(self, el: layout.Doc):  # noqa: F811 (plum multiple dispatch)
        title = _node_title(el)
        if title is None:
            return super().render_header(el)
        return f"{'#' * self.crnt_header_level} {title} {{ #{el.obj.path} }}"

    @dispatch
    def render(self, el: layout.DocClass):  # noqa: F811 (plum multiple dispatch)
        text = super().render(el)
        if not el.obj.path.startswith("molecularnodes.nodes."):
            return text
        text = _tag_annotation_links(text)
        title = _node_title(el)
        if title:
            # the docstring summary line duplicates the heading; drop it
            text = re.sub(rf"^{re.escape(title)}$\n?", "", text, count=1, flags=re.M)
        href = self._links.get(el.obj.path)
        if href is None:
            return text
        block = f"\n{MARKER}\n\nGUI reference: [{title or el.obj.name}]({href}).\n"
        # insert before the first docstring section (Parameters), so the link
        # reads as part of the class description
        match = re.search(r"\n## .*\{\s*\.doc-section", text)
        if match:
            index = match.start()
            return text[:index] + "\n" + block + text[index:]
        return text + "\n" + block
