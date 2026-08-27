"""Custom quartodoc renderer for the node class pages.

Injects the extras from `docs/nodes.yml` (long-form prose and demo videos of
the node in use, which cannot live on the nodes inside the .blend file) into
the rendered page for each node class — after the docstring description,
before the Parameters section.

Configured in `_quarto.yml` via `quartodoc: renderer: style: _renderer.py`.
"""

import pathlib
import re
import yaml
from plum import dispatch
from quartodoc import MdRenderer, layout
from quartodoc._griffe_compat import docstrings as ds
from quartodoc.renderers.md_renderer import ParamRow

DOCS_FOLDER = pathlib.Path(__file__).resolve().parent
MARKER = "<!-- injected from docs/nodes.yml -->"


def _class_name(name: str) -> str:
    "Mirrors nodebpy.assets._codegen._class_name (avoids importing bpy here)."
    cleaned = "".join(c if c.isalnum() or c.isspace() else " " for c in name)
    parts = cleaned.split()
    cleaned = "".join(p[:1].upper() + p[1:] for p in parts)
    if cleaned and cleaned[0].isdigit():
        cleaned = "_" + cleaned
    return cleaned or "AssetGroup"


def _extras_by_class() -> dict[str, dict]:
    extras = yaml.safe_load((DOCS_FOLDER / "nodes.yml").read_text())
    return {
        _class_name(name): extra
        for name, extra in extras.items()
        if not extra.get("custom")
    }


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
                ParamRow(name.strip(), "", annotation=f"`{annotation.strip()}`")
            )
    return rows


def _extra_markdown(extra: dict) -> str:
    text = f"\n{MARKER}\n"
    if extra.get("description"):
        text += f"\n{extra['description']}\n"
    for url in extra.get("videos", []):
        if not url.endswith(".mp4"):
            url += ".mp4"
        text += f"\n![]({url})\n"
    return text


class Renderer(MdRenderer):
    style = "markdown_node_extras"

    _extras = _extras_by_class()

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
    def render(self, el: layout.DocClass):  # noqa: F811 (plum multiple dispatch)
        text = super().render(el)
        if not el.obj.path.startswith("molecularnodes.nodes."):
            return text
        extra = self._extras.get(el.obj.name)
        if extra is None:
            return text
        block = _extra_markdown(extra)
        # insert before the first docstring section (Parameters), so the prose
        # and videos read as part of the class description
        match = re.search(r"\n## .*\{\s*\.doc-section", text)
        if match:
            index = match.start()
            return text[:index] + "\n" + block + text[index:]
        return text + "\n" + block
