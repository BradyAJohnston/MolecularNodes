"""Build the documentation site with Great Docs.

    uv run docs/build.py build            # full build -> great-docs/_site/
    uv run docs/build.py build --no-refresh

This is a thin wrapper around the `great-docs` CLI that first registers a
`docstring_parsed` hook. The hook injects the extras from `docs/nodes.yml`
(long-form prose and demo videos that cannot live on the nodes inside the
.blend file) into the API reference page of each generated node class, after
the class description and before the Parameters section. Everything else is
plain `great-docs build`; `great-docs preview` needs no wrapper.
"""

import pathlib
import griffe as gf
import yaml
from great_docs.cli import main
from great_docs.hooks import on_docstring_parsed

DOCS_FOLDER = pathlib.Path(__file__).resolve().parent
NODE_MODULES = "molecularnodes.nodes."


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


EXTRAS = _extras_by_class()


def _extra_markdown(extra: dict) -> str:
    text = ""
    if extra.get("description"):
        text += f"{extra['description'].strip()}\n"
    for url in extra.get("videos", []):
        if not url.endswith(".mp4"):
            url += ".mp4"
        text += f"\n![]({url})\n"
    return text.strip()


@on_docstring_parsed
def inject_node_extras(obj, sections):
    """Insert the nodes.yml prose and videos into a node class docstring."""
    if not obj.path.startswith(NODE_MODULES) or not obj.is_class:
        return sections
    extra = EXTRAS.get(obj.name)
    if extra is None:
        return sections
    block = gf.DocstringSectionText(_extra_markdown(extra))
    # after the leading description, before Parameters
    index = 1 if sections and isinstance(sections[0], gf.DocstringSectionText) else 0
    return [*sections[:index], block, *sections[index:]]


if __name__ == "__main__":
    main()
