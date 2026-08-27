"""Regenerate the typed node classes in `molecularnodes/nodes/`.

Introspects the asset-marked node groups in
`molecularnodes/assets/node_data_file.blend` and writes one module per tree
type (`geometry.py`, `shader.py`, ...) with full numpy-style docstrings pulled
from the node group and socket descriptions in the .blend file. Those
docstrings are the single source of truth for the node documentation — they
feed editor tooltips and the quartodoc API reference.

Run after any change to the .blend file:

    uv run generate_node_classes.py

Uses nodebpy's internal codegen because the public
`generate_asset_modules` does not yet expose the `docstrings` option.
"""

import pathlib
import subprocess

from nodebpy.assets import _codegen
from nodebpy.builder import PackageLibrary

ROOT = pathlib.Path(__file__).resolve().parent
NODES_DIR = ROOT / "molecularnodes" / "nodes"

# The anchor file only serves to resolve the relative .blend path the same way
# the generated modules will at runtime (relative to molecularnodes/nodes/).
library = PackageLibrary(
    str(NODES_DIR / "geometry.py"), "../assets/node_data_file.blend"
)

classes = _codegen._introspect(library, None)
written = []
for tree_idname, module in _codegen._TREE_MODULES.items():
    tree_classes = [c for c in classes if c.tree_idname == tree_idname]
    if not tree_classes:
        continue
    path = NODES_DIR / f"{module}.py"
    path.write_text(
        _codegen._render_module(tree_classes, nodebpy_pkg="nodebpy", docstrings=True),
        encoding="utf-8",
    )
    written.append(path)
    print(f"wrote {path.relative_to(ROOT)}: {len(tree_classes)} classes")

subprocess.run(["ruff", "format", *(str(p) for p in written)], check=True)
