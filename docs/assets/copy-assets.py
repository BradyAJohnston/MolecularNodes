"""Quarto pre-render hook: copy docs/assets/* into the build directory root.

Great Docs only auto-copies a root-level `assets/` directory; ours lives in
`docs/assets/` (see posit-dev/great-docs#329). `mn.css` and the logo are
wired through `site.css` / `logo` in great-docs.yml, so this is only needed for
files Quarto must find next to `_quarto.yml`, i.e. `filters.lua`.
"""

import os
import shutil
from pathlib import Path

build_dir = Path(os.environ.get("QUARTO_PROJECT_DIR", os.getcwd()))
repo_root = next(
    p for p in (build_dir, *build_dir.parents) if (p / "pyproject.toml").exists()
)

for src in (repo_root / "docs" / "assets").iterdir():
    if src.is_file() and src.suffix != ".py":
        shutil.copy2(src, build_dir / src.name)
