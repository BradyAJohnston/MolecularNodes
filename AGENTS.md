# Molecular Nodes: notes for coding agents

You are working in the Molecular Nodes repository, a Blender extension and Python package for importing, styling, animating and rendering structural biology data. The person driving you is the author of any contribution; read `AI_POLICY.md` before doing anything that will reach GitHub.

## Rules that apply to you

- **Never post to GitHub on the user's behalf.** Do not open issues or pull requests, comment, review, or reply to reviewers. Summarise in the conversation and let the user write their own words. `gh pr create` is acceptable only when the user has explicitly asked for it and has read the body.
- **Do not draft the user's AI disclosure.** Every PR has an "AI disclosure" section that the human writes. Do not write it for them and do not paste this session's summary into it.
- **Do not add yourself as an author.** No `Co-authored-by`, `Assisted-by`, `Generated-by` or similar trailers, and no "Generated with ..." footers, in commit messages or PR bodies. Disclosure goes in the PR's disclosure section, in the human's words.
- **Do not invent.** No made-up APIs, node names, socket names, test results or citations. If you did not run something, say so explicitly.
- **Do not alter, weaken or delete tests to make a change pass.** Do not touch `tests/__snapshots__` unless the user asked for a snapshot update and has looked at the result.
- **Keep diffs minimal.** Do not reformat, rename or refactor code you were not asked to touch. Models love to; reviewers hate it.
- **Do not use AI on `good first issue` tickets.** They are for humans to learn on.

## Where to look

| Question | Read |
| --- | --- |
| How the project is structured, how to build, test and write docs | `CONTRIBUTING.md` |
| The AI policy in full | `AI_POLICY.md` |
| Writing, building and testing node groups | `skills/mn-nodes/SKILL.md` |
| Rendering from a script: Canvas, framing, styles, materials, animation | `skills/mn-render/SKILL.md` |
| Design notes, surveys, in-progress plans | `docs/dev/` |
| User-facing docs, API reference sources, node prose | `docs/` (`docs/api/*.qmd`, `docs/nodes.yml`, `docs/attributes.qmd`) |
| Changelog | `docs/changelog.qmd` (add a line for user-visible changes) |
| Build, dependency and tool config | `pyproject.toml` (`[tool.extbpy]`, `[tool.nodebpy.assets]`, `[tool.ruff]`) |
| CI | `.github/workflows/ci.yml` |
| nodebpy (the node-building library) | https://bradyajohnston.github.io/nodebpy |

## Layout

```
molecularnodes/
  entities/      Molecule (biotite / MDAnalysis backed), density, ensemble, mvs, reload
  nodes/         node groups as nodebpy Python: geometry/, geometry/_shared/, shader/, materials/, compositor/
  assets/        resources.blend (tracked); nodes.blend is BUILT from nodes/ and gitignored
  scene/         Canvas, camera, engines, world, compositor, recorder
  ui/            Blender operators, panels, properties, preferences
  blender/       low-level bpy helpers (mesh, collections)
  material.py, color.py, framing.py, download.py, session.py, handlers.py
tests/           pytest suite; snapshots in tests/__snapshots__; small structures cached in tests/data
docs/            Quarto site; docs/generate.py regenerates node pages from the built blend
skills/          agent skills (also symlinked from .claude/skills/)
```

Key facts:

- **Node groups are Python.** `molecularnodes/nodes/**/*.py` are the source of truth; `assets/nodes.blend` is a build product. Edit the `.py`, then `build`, then `dump`, and read the diff. Comments inside `_build_group` do not survive `dump`. Details in the `mn-nodes` skill.
- **Positions are Å times 0.1.** `Molecule.world_scale = 0.1`, so 1 Å is 0.1 Blender units.
- **There is no `mn.Trajectory`.** MD topologies and trajectories load into `mn.Molecule`.
- **Operators are thin.** Logic lives in plain functions that operators call, so everything works from a script as well as the GUI.
- **Python is pinned** to what the target Blender ships (`requires-python` in `pyproject.toml`). Always run things through `uv run`, never bare `python`.

## Commands

```bash
uv sync --all-extras                      # dev environment (bpy, test and docs extras)
uv run -m nodebpy.assets build            # nodes/*.py -> assets/nodes.blend (~1 min)
uv run -m nodebpy.assets dump             # assets/nodes.blend -> nodes/*.py
uv run -m nodebpy.assets check            # CI check: build -> dump reproduces sources
uv run pytest tests/test_nodes.py -k name # run only what you changed; the suite runs `ensure` itself
uv run pytest -n auto                     # full suite (slow; renders happen on CPU)
uv run ruff format && uv run ruff check   # formatting and lint (pre-commit runs these)
uv run docs/generate.py                   # regenerate node docs pages from the built blend
uvx extbpy sync                           # make molecularnodes/ loadable by Blender from source
uvx extbpy build -p current               # build the extension zip for this platform
```

If `biotite` complains about CCD files under `uv run`, run `uv run python -m biotite.setup_ccd`.

## Gotchas

- In any script that imports `bpy` directly, set `BLENDER_USER_EXTENSIONS` to a temp directory before the import, or an installed MN extension's bundled wheels shadow the venv. `tests/conftest.py` shows the pattern.
- Create `mn.Canvas()` before loading molecules so the scene preset loads. Blender 5 defaults the compositor to GPU; set `canvas.compositor.device = "CPU"` on machines without one.
- `mathutils` imports only after `bpy`.
- Nodes must be created inside a tree context (`with mol.tree.reset() as (atoms, join):` in tests).
- Golden-image tests compare decoded pixels with Blender's own tolerances. Cross-platform mismatches on whole-structure orientation are usually a degenerate principal-component alignment in the tree, not render noise. Fix the geometry, do not loosen tolerances.
- Multi-input socket link order is not preserved by the asset build when the source is a group input. Route it through another node first.
- Keep scratch scripts and downloaded files outside the repo.

## When you are done

Tell the user, in plain terms: what you changed, what you ran and what its output was, what you could not run, and anything you are unsure about. They need that to write the PR description and disclosure themselves.
