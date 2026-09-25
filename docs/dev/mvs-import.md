# MolViewSpec import

Plan for supporting `.mvsj` / `.mvsx` import ([#787](https://github.com/BradyAJohnston/MolecularNodes/issues/787)),
translating [MolViewSpec](https://molstar.org/mol-view-spec-docs/) scene descriptions
into Molecular Nodes entities and node setups.

## What MVS is, structurally

An `.mvsj` file is a JSON tree (`.mvsx` is a zip with `index.mvsj` plus
relative-URI assets). Nodes inherit context from their ancestors:

```
download (url)
└── parse (format)
    └── structure (model / assembly / …)
        └── component (selector)
            └── representation (type)
                ├── color (+ optional sub-selection)
                └── opacity
```

with scene-level `canvas`, `camera` / `focus`, and (later) `label` / `tooltip` /
`primitives` / `volume` nodes. The official `molviewspec` PyPI package provides
the pydantic models for the whole schema, validation, and MVSX packing.

## Central design idea: translate statically into attributes + style branches

We do not mirror MVS's tree with node-tree machinery. Instead:

1. **Each `structure` node → one `Molecule` entity** (`Molecule.load` on the
   downloaded file).
2. **Each `representation` under a `component` → one style branch**, built by
   `add_style`.
3. **Selections and colors are resolved in Python with numpy, not in nodes.**
   A `ComponentExpression` (or static selector like `polymer` / `ligand`)
   becomes a boolean mask computed against the universe, stored as a named
   boolean attribute; `add_style` accepts an attribute name as `selection`.
   MVS `color` nodes (several per representation, each with an optional
   sub-selection) are baked into a per-representation RGBA attribute which the
   style branch reads as its color, the same idiom as the import-time `Color`
   attribute. `opacity` is baked into that attribute's alpha channel, which the
   default material already wires into the shader's alpha.

This keeps the translation deterministic and testable without evaluating node
trees, and requires no node-asset changes. Dynamic alternatives (building
compare-node chains per expression) remain possible later but are not needed
for fidelity.

## Mapping table

| MVS | MolecularNodes |
|---|---|
| `download.url` | URL download with hashed cache name (plus relative paths, resolved against the `.mvsj` location) |
| `parse.format` (bcif / mmcif / pdb) | existing readers |
| `structure` model | `Molecule.load` |
| `structure` assembly | `add_style(assembly=True)` |
| `component` static selectors | import-time attributes: `is_peptide`, `is_nucleic`, `is_solvent`, `is_carb`, `is_hetero` (`polymer` = peptide ∨ nucleic, `ligand` ≈ hetero ∧ ¬solvent ∧ ¬carb) |
| `ComponentExpression` (auth/label asym / seq / comp ids, ranges, atom ids) | numpy masks → named boolean attribute; `label_*` fields resolved by a secondary biotite parse with extra fields, falling back to `auth_*` with a warning for formats without them |
| `representation`: cartoon, ball_and_stick, spacefill, surface | `StyleCartoon`, `StyleBallAndStick`, `StyleSpheres`, `StyleSurface` |
| putty, backbone, line, carbohydrate | nearest style + warning (putty/backbone ≈ ribbon, line ≈ sticks, carbohydrate ≈ ball_and_stick); dedicated support later |
| `color` (X11 names, hex) + sub-selection | baked RGBA attribute per branch |
| `color_from_source` / `color_from_uri` + palettes | phase 2: source column → attribute; categorical palettes → `custom_color_iswitch`, continuous → color ramp |
| `opacity` | alpha channel of the baked color attribute |
| `canvas.background_color` | `Canvas.world` background |
| `camera` / `focus` | `Canvas` camera + `look_at` |
| `transform` / `instance` | phase 2: Blender object transforms |
| `label` / `tooltip` | deferred (MN annotations are the natural target) |
| `volume` + isosurface / slice | phase 2/3 via `Grids` + density styles |
| multi-state `snapshots` / animation | phase 3 |

## Architecture

```
molecularnodes/entities/mvs/
  __init__.py   # public load()
  selectors.py  # static selector / ComponentExpression -> numpy boolean mask
  importer.py   # tree walk: inherited context -> per-structure recipes ->
                # Molecule entities + add_style branches; scene nodes -> Canvas
```

- **Depend on `molviewspec`** rather than re-typing the schema: small,
  pure-Python, officially maintained, pydantic-based. Pin it; the spec is
  versioned and evolving.
- **Fidelity policy: warn, never fail.** Unsupported node kinds and parameters
  are collected and reported once at the end of import, so any valid file
  yields a best-effort scene. This mirrors Mol*'s forward compatibility.
- **Entry points**: `mn.entities.mvs.load("scene.mvsj")` for scripting; a
  `FileHandler` for drag-and-drop of `.mvsj` into the viewport; a menu entry.
- **Testing**: `.mvsj` files are generated in the tests with the `molviewspec`
  builder against local structure files (relative URIs), so tests stay
  offline. The `molstar/mol-view-spec` repo's `test-data/` examples serve as
  reference material.

## Phases

1. **MVP** (this phase, implemented): `.mvsj` single-state → download (URL,
   relative and `file://` paths) / parse (bcif, mmcif, pdb) / structure
   (model + assembly) / components (static selectors + expressions) / the four
   core representations with nearest-style fallbacks / flat colors with
   sub-selections / opacity / canvas background / camera + focus.
2. `color_from_source` / `color_from_uri` + palettes, putty / carbohydrate,
   transforms & instances, `.mvsx` archives, volumes via `Grids`,
   `model_index` / `block_index` selection.
3. Labels / tooltips (on top of MN annotations), primitives, multi-state
   snapshots / animation (per-snapshot entity sets plus scene frame ranges).
