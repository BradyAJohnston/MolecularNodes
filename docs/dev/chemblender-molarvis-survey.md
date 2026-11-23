# ChemBlender and molar_vis: features worth borrowing (2026-10-17)

A survey of two neighbouring projects for features Molecular Nodes lacks.

- **ChemBlender.** The original is [sphereman/ChemBlender_2_0](https://github.com/sphereman/ChemBlender_2_0) (HEAD 9077096, about 5.9k lines of Python plus about 85 `CH_*` node groups). The fork is [psiQAQ/ChemBlender_2_x](https://github.com/psiQAQ/ChemBlender_2_x) (v2.4.0, HEAD 54ecf4c, about 195k lines). Both are GPL-3.0.
- **molar_vis.** [yesint/molar_vis](https://github.com/yesint/molar_vis) (HEAD 6e81433) is a standalone Rust/wgpu viewer built on MolAR, licensed Artistic-2.0. Only its algorithms are relevant here.

Everything below comes from reading the source. No code from either project was run. Gaps marked † were checked against MN `main` at 5f759d40.

## ChemBlender

The fork carries MN as a git submodule, but only as reference: no Python file imports `molecularnodes`. Its docs say it borrowed ideas from `entities/density/grids.py`, `nodes/geometry.py` and `session.py`, and rejected the pickle-based `.MNSession` storage. One fork doc (`docs/quantum-visualization/2.3.0/references.md`) wrongly lists MN as MIT-licensed.

Most of the fork is process scaffolding: 43 architecture decision records (ADRs), a `.cbq/` sidecar project folder of NPY arrays and a JSON manifest, provenance, and quality badges. The parts worth borrowing are a few readers and its volume and vector visualisation.

### Worth doing

| Feature | Where in ChemBlender | MN today | Notes |
|---|---|---|---|
| Coordination polyhedra | `CH_Coord Polyhedra` + `CH_Remove Coplanar Edges` (original, `Chem_Nodes_En.blend`); `AddCoordPolyhedraButton` | none | For each centre, Proximity ligands within an RMin–RMax shell, then a ConvexHull, inside a For-Each zone. About 25 nodes. Suits metalloprotein sites (Fe-S, Zn) as well as crystals. |
| Colour a style by a second volume | fork `surface_view.py::ensure_property_surface_cache`: GridToMesh(density) → SampleGrid(property) → StoreNamedAttribute | †no `SampleGrid` anywhere in `molecularnodes/nodes` | The MN form is APBS `.dx` electrostatics on `Style Surface`; we already read `.dx`. The NCI variant colours an RDG isosurface by sign(λ2)ρ. |
| Isosurface fixes | fork `create_signed_isosurfaces` negates the grid to get the negative lobe | †`density_style_iso_surface.py:181` clamps `Threshold` to 0–1; line 256 meshes the same grid at `-threshold` | Lift the clamp for orbital and ESP values. Check the normals on the negative lobe; this is untested. Add a default isovalue as a fraction of max\|v\| (`grid_semantics.py::default_grid_isovalue`). |
| Gaussian Cube reader + full affine grid transform | fork `core/cube.py` (~150 lines, no dependencies); `grid_cache_service.py::_transform_matrix` → `openvdb.createLinearTransform` | †`grids.py:286` only `preScale` + `postTranslate` | Non-orthogonal grids come out wrong today. |
| Structure export | fork `core/exporters/` (xyz, pdb, pqr, mol2, cif, poscar, cube, MOL/SDF) with loss previews | †no structure writer anywhere | `Molecule.write(path, selection=)` via biotite or MDAnalysis writers, using current positions. |
| Normal-mode animation | fork `vibration_view.py` (Python frame handler, `co = ref + d·A·sin φ`); `phonon_frames.py` | none | Do it in geometry nodes instead: a stored displacement attribute plus `Position + d·A·sin(2π·frame/period)`, with no handler. ProDy/MDAnalysis ANM modes are the biology source. |
| Per-atom vector arrows | fork `dataset_view.py::_ensure_vector_arrow_group` (AlignEulerToVector + cone instances) | †`ts.forces` / `ts.velocities` not exposed | Expose them as optional per-frame attributes and add a vector style built on `primitive_arrow`. |
| Colour-bar legend | fork `render_annotations.py::create_render_annotations` | †none; only the Rg chart template | A legend annotation for attribute-coloured figures. |
| PubChem / SMILES import | original `read.py::download_sdf_from_pubchem`, `scaffold.py` (RDKit ETKDGv3 + UFF); fork `core/derivations/smiles_3d.py` | †no RDKit or PubChem | PubChem fetches SDF, which our existing reader handles. SMILES → 3D needs RDKit as an optional dependency. |

### Maybe

- **Readers.** extXYZ cell and per-atom properties, which the MDAnalysis XYZ reader drops (fork `core/formats/extxyz.py`). A POSCAR reader. MOL2/PQR/GRO in `read_structure` and the FileHandler (`ui/ops.py:42`).
- **Bond drawing.** Dashed aromatic bonds offset towards each ring's centre (a precomputed `ring_center` attribute), plus a per-bond `dashed` flag for H-bond or coordination bonds in Ball & Stick and Sticks.
- **Periodic Array.** Clipping to fractional bounds and keeping molecules whole across the boundary (original `CH_Supercell`, `CH_Cell Expand`).
- **Lattice-direction framing** along a, b, c, [uvw] or (hkl). The original's (hkl) uses the direct lattice vectors where it needs the reciprocal ones, so it is wrong for non-orthogonal cells; don't copy it.
- **Volume rendering.** A two-colour Principled Volume material for signed fields (fork `scientific_materials.py::volume_material`) and a colour-mapped slice plane (`grid_sampling.py`).
- **Thermal ellipsoids,** which need ADPs parsed in `corecif`.
- **Import tooling.** A pre-import info popup showing atoms, frames, cell, and whether bonds were read or guessed. A reader registry that sniffs content rather than matching on suffix (fork `reader_api/`). Session storage as NPY plus a JSON manifest instead of a pickle.

### Skip

Interactive building and minimisation, crystallographic symmetry editing, the provenance and quality system, the separate worker process for optional dependencies, band structure, DOS, Fermi surfaces, critic2, in-scene spectra, and the generic mesh generators. Measurements, selection, render setup and element colours are already covered as well or better in MN.

## molar_vis: framing and camera animation

`framing.py` solves how close the camera should be for a given orientation and deliberately does not choose the orientation. molar_vis's main contribution is choosing it.

### Unobstructed view (`crates/molar_vis_core/src/unobstructed.rs`)

`vis.unobstructed_view(rep, zoom_out=1.0, resolution=256)` turns the camera so that the most target atoms are directly visible.

- **Visibility test.** Spheres are projected orthographically. A target atom is hidden when a nearer sphere covers its projected centre. Nearer target atoms count too, so self-overlap is penalised. A 2D cell grid with cells of 2× the largest radius keeps the test cheap.
- **Scorers.** Mesh representations occlude with their real triangles, tested by rays against a BVH. Small targets switch to a depth rasteriser that counts visible *area* (`unobstructed/pixels.rs`). A GPU version of the scorer also exists.
- **Search.** A 256-direction Fibonacci sphere, then three refinement rings of 16 directions at 10°, 4° and 1.5° around the best direction (`search_directions_with_seeds`).
- **Orientation.** `look_along_quat` keeps world-up upright, falling back to another axis when the direction is near-vertical.
- **Framing.** `Camera::focus_visual_bounds` fits centre ± radius in the new basis at the viewport aspect ratio. `zoom_out` widens the framing about the target to show its surroundings.

### Animation (`animation.rs`, `docs/ANIMATION.md`)

`vis.movie()` takes the steps `rotate`, `zoom`, `focus`, `unobstructed_view`, `go_to_view`, `hold` and `move_view`, plus `movie.trajectory(mol, start, end, fps)`. That shape matches what the canvas-timeline prototype (#1238) already has. Interpolation details worth checking against #1238:

- **Rotation.** Orientation is slerped; #1238 keys Euler angles with `compatible`. Euler interpolation between arbitrary viewpoints does not follow the great-circle path and can wobble. This is untested in MN.
- **Distance.** Distance is interpolated in log space, so zoom speed feels constant. #1238 is linear.
- **Clipping.** During a move, molar_vis keeps the wider of the start and end clip ranges. #1238's `_PoseMove` only raises `clip_end` when the end pose needs it. Could a move between two tight framings clip partway through? Untested.
- **Orbits.** They are evaluated from the start pose, so a 720° orbit keeps both turns. #1238's single-Euler-component `Orbit` already does this.
- **Rotation centre.** It can be separate from the frame centre (VMD's `c` key, `Camera::set_center`). #1238's `orbit` could take a `pivot=` argument.
- **Saved poses.** `view_json()` / `go_to_view(json)` store a pose without touching lighting.

### Other methods (#1171)

Issue #1171 ("Automatic viewpoint selection for look_at()") lists:

- Viewpoint entropy (Vázquez et al. 2003).
- Perceptual preference (Secord et al. 2011).
- PCA of the point set.
- The Bonaventura et al. 2018 survey.
- `framing.enclosing_sphere()` for turntables that don't "breathe" (it is implemented but unused).

#1227 ("look_at with a point array and a custom viewpoint frames against the stale camera basis") affects any automatic viewpoint that frames a `get_view()` array.

## Prototype: `molecularnodes/viewpoint.py` (branch `unobstructed-view`)

The prototype is a pure-numpy module that is independent of `framing.py`. Every method returns a basis as rows `(right, up, forward)`, so it composes with `fit_camera_to_points`. All the searched methods share `search()`, a 256-direction Fibonacci sphere followed by three refinement rings around the best 3 seeds (molar_vis refines only 1).

| Method | Scorer | Notes |
|---|---|---|
| `pca` | closed form: look along the least-variance axis, longest axis horizontal, signs fixed (look in from −Y, up +Z) | Blind to occlusion. Unstable for near-degenerate spreads (see the oxDNA planar-PCA issue). |
| `unobstructed` | molar_vis centre test: count target atoms whose projected centre no nearer sphere covers. Prefilters occluders to the target's footprint, then does a pairwise comparison per direction | Ties are broken towards the largest projected target bounding box. |
| `visible_area` | orthographic depth raster (`rasterize`) over the target's enclosing sphere at a fixed world scale; counts winning target pixels | About 40× slower than `unobstructed` in numpy. |
| `entropy` | Vázquez viewpoint entropy over labelled parts (per atom by default), background included | Needs no target, so it suits whole structures. |
| `raycast` (Canvas only) | `BVHTree` ray from each target atom towards the camera against the evaluated meshes | Occludes with the real cartoon. Point-cloud styles are not meshes and are ignored. |

Wired up as `canvas.look_at(target, viewpoint="unobstructed" | "visible_area" | "raycast" | "pca" | "entropy", occluders=None)`. Occluders default to every visible mesh, point-cloud or curves object.

Occluder handling:
- Rendered vertices are thinned to one per 1.7 Å cube.
- Points within 2 Å of the target are dropped. Without this, the target's own rendering hides every atom from every direction.

Also included is the fix for #1227: `frame_points` now updates the view layer before reading the basis.

### Results (cartoon protein + ball-and-stick ligand, Cycles 16 samples, 480×360, margin 0.1)

Visible = ligand atoms with a clear ray to the camera through the rendered cartoon. This is the `raycast` objective, so it is biased towards that method.

| Structure | Method | Time (s) | Visible |
|---|---|---|---|
| 1STP (biotin) | front / pca / unobstructed / visible_area / raycast | 0.01–1.2 | 16/16 for all |
| 1STP | entropy (whole structure) | 2.4 | 13/16 |
| 1cbs (retinoic acid in a β-barrel) | front | 0.01 | 10/22 |
| 1cbs | pca | 0.01 | 14/22 |
| 1cbs | unobstructed | 0.03 | **22/22** |
| 1cbs | visible_area | 1.3 | **22/22** |
| 1cbs | raycast | 0.02 | **22/22** |
| 1cbs | entropy (whole structure) | 4.2 | 18/22 |

- **Recommendation (two structures only; superseded by the gallery below).** `unobstructed` is as good as the more faithful scorers and fast (about 30 ms), so it should be the default. `raycast` is the best fit when the scene is all meshes. `visible_area` adds nothing here at about 40× the cost; keep `rasterize` for `entropy`. `pca` is a cheap first guess for whole structures.
- **1STP.** Biotin is visible from every direction, so the methods differ only by tie-break.

### Gallery: six structures (`docs/dev/viewpoint_gallery.py`)

The script renders pocket and whole-protein framing for 1cbs, 1STP, 1ATN, 8H1B, 4PNO and 6ZPA. Results below are in pocket mode; "clear" counts ligand atoms with a clear ray to the camera through the rendered cartoon.

- **`raycast`** has every atom clear on all six structures, and it also looks best by eye.
- **`unobstructed`** reaches every atom only on 1cbs and 1STP. It falls behind on 1ATN (23/31), 8H1B (23/27), 6ZPA (27/31) and 4PNO (27/29). Using each cartoon vertex as a 1.7 Å sphere makes the ribbon look thicker than it renders. Try a smaller occluder radius, or vertex normals, before making it the default.
- **`visible_area`** sits between the two (26–29 clear on the harder cases) at about 40× the cost.
- **`pca`** and **`entropy`** often hide the ligand. They suit whole-structure framing, not pockets.

Golden snapshots: `test_render_viewpoint[front|pca|unobstructed|visible_area|raycast|entropy]` in `tests/test_render_images.py` (1cbs, 128 px).

### Open questions

- **API shape.** Should automatic methods be more `viewpoint=` strings, or a separate `viewpoint="auto", method=...`? Should `occluders` be explicit or default to the scene?
- **Framing radii.** `get_view()` returns atom centres, so selection framing ignores radii. molar_vis frames centre ± radius.
- **Context.** There is no `zoom_out` for pocket context yet. `margin` works, but a pocket selection is closer to what molar_vis does.
- **Timeline.** Add a timeline clip (#1238 `LookAt`) that solves `unobstructed` at play time.
- **Turntables.** Use `enclosing_sphere` framing so they don't breathe (#1171).
- **Perceptual preference** (Secord 2011) is not prototyped. It needs fitted weights over several features.
- **Ensembles.** An entropy label per residue or chain would suit ensembles better than per atom.
