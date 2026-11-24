# ProteinMotion vs Molecular Nodes: what to borrow

Local working report, 2026-09-18. Not committed. Sources: ProteinMotion v0.9.1 at commit
`1d91186` (2026-09-17, https://github.com/pdpppd/proteinmotion, MIT) read from a local clone,
and Molecular Nodes on `main` at `8af4b253`.

## 1. Summary

ProteinMotion (PM) is a ~7.6k-line Python package that renders protein movies with its own
wgpu/Metal renderer, driven by a Manim-style `add` / `play` / `wait` scene script. It has an
optional Blender EEVEE backend that streams meshes into a headless Blender for depth of field.
It is protein-only, file-only (no fetch, no DSSP, no nucleic acids, no GUI), and verified on
Apple silicon. Its strengths are the authoring model, first-class measurements (hydrogen bonds,
screened Coulomb contacts, distance rulers), per-residue staggered animation, synchronized 2D
plots, labels that avoid each other, contact-guided morphs between different proteins, and an
unusually disciplined testing and documentation culture.

Molecular Nodes (MN) is far broader in what it can load and draw, and Blender gives it real
materials, Cycles, keyframes and a GUI. What MN lacks, relative to PM, is almost entirely on the
*authoring* and *measurement* side: there is no Python way to say "fade the context, focus on
this helix, label it, draw its hydrogen bonds, and show the distance trace beside it" without
hand-building node trees, annotations and keyframes.

Three things are worth acting on first:

1. **Measurement nodes**: Hydrogen Bonds and Coulomb Contacts, with a dashed-ruler style. The
   Charge node from #1225 and Find Bonds already give most of the plumbing.
2. **Per-residue staggering as a field**: one node that turns frame range plus residue delay
   into a 0..1 factor per atom. It upgrades Set Color, opacity, Animate Frames and the symmetry
   `Factor` inputs at once.
3. **A timeline layer on `Canvas`**: eased keyframes for camera orbit, zoom, focus (depth of
   field on a selection), style parameters and annotation opacity, so a script reads as a
   storyboard and Blender's timeline provides exact seeking.

PM's MIT licence means algorithms can be ported with attribution in the module docstring.

## 2. What ProteinMotion is

| | ProteinMotion 0.9.1 |
| --- | --- |
| Renderer | wgpu (Metal on macOS); one WGSL shader with ribbon, bond, sphere, mesh, density entry points; MSAA 1 or 4; weighted-blended OIT for transparency; distance-fog depth cue; fixed studio lighting. Optional EEVEE via a persistent background Blender process. |
| Output | H.264/HEVC through PyAV (VideoToolbox on macOS, libx264 elsewhere); NV12 conversion on the GPU. A 24 s, 1080p60 NMR cartoon film exports in 3.7 s on an M3 Max. |
| Input | PDB/mmCIF via Gemmi, all models kept; `.npy` and MDAnalysis trajectories (`md` extra); MRC/CCP4 density. No fetch. Water and hydrogens dropped by default. |
| Representations | cartoon, ribbon (same sweep; cartoon↔ribbon is a shader blend), ball-and-stick (analytic spheres, bonds clipped at atom surfaces), surfaces vdW/SAS/approximate SES from a local SDF grid plus marching cubes. |
| Secondary structure | HELIX/SHEET records only; no DSSP. Override with a string. |
| Selection | `protein.select(chain=, residues=, atoms=)`, ANDed keyword filters, inclusive author-number ranges, `|` to union. Not a query language. |
| Styling | per-atom colour and opacity tracks; `Colorize` and `SetOpacity` with `residue_delay` (N→C stagger, reversible); `color_by(values, scale, thickness=)` maps per-residue values to colour and cartoon cross-section. |
| Animation | `Rotate`, `Morph`, `Deform`, `BackboneMorph`, `PlayTrajectory`, `Representation`, `Focus`, `FocusPull`, `FadeIn/Out`, `Write/Unwrite`; rate functions `linear`, `smooth` (mirrored quintic), `ease_in_out_sine`, `there_and_back`. Timeline is compiled: `seek(t)` restores snapshots and replays, so backward seeking is exact; two animations writing the same channel of one object in one `play` raise. |
| Measurements | `Distance` (2D overlay or depth-tested 3D dashed ruler with a caption gap), `HydrogenBonds` (D–A ≤ 3.5 Å, D–H–A ≥ 150°, virtual amide H, templates for side chains), `Electrostatics` (screened Coulomb, `formal` or PQR charges), all with `.highlight()` using a bounded pool of ruler slots. |
| Labels | `Text`, `Callout` (screen-fixed text with a leader to a region centroid), `ResidueLabel(s)` with automatic overlap avoidance over 10 candidate offsets, `Write` animation adapted from Manim. Sizes are 1080p design pixels. |
| Plots | `TimeSeriesPlot` (cursor follows the interpolated trajectory frame, `live_value` closures), `ContactMap` (live Cα map), `SequenceTrack`, `ColorLegend`. Drawn by the same vector overlay as text. |
| Density | `DensityMap.from_file`, `isosurface(level, units="sigma")`, `slice(axis, position)`, `crop(region)` keeping the parent mean/std; animatable level and slice. |
| Morphing | `Morph` by atom identity with Kabsch alignment; `BackboneMorph` for different proteins: contact-map matching (soft contacts, fingerprint seeds, monotone DP, bitset branch-and-bound clique search with a deadline and an honest optimality report), per-residue rigid translation, N→C delay, fades for unmatched residues. |
| Agent skill | `skills/proteinmotion-movies/` shipped inside the wheel with `proteinmotion install-skill`; gated references; explicit scientific-honesty rules. |
| Testing | 103 tests: geometry against exhaustive or analytic references, chirality checks, pixel-identical backward seeks, boundary continuity audits of a 100 s film, codec metadata via ffprobe. `docs/VALIDATION.md` is a dated log with hardware and scope for every number. |

## 3. Compare and contrast

| Dimension | ProteinMotion | Molecular Nodes | Verdict |
| --- | --- | --- | --- |
| Rendering | Own rasteriser, fast, no shadows/AO/refraction; EEVEE optional | Cycles and EEVEE, materials, world, compositor, passes | MN ahead on quality; PM ahead on iteration speed |
| Authoring model | Scene class, `play()` clips, rate functions, channel conflicts, exact seek | Node trees plus Blender keyframes; `Canvas.record()` loop; no easing helpers, no timeline API | PM ahead; the gap is an API layer, not capability |
| Geometry generation | CPU records + GPU sweep; Catmull–Rom cartoon; sheet arrowheads; cartoon↔ribbon blend | Geometry nodes; cartoon, ribbon, surface, sticks, spheres; nucleic cartoons; oxDNA; density styles | MN broader; PM's representation crossfade is the one idea to take |
| Secondary structure | File records only | File, MDAnalysis DSSP per frame with averaging, pure-GN Topology DSSP | MN ahead |
| Selection | Keyword filters on chain/residue/atom | MDAnalysis language, 30 Select nodes, callables, managed attributes | MN ahead; PM's inclusive author ranges and `|` union are nice ergonomics |
| Per-residue staggering | Built into every colour/opacity/morph animation | Only `Factor`/`Stagger` in the symmetry nodes and manual field maths | PM ahead; cheap to add as a node |
| Opacity | Per-atom alpha, OIT, bonds take `min(a,b)`, EEVEE layered approximation | `Default` reads per-point alpha from the Color attribute; `Flat` and `AmbientOcclusion` do not; no Set Opacity node | PM ahead; the gap is a node and two material recipes |
| Camera | Spherical target/distance, `frame()` by bounding sphere, `focus(region, follow=True)`, DOF on a selection, `FocusPull` | `look_at` with viewpoints and margin, `Camera.frame_points`; no focus tracking, no DOF helper | PM ahead on focus and DOF; MN's framing is now comparable |
| Measurements | Distance rulers, H-bonds, Coulomb, all animatable | `COMDistance`, dihedral annotations, in-node distances and angles; H-bond energy only inside DSSP | PM ahead |
| Labels | Callouts with leaders, residue labels with overlap avoidance, `Write` | `Label2D/3D`, `AtomInfo`, styled text and lines, PIL render path | MN has the primitives; missing callout, auto-placement and write-on |
| 2D plots | Four synchronized plot types built in | matplotlib templates only, not a dependency | PM ahead |
| Density | Iso, slice, crop, sigma units, animatable | Surface, wire, ISO with contours and slice; VDB volumes | Comparable; MN lacks `crop` around a region and sigma units in the API |
| Trajectories | Lazy MDAnalysis reader, `aligned()` view, 3-frame cache | Full MDAnalysis integration, subframes, interpolation, averaging, periodic correction, streaming IMD | MN ahead |
| Morphing | Same-topology `Morph`; different-protein `BackboneMorph` with contact matching | Animate Frames (same topology), Animate Dihedrals | PM ahead on different-protein morphs |
| Numerical properties | `ResidueValues` (B factor, streaming RMSF), `ColorScale`, thickness mapping | `Color Attribute Map`, `b_factor`, pLDDT; no RMSF | Comparable colouring; MN lacks RMSF and thickness mapping |
| Loading breadth | Proteins, local files | Proteins, nucleic acids, ligands, ensembles, cellPACK, starfile, oxDNA, density, fetch from RCSB/AlphaFold | MN far ahead |
| Platform | Apple silicon verified; other wgpu adapters untested | Linux, macOS, Windows | MN ahead |
| GUI | None | Full Blender add-on | MN ahead |
| Docs and testing culture | Writing guide, validation log, executable examples with verification JSON, continuity audits | Quarto docs, golden renders, snapshot tests | PM ahead on rigour; MN ahead on volume |
| Agent skill | Shipped in wheel, installer, gated references, honesty rules | `skills/` in #1226 | Adopt PM's conventions |

## 4. Nodes we could create

Ordered roughly by value over effort. "Plumbing" names existing MN pieces to reuse.

### 4.1 Hydrogen Bonds

- **What**: edges between donor and acceptor atoms with `hbond_distance` and `hbond_angle`
  attributes, using PM's criteria: D–A in (1.5, 3.5] Å, D–H–A ≥ 150°, different residues,
  most linear angle wins per pair. Inputs: Selection, Max Distance, Min Angle, Hydrogens
  menu (Auto / Explicit / Backbone), Include Side Chains.
- **Virtual amide H**: place H 1.01 Å from N along the normalised sum of unit vectors N→C(prev)
  and N→CA, skipping PRO and chain starts. This is its own small helper worth exposing
  (`Backbone Amide H`) because hydrogen-free structures are the norm.
- **Plumbing**: `find_bonds.py` already does a cutoff pair search in GN; `topology_dssp.py`
  has `HBondEnergy` and `CheckHBond`; `_shared/hydrogen_bonding_partner.py` knows side-chain
  partners; `backbone_n/c/ca` nodes give the atoms.
- **Pairs with**: a `Style Dashed Bonds` (4.3) so the output renders directly.

### 4.2 Coulomb Contacts

- **What**: pairs within a cutoff with `E = 332.06 q_i q_j exp(-r/λ) / (ε r)` kcal/mol
  as an edge attribute, sign for attractive/repulsive, `Min Energy` threshold. Inputs:
  Selection, Cutoff (12 Å), Dielectric (80), Screening Length (8 Å, 0 for unscreened),
  Min Energy, Exclude Same Residue, Exclude Bonded.
- **Plumbing**: the `charge` attribute now comes from the file or the lookup table (#1225),
  so the node reads `Charge` directly; PM's `formal` scheme (−1 over ASP/GLU carboxyl O,
  +1 on LYS NZ, ARG NH1/NH2) is a useful fallback menu option when charges are all zero.
- **Colouring**: attractive blue, repulsive red, magnitude to width, via the same dashed style.

### 4.3 Style Dashed Bonds / Distance Ruler

- **What**: renders edges as dashed cylinders with `Dash Count`, `Dash Ratio`, `Radius`,
  optional gap in the middle sized for a label, colour from an edge attribute, opacity as
  `min` of endpoint opacities. Works on H-bond and Coulomb outputs and on a user-picked
  pair of points.
- **Plumbing**: `points_of_edge`, `edge_length`, `curve_custom_profile`; `COMDistance`
  annotation for the caption.

### 4.4 Animate Stagger (per-residue factor)

- **What**: a field node returning 0..1 per atom from `Frame Start`, `Frame End`,
  `Residue Delay` (frames), `Reverse`, `Easing` menu (Linear / Smooth / Sine / There and
  Back). Rank residues along the chain (unique residue id order), `begin = start +
  rank·delay`, `span = length − (n−1)·delay`. Errors if span ≤ 0.
- **Why**: PM applies this to colour, opacity, morph progress and cartoon thickness. In MN
  one node feeds `Set Color` (mix factor), an alpha attribute, `Animate Frames` (`Frame`
  input per point), and the symmetry nodes' `Factor`.
- **Plumbing**: `animate_value.py` (SceneTime, smoother step), `unique_residue_id`,
  `relative_index`, `attribute_run`. Add the four rate functions as a shared
  `Rate Function` helper used by Animate Value too; the mirrored quintic avoids `1 + ε`
  at clip ends.

### 4.5 Style Highlight

- **What**: transparent sphere / box / per-atom halo around a selection. Inputs: Selection,
  Shape menu, Padding, Colour, Opacity, Line Width (box edges as tubes), Atom Scale.
  Sphere radius from the selection's enclosing sphere, box from the bounding box.
- **Plumbing**: `select_sphere` and `select_cube` do the inverse already; `centroid`,
  bounding box nodes; Transparent material.

### 4.6 Morph To Attribute

- **What**: moves `position` toward a stored vector attribute (`target_position`) with a
  factor, per-residue delay (4.4), and fades atoms flagged `unmatched` in or out. Same
  topology is Animate Frames' job; this is for different proteins where Python has computed
  the correspondence.
- **Python side**: port `match_backbones` (scipy only, MIT) as `mn.analysis.match_backbones`
  or an `mol.morph_to(other)` that writes `target_position` and `unmatched` attributes and
  applies the Kabsch alignment. Keep PM's honesty: the match report carries
  `search_completed`, `cardinality_proved`, `global_optimal`.

### 4.7 Cartoon thickness and width by field

- **What**: `Style Cartoon` gains `Thickness Factor` and `Width Factor` float fields
  (default 1) multiplied into the helix, sheet and loop dimensions, so B factor, RMSF or
  pLDDT can drive cross-section the way PM's `color_by(thickness=)` does.
- **Plumbing**: `CAToHelix`, `CAToSheet`, `CAToLoops` take scalar widths today; sample the
  field at the CA before the curve conversion.

### 4.8 Style Ribbon ↔ Cartoon blend

- **What**: a `Ribbon Factor` input on Style Cartoon that lerps helix and sheet cross
  sections toward the loop tube, so a representation transition is a shape morph rather
  than a cut. PM's cartoon and ribbon share one sweep for exactly this reason.
- **Difficulty**: medium; arrowheads and sheet smoothing need to fade with the factor.

### 4.9 Density Slice and Density Crop

- **What**: `Density Slice` as a coloured plane (axis, position 0..1, colour ramp, resolution)
  separate from the ISO node's slice-of-contours; `Density Crop` boxing a volume around a
  selection with padding while keeping the parent map's mean and standard deviation so
  sigma levels stay comparable. Sigma units as a menu on the existing density styles.
- **Plumbing**: `density_style_iso_surface.py` already has Slice Width/Center; VDB volume
  sampling nodes.

### 4.10 Depth cue

- **What**: distance fog toward the background: `smoothstep(start, end, depth)` mixed into
  the shader or the compositor mist pass. A `Canvas.depth_cue = 0.65` style property that
  wires the Mist pass in the compositor is probably enough; a shader-side version in the
  preset materials is the higher-quality option.

## 5. Feature sets we are missing

These are Python-API and annotation features rather than nodes.

### 5.1 A timeline layer on `Canvas`

PM's whole appeal is that a script reads as a storyboard. Blender already has the engine
(keyframes, F-curve interpolation, exact seeking, no conflict problem because a property has
one F-curve). What is missing is the sugar:

```python
with canvas.timeline(fps=60) as t:
    t.play(canvas.camera.orbit(theta=0.3), run_time=2)  # keyframes on the camera
    t.play(style.set("opacity", 0.06), run_time=1)  # keyframes on a node input
    t.play(
        canvas.camera.focus(helix, fstop=4), run_time=1.5
    )  # DOF focus object tracks selection
    t.wait(0.5)
    t.play(label.write(), run_time=0.8)  # annotation progress property
```

Pieces: rate functions applied through F-curve easing (Blender's Bezier handles or
`interpolation="SINE"` etc.), `orbit`/`zoom`/`shift` on `Camera` that keyframe location and
rotation, `set()` on node-group inputs, a `focus` that creates an empty parented to the
selection's centre of mass (the `COM` annotation already computes it) and sets
`camera.data.dof.focus_object` with an f-stop. Blender's DOF is better than PM's EEVEE layering
trick, so this is a place MN can immediately be ahead.

### 5.2 Measurements as Python objects

`mol.measure.distance(a, b)`, `mol.measure.hbonds(...)`, `mol.measure.contacts(...)` returning
arrays plus a `.show()` that adds the style from 4.3. MDAnalysis has
`HydrogenBondAnalysis` and `distances` for the numerical side; the node versions cover the
per-frame visual side. Salt bridges fall out of the Coulomb contacts.

### 5.3 Labels

- **Callout**: `Label2D` anchored to a selection's projected centroid with a leader line,
  `tip` dot/arrow, hidden when the anchor leaves the view, opacity following the target. The
  annotation system has every primitive (`draw_line_2d` with pointers, `draw_text_2d`).
- **Residue labels with overlap avoidance**: `mol.annotations.add_residue_labels(selection)`
  trying a handful of pixel offsets and choosing the least overlapping. PM's 10-candidate,
  minimise-summed-overlap heuristic is small and effective.
- **Write-on**: a `progress` property on text annotations, keyframeable, drawing outline then
  fill. PM adapted Manim's lag timing (MIT).
- **Design-pixel convention**: PM sizes text in 1080p pixels scaled by `height/1080`. MN's
  `text_size` should state its unit and scale with resolution, or renders at different sizes
  come out inconsistent.

### 5.4 Synchronized plots as annotations

Promote the matplotlib templates to annotation classes: `TimeSeries` (cursor mapped through
the entity's current universe frame, so subframes and interpolation are honoured),
`ContactMap` (live Cα map, `min_separation`), `SequenceTrack` (current residue colours, so it
animates with colour changes), `ColorLegend` (for `Color Attribute Map`). matplotlib as an
optional extra, PIL-drawn fallbacks for the simple ones.

### 5.5 Numerical properties

`mol.rmsf(selection="name CA", align=True)` via MDAnalysis, written as a per-residue attribute;
`mol.residue_values({(chain, resid): v}, name=)` for imported per-residue data; a `ColorScale`
helper that mirrors `Color Attribute Map` so the legend and the node agree. Thickness mapping
via 4.7.

### 5.6 Opacity as a first-class channel

`Default` already wires the `Color` attribute's alpha into the BSDF through the `MN Color`
shader group (`nodes/materials/default.py:19`); `Flat` and `AmbientOcclusion` do not read it.
So the missing pieces are a `Set Opacity` node writing alpha (with the stagger factor from 4.4
as its natural input), alpha support in the other two presets, and a documented blend-mode
expectation per engine. PM's rule that a bond takes `min(alpha_a, alpha_b)` and that
sticks are clipped at sphere surfaces so interiors never show during fades are both worth
copying into Style Ball and Stick.

### 5.7 Trajectory alignment view

`mol.aligned(reference=0, selection="name CA")` applying a Kabsch fit per frame on load, the
equivalent of PM's `Trajectory.aligned()`. MDAnalysis `AlignTraj` does the work; MN needs the
switch on the entity so it composes with subframes and interpolation.

### 5.8 `doctor` and skill installation

- `molecularnodes doctor` (or `uv run -m molecularnodes doctor`): Blender/bpy version,
  GPU and Cycles devices, whether an installed extension's wheels shadow the venv
  (the `BLENDER_USER_EXTENSIONS` trap), compositor device, ffmpeg presence.
- Ship `skills/` inside the wheel with an `install-skill` command that copies to
  `~/.claude/skills` or `~/.codex/skills`, idempotent and atomic, as PM does.

## 6. Where MN should not follow PM

- **A second renderer.** PM's renderer exists because it has no Blender; MN's answer to
  "fast preview" is EEVEE and lower samples, not a rasteriser.
- **Keyword-only selection.** MDAnalysis phrases and Select nodes are strictly more capable.
  Adopt only the ergonomics: inclusive author-numbered ranges in `select_res_id_string`
  (already close) and `|` on selection objects.
- **Compiled timeline with snapshots.** Blender keyframes already give exact seeking.
- **File records for secondary structure.** MN's DSSP paths are better.
- **Protein-only assumptions.** Every node above must handle nucleic acids and ligands or
  state clearly that it does not.

## 7. Process and documentation practices to adopt

- **A validation log.** `docs/VALIDATION.md` records, per version, the machine, inputs, frame
  counts, timings and the scope of each claim, and keeps negative results. MN's golden renders
  prove images are stable, not that numbers are right; a dated log for the analysis features
  (DSSP agreement, H-bond counts on an ideal helix, symmetry RMSD to deposited assemblies)
  would make PR claims checkable.
- **Reference fixtures.** PM's 16-residue ideal helix that must yield exactly 12 i→i+4
  hydrogen bonds, with a reversed-H negative control, is the pattern for 4.1.
- **Continuity audits.** For animation features, assert frames ±1 µs across a clip boundary
  are pixel-identical. Cheap with the 32 px render fixture.
- **A writing guide.** PM's `docs/writing-guide.md` (Google developer style, simplified
  technical English, no slogans, benchmarks with hardware and scope) is short and would suit
  the great-docs migration.
- **Executable docs examples with verification artefacts.** Every PM docs example is rendered
  by a script and its numbers are tied to JSON produced by real runs.
- **Epistemic labelling.** PM states everywhere that morphs and NMR interpolation are
  illustrations, hydrogen bonds are geometric, Coulomb is screened and charge-dependent, and
  the morph search reports whether it proved optimality. MN's docs for the same features
  should carry the same sentences.
- **Skill conventions.** Gate reference files ("read only the relevant reference"), forbid
  invented API, warn against keeping a starter template's residue numbers on a different
  protein, and state units. Our two skills in #1226 should adopt the reference-gating layout
  as they grow.

## 8. Suggested order

| Step | Items | Why first |
| --- | --- | --- |
| 1 | 4.4 Animate Stagger, 4.3 Dashed style, 4.1 Hydrogen Bonds | Small nodes, immediate visual payoff, reuse existing helpers, unblock 4.2 and 5.2 |
| 2 | 5.1 timeline layer with camera orbit/zoom/focus and DOF | Largest gap in how MN scripts read; Blender's DOF beats PM's approximation |
| 3 | 4.2 Coulomb Contacts, 5.3 Callout and residue labels | Completes the "explain an interaction" story |
| 4 | 4.7 thickness by field, 5.5 RMSF and residue values, 5.4 plots | Data-to-picture features |
| 5 | 4.6 Morph To Attribute plus `match_backbones` port, 4.8 ribbon blend | Bigger, benefits from the stagger and opacity work |
| 6 | 5.8 doctor and install-skill, section 7 practices | Cheap, improves everything else |

## 9. Verification notes

- All PM claims come from reading the clone (`src/proteinmotion/*.py`, `docs/*.md`,
  `skills/`, `tests/`); none from running it. The EEVEE backend, matching and rate functions
  were read directly; the rest through a structured survey of the source.
- MN claims come from the current `main`, including the node catalogue (~253 assets),
  annotation classes, `DSSPManager`, `Canvas` and `FrameRecorder`. Alpha handling in 5.6 was checked
  by reading the material recipes, not by rendering.
