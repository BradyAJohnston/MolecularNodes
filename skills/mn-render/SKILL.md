---
name: mn-render
description: How to render images and animations of molecules from a Python script with the molecularnodes package (Canvas, engines, framing, styles, materials, colour, trajectories, deterministic renders for tests). Load before writing or debugging any script that calls mn.Canvas, snapshot, animation or record.
---

# Rendering with Molecular Nodes from a Python script

## 1. Mental model

- **One `Canvas`, many entities.** `mn.Canvas` wraps the active Blender scene: engine,
  resolution, camera, world lighting, compositor and output. Molecules are
  `mn.Molecule` entities that you load, style and then frame with `canvas.look_at`.
- **Create the Canvas before loading anything.** It loads the Molecular Nodes scene
  preset (camera, lights, world, compositor) into an empty scene. Calling `mn.Canvas()`
  again is safe: it binds to the existing scene and only reloads the preset when the
  scene holds no molecules.
- **Styles are node trees.** `add_style` builds a small node branch per call; the
  render shows whatever the evaluated geometry is. `look_at` frames that evaluated
  geometry, so a molecule with one chain styled is framed on that chain.
- **Units.** Positions are Å times `world_scale = 0.1`, so 1 Å is 0.1 Blender units.
- There is no `mn.Trajectory` any more; MD topologies and trajectories load into
  `mn.Molecule` as well.

## 2. Environment for a script

Always via `uv run` (never bare `python`). In a script that imports `bpy` directly:

```python
import os, tempfile

os.environ.setdefault(
    "BLENDER_USER_EXTENSIONS", tempfile.mkdtemp()
)  # before import bpy
import bpy  # mathutils is importable only after this
import molecularnodes as mn

canvas = mn.Canvas(mn.scene.Cycles(samples=64, device="CPU"), resolution=(1200, 900))
canvas.compositor.device = (
    "CPU"  # Blender 5 defaults the compositor to GPU; aborts without one
)
```

Without the `BLENDER_USER_EXTENSIONS` line, an installed MN extension's bundled wheels
shadow the venv. `tests/conftest.py` does the same thing. Keep scratch scripts outside
the repo.

## 3. Minimal still

```python
mol = mn.Molecule.fetch("4ozs")  # .bcif from RCSB, cached
mol.add_style("cartoon", material=mn.material.AmbientOcclusion())
canvas.look_at(mol, viewpoint="front")
canvas.snapshot("4ozs.png")  # returns an IPython Image for PNG/JPEG
```

`snapshot(path=None, frame=None, file_format="PNG", render_scale=100)`. With no path
it renders to a temp file and returns the image (notebooks display it). `frame=` renders
one scene frame and restores the current frame afterwards. TIFF/EXR write the file but
return `None`. The canonical docs versions are in `docs/api/index.qmd` and
`docs/api/rendering.qmd`.

## 4. Canvas settings

| member | meaning |
| --- | --- |
| `engine` | `mn.scene.EEVEE(samples=64, raytracing=True)` or `mn.scene.Cycles(samples=256, device="GPU", denoise=True, denoise_gpu=True)`; also accepts `"EEVEE"` / `"CYCLES"`. Constructor values are written to the scene on assignment. |
| `resolution` | `(x, y)` pixels. `render_scale` is the percentage. |
| `samples` | samples on the active engine. |
| `transparent` | film transparency (alpha background). |
| `background` | world background RGBA; shortcut into the `MN_world_shader` node. |
| `world.hdri_strength` | lighting strength of the preset HDRI. |
| `view_transform`, `exposure`, `gamma`, `look` | colour management. Blender defaults to AgX; use `"Standard"` when colours must come out as specified. |
| `passes` | render passes: combined, z, mist, normal, position, vector, diffuse_color, emit, environment, ambient_occlusion. |
| `frame`, `frame_start`, `frame_end`, `frame_range`, `fps` | timeline. Setting `frame` uses `frame_set`, so trajectories update. |
| `camera` | `lens`, `clip_start`, `clip_end`, `rotation` (XYZ Euler in degrees), `basis`, `set_viewpoint(...)`, `frame_points(points, margin)`. |
| `compositor` | `device`, `precision`, denoise settings, `reset()`, `clear()`, `add_annotations()`. |

- `canvas.clear()` deletes every object except cameras and lights (including the
  preset backdrop), purges orphans recursively, and keeps engine, world, compositor and
  render settings. `canvas.load_preset()` reloads the whole preset.
- `Cycles(device="GPU")` tries OPTIX, CUDA, METAL, HIP, ONEAPI and falls back to CPU
  with a warning. EEVEE needs a GPU; the GitHub runners do not have one, so CI-facing
  scripts use Cycles CPU.
- `world.reset()` and `compositor.reset()` remove the preset nodes: after them
  `canvas.background` raises `ValueError` and the annotation overlay is gone
  (`compositor.add_annotations()` restores it).

## 5. Loading entities

```python
mol = mn.Molecule.fetch("9MD2")  # format=".bcif", database="rcsb"
af = mn.Molecule.fetch("Q8W3K0", database="alphafold")
pdb = mn.Molecule.load("path/protein.pdb", style="cartoon")  # single structure
traj = mn.Molecule.load("topol.tpr", "traj.xtc", name="md")  # MD topology + coordinates
u_mol = mn.Molecule(u)  # from an MDAnalysis Universe
```

`fetch(code, format=".bcif", cache=download.CACHE_DIR, database="rcsb")`. `load` routes
single files through biotite and topology plus coordinates through MDAnalysis; `style=`
is optional and defaults to no style, leaving the tree empty. In tests use
`cache=data_dir` so downloads land in `tests/data`.

## 6. Styles, materials, colour

```python
mol.add_style(
    "surface", selection="protein", material=mn.material.Flat(), color="common"
)
mol.add_style("ball_and_stick", selection="not protein", color=(1.0, 0.5, 0.0, 1.0))
mol.add_style(
    lambda: mg.StyleCartoon(quality=5, loop_radius=0.6), color=lambda: mg.ColorRainbow()
)
```

- **Style names:** `spheres`, `cartoon`, `ribbon`, `surface`, `sticks`, `ball_and_stick`.
  Extra kwargs go to the style node's inputs (`quality`, `scale`, `loop_radius`,
  `sphere`, ...); unknown names raise `TypeError`, unknown style strings `ValueError`.
  `add_style` returns the molecule, so calls chain.
- **Callable style** (`lambda: mg.StyleX(...)`) builds the node itself, so it cannot be
  combined with `selection`, `material` or kwargs (`TypeError`); `color`, `assembly`
  and `name` are fine.
- **Selection** is an existing boolean attribute name first, then an MDAnalysis phrase
  (stored as a managed `sel_N` attribute), an `AtomGroup`, or a callable returning a
  boolean socket (`lambda: mg.IsPeptide() & mg.IsSideChain()`).
- **`assembly=True`** instances the style over the file's biological assembly.
- **Materials.** Presets in `mn.material`: `Default(roughness, ao_distance, ao_exponent)`,
  `AmbientOcclusion(distance, exponent)`, `Flat(outline, threshold, thickness)`,
  `Squishy(subsurface_scale, roughness)`, `Transparent(transparency, fresnel, outline_color)`.
  Each instantiation builds an independent datablock (`"Flat.001"`), and its parameters
  are live: `mat.distance = 0.15` changes the next render. Passing the string
  `material="Flat"` reuses one shared datablock, so edits affect every style using it.
  `"Flat Outline"` exists only as `mn.material.append_material("Flat Outline")`.
  `Transparent` needs the blended surface method, which its recipe sets. To swap the
  material on an existing style: `style.i.material.default_value = mat.material`.
- **Colour.** `color=` takes `"common"`/`"default"` (elements, carbons random per
  chain), `"plddt"`, the name of an existing *colour* attribute, an RGBA sequence, or a
  callable returning a colour socket. Any other string warns and applies nothing (this
  replaced a silent black render). Arrays go on directly with `mol["Color"] = rgba_array`.
- **Spheres in EEVEE.** Point clouds are ray-traced only by Cycles; on other engines
  `add_style` switches `StyleSpheres` from `sphere="Point"` to `"Instance"` for you.
  If you build the tree by hand, set `sphere="Instance"` or `"Mesh"` yourself, and use
  instances whenever point clouds with different materials are joined.

## 7. Framing

```python
canvas.look_at(mol)  # keep current direction, fit the subject
canvas.look_at(
    mol, viewpoint="top", margin=0.15
)  # front, back, top, bottom, left, right, default
canvas.look_at(mol, viewpoint=(90, 0, 45))  # XYZ Euler degrees
canvas.look_at(mol.get_view("chainID A and resid 1-40"))
canvas.look_at(points_xyz)  # any (N, 3) array in world units
```

- `look_at(target, viewpoint=None, margin=0.05)` fits the camera to the evaluated
  geometry without changing where it points; `margin=0` is exact, negative crops.
- `camera.lens` changes after `look_at` invalidate the fit: set the lens first, or call
  `look_at` again.
- **Stale camera basis.** `set_viewpoint` writes `rotation_euler`, but the fit reads
  `matrix_world`, which only updates with the depsgraph. `look_at(entity_or_object,
  viewpoint=...)` is safe because reading the evaluated geometry updates it first.
  `look_at(points_array, viewpoint=...)`, or `set_viewpoint(...)` followed by a separate
  `look_at`, fits against the previous orientation and renders empty or off-centre
  (measured 2026-09-18). Call `bpy.context.view_layer.update()` between the two until
  `look_at` does it itself.
- `frame_points` extends `clip_end` when the subject would be clipped; if a hand-placed
  camera renders nothing, check `camera.clip_end` first.
- Rotating the subject instead of the camera: leave the camera at `"default"` and
  rotate the object, or use the `viewpoint` tuple in a loop for orbits (section 9).

## 8. Trajectories

Per-molecule playback properties, all stored on the Blender object: `frame`,
`subframes`, `offset`, `average`, `correct_periodic`, `interpolate`. A trajectory with
`n` universe frames and `subframes = s` spans `n * (s + 1) - 1` scene frames.

```python
traj = mn.Molecule(u)
traj.dssp.init()  # secondary structure per frame, if cartoon/ribbon needs it
traj.subframes, traj.interpolate, traj.average = 1, True, 1
traj.add_style("cartoon", color="common")
canvas.fps = 24
canvas.frame_range = (0, u.trajectory.n_frames * 2 - 1)
canvas.frame = 0
canvas.look_at(traj, viewpoint="top", margin=0.35)
```

`Molecule.frames_to_collection(start=0, stop=None, step=1)` bakes frames into a
collection for the *Animate Frames* node when the timeline route is not wanted.
`set_frame` is what the frame-change handler calls; user code sets `canvas.frame`.

## 9. Animations and image sequences

Three routes, pick by what changes between frames:

1. **Timeline only changes** (trajectory playback, keyframed values, `SceneTime` in a
   tree): `canvas.animation(path="out.mp4", frame_start=None, frame_end=None,
   render_scale=100, fps=None, format=None)`. Format is inferred from a `.gif` suffix;
   MP4 is H.264; GIF needs pillow. It locks the interface, renders to a temp directory
   and restores frame, fps, scale and output settings afterwards.
2. **Camera, style or colour change per frame:**
   ```python
   with canvas.record("orbit.mp4", fps=12) as movie:
       for i in range(36):
           canvas.look_at(mol, viewpoint=(90, 0, i * 10))
           movie.render(render_scale=50)
   ```
   `record(path=None, fps=None, render_scale=100, frames_dir=None, overwrite=True)`
   returns a `FrameRecorder`; `render()` writes `%05d.png`, `finalize(path, fps, format)`
   assembles them. The context manager finalises only on a clean exit and only when a
   path was given. `frames_dir=` plus `overwrite=False` resumes an interrupted run.
   All frames must share one resolution.
3. **A few stills at chosen frames:** `canvas.snapshot(path, frame=f)` in a loop.

## 10. Deterministic renders and tests

- Same-seed Cycles CPU renders are pixel-identical on one platform. Recipe:
  ```python
  canvas.engine = mn.scene.Cycles(samples=32, device="CPU", denoise=False)
  canvas.compositor.device = "CPU"
  ```
  Compare decoded pixels, not PNG bytes (metadata differs). Exactness does not hold
  across platforms or Blender versions.
- **Golden tests** live in `tests/test_render_images.py` with the `golden_canvas`
  fixture (128 x 128, Cycles 256 samples CPU, no denoise, compositor CPU). Pattern:
  fetch, `add_style`, `golden_canvas.look_at(mol, viewpoint=...)`, then
  `assert image_snapshot == _render(golden_canvas, tmp_path)`. Use
  `assembly_image_snapshot` (15 % failing pixels allowed) for full-frame assemblies and
  thin ribbons; the default is 1 % of pixels over 4/255, mirroring Blender's own
  `render_report`. Goldens are in `tests/__snapshots__/test_render_images/`; update
  with `uv run pytest tests/test_render_images.py -k <name> --snapshot-update` and look
  at the PNGs before committing. Failures write `<test>.received.png` and an amplified
  `<test>.diff.png` to `tests/image_failures/` (CI uploads them).
- For tests that only need a render to succeed, copy the `render_canvas` fixture in
  `tests/test_canvas.py`: 32 x 32, Cycles with 1 sample on CPU.
- Whole-structure rotations between platforms are usually a degenerate orientation in
  the tree (principal-component alignment with equal eigenvalues), not render noise.
  Fix the geometry; do not loosen tolerances.

## 11. Renders for pull requests

The 128 px goldens are too small to show a change. Write a script outside the repo with
`mn.Canvas(mn.scene.EEVEE(), resolution=(1200, 900))` (or Cycles CPU on a machine
without a GPU), style, `look_at`, `snapshot(path)`, and attach the PNGs with
`gh pr create --attach` as described in the `mn-nodes` skill, section 11.

## 12. Where to read more

- `docs/api/canvas.qmd`, `docs/api/rendering.qmd`, `docs/api/materials.qmd`,
  `docs/api/styles.qmd`, `docs/api/trajectories.qmd`, `docs/api/mdanalysis.qmd`
  (full trajectory to MP4 walkthrough), `docs/api/blender.qmd` (colour management and
  interface locking).
- Source: `molecularnodes/scene/base.py` (Canvas), `scene/camera.py`, `scene/engines.py`,
  `scene/recorder.py`, `scene/world.py`, `scene/compositor.py`, `material.py`,
  `entities/molecule/base.py` (`add_style`, `fetch`, `load`, frame properties).
- Tests as examples: `tests/test_canvas.py`, `tests/test_render_images.py`,
  `tests/test_add_style.py`, `tests/test_material.py`.
