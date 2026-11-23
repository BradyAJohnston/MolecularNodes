# Illustrate Rendering Pipeline - Technical Reference

Pseudocode summary of David Goodsell's Illustrate renderer (derived from the original Fortran
source). Non-photorealistic renderer for biomolecular illustration. No traditional lighting model
(no Phong, no Lambertian). Depth and form come entirely from conical soft shadows, depth fog, and
edge-detection outlines. Colors are flat per-atom-type assignments.

---

## 1. Shading Model: Flat / Unlit

No surface normals calculated and no light sources. Each atom gets a flat RGB color (0.0–1.0 per
channel) based on its type. 3D appearance comes from:
- Pre-rasterized sphere shapes (implicit curvature from the projection itself)
- Shadows and outlines doing all the heavy lifting for depth perception

**Blender equivalent:** Use an Emission or flat Diffuse shader with no lighting contribution, or a
Shader-to-RGB node to fully control shading.

---

## 2. Conical Soft Shadows (AO Substitute)

Primary depth-cueing technique. Functionally behaves like screen-space ambient occlusion.

### Algorithm

- For each pixel, sample a neighborhood (±50 pixels, stepping by 5)
- For each sample, compare Z-depth of the sample vs. the current pixel
- If the sample is above the current pixel (closer to camera) by more than a threshold (`rcone`),
  and the horizontal distance × cone angle < the vertical distance, it contributes shadow
- Each qualifying sample subtracts a small constant (`pcone`) from a shadow accumulator (starts at 1.0)
- The accumulator is clamped to a minimum (`pshadowmax`), preventing full blackness

### Parameters

| Param | Typical | Effect |
|-------|---------|--------|
| `pcone` | 0.0023 | Shadow contribution per occluder — higher = darker |
| `coneangle` | 2.0 | Cone tightness — higher = more localized shadows |
| `rcone` | 1.0 Å | Minimum Z-gap to cast shadow — removes artifacts in tight crevices |
| `pshadowmax` | 0.2–0.7 | Maximum darkening clamp |

### Blender Equivalent

Use the AO node in shader nodes, or a screen-space AO pass in compositing. `coneangle` maps to AO
radius/distance, `pcone` to intensity, and `pshadowmax` to a clamp on the AO factor. A Pointiness or
Cavity map baked from geometry could also work.

---

## 3. Outlines (3 Independent Layers)

All outlines use screen-space post-processing on the Z-buffer or ID buffers, not geometry-based.

### A. Contour Outlines (Depth-Edge Detection)

Applies a second-derivative (Laplacian) kernel on the Z-buffer. Four kernel options available:

1. **Kernel 1** — 3×3 Laplacian:
   ```
   -0.8  -1.0  -0.8
   -1.0  +7.2  -1.0
   -0.8  -1.0  -0.8
   ```
2. **Kernel 2** — 5×5 extended Laplacian (smoother)
3. **Kernel 3** — 3×3 Z-difference accumulation (accumulates normalized |Z_center - Z_neighbor| values)
4. **Kernel 4** — 5×5 Z-difference accumulation (excludes corners)

Result is thresholded with a ramp (`l_low` to `l_high`) to convert to opacity. An additional
smoothing pass: if ≥6 of 9 neighbors in a 3×3 grid have outline signal, average them; otherwise use
center value only.

**Parameters:**

| Param | Typical | Effect |
|-------|---------|--------|
| `l_low` | 3.0 | Lower threshold for outline opacity |
| `l_high` | 10.0 | Upper threshold — narrower range = jaggier outlines |
| `ikernel` | 1–4 | Kernel selection |
| `l_diff_min` | 0.0 Å | Minimum Z-height difference used in derivative |
| `l_diff_max` | 5.0 Å | Maximum Z-height difference — wider range detects larger features |

### B. Subunit/Chain Outlines

- Samples a 5×5 neighborhood of the object/chain ID buffer
- Counts how many neighbors belong to a different chain or biological assembly copy
- Thresholded with ramp (`r_low` to `r_high`, typically 3.0–10.0) to get opacity

### C. Residue Outlines

- Same 5×5 neighborhood sampling but on residue number
- Counts neighbors whose residue number differs by ≥ `resdiff` (typically 6000)
- Thresholded similarly (`g_low` to `g_high`, typically 3.0–8.0)

### Outline Combination

All outline layers combined via `max(contour_opacity, residue_opacity)`, with subunit outline
handled separately. Outlines darken the pixel:

```
final_color = (1 - outline_opacity) × shaded_color
```

### Blender Equivalent

Use the Compositor with:
- A Depth pass → Sobel/Laplacian filter for contour outlines
- Object Index or Material Index passes → edge detection for chain boundaries
- Combine with Mix nodes using the outline as a factor to multiply in black
- Alternatively, use Freestyle lines with different linesets for contour vs. chain boundaries

---

## 4. Depth Fog (Atmospheric Perspective)

Linear fog blending between the atom color and a fog color based on Z-depth:

```
fog_factor = front_fog - (z_max - z_pixel) / z_range × (front_fog - back_fog)
final = fog_factor × shadow × base_color + (1 - fog_factor) × fog_color
```

| Param | Typical | Effect |
|-------|---------|--------|
| `pfogh` | 1.0 | Fog factor at front (1.0 = no fog) |
| `pfogl` | 1.0 | Fog factor at back (0.0 = fully fogged) |
| `rfog(3)` | configurable | Fog RGB color |
| `rback(3)` | configurable | Background RGB color |

**Blender equivalent:** Use a Map Range node on the Z-depth pass in compositing, then Mix the render
with the fog color. Or use the Mist pass with appropriate start/end distances.

---

## 5. Final Pixel Composition Order

For each pixel, in order:

1. **Base color** — flat RGB from atom type lookup
2. **× Shadow factor** — conical soft shadow multiplier (0.2–1.0)
3. **Fog blend** — linear interpolation with fog color based on depth
4. **× (1 - outline_opacity)** — outlines darken toward black
5. **Alpha** — `max(is_atom ? 1 : 0, outline_opacity)` — outlines extend the alpha mask

A separate opacity image is written for compositing.

---

## 6. What's NOT Present

- No specular highlights
- No diffuse lighting / light direction
- No surface normals
- No anti-aliasing or supersampling
- No transparency/translucency (atoms are fully opaque)
- No reflections
- No subsurface scattering

---

## 7. Blender Recreation Strategy

1. **Material:** Emission shader or Diffuse BSDF with a flat color, no scene lights (or use
   Shader-to-RGB for full control)
2. **AO/Shadows:** Enable AO in World settings or use an AO node in the shader, with distance
   matching `coneangle`, intensity matching `pcone`, and a clamp node for `pshadowmax`
3. **Outlines:** Freestyle (contour + material boundary + edge marks for chains) or
   compositor-based edge detection on Depth + Object Index passes
4. **Fog:** Compositor Mist/Z-depth pass → Map Range → Mix with fog color
5. **Compositing order:** Render → AO multiply → Fog mix → Outline overlay → Output with alpha

---

## 8. Other Implementation Notes

- **Coordinate system:** Origin at upper-left, +X down, +Y left-to-right, +Z toward viewer
- **Clipping:** Molecules clipped at Z=0
- **Framebuffer:** Max 3000×3000 pixels
- **Output format:** PPM (Portable Pix Map) text format, 8-bit color
- **Biological units:** Supports BIOMT transformation matrices from PDB files
- **Transformation order:** Rotation → Centering → Translation
