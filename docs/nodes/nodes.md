---
title: "Nodes"
editor: visual
number-sections: true
number-depth: 2
toc-depth: 3
---

## Styling

### Atoms

Style to apply the traditional space-filling atomic representation of atoms. Spheres are scaled based on the `vdw_radii` attribute. By default the `Point Cloud` rendering system is used, which is only visible inside of Cycles.

![`MN_style_atoms`](videos/style_atoms.mp4)

#### Inputs:

| Name                | Data Type           | Default              | Description                                                                                                    |
|------------|------------|------------|--------------------------------------|
| Atoms               | Geometry            |                      | Vertices + edges describing atomic data to apply the style to.                                                 |
| Selection           | Boolean \<*Field*\> | `True`               | Boolean field, defining whether or not to apply the style.                                                     |
| Eevee / Cycles      | Boolean \<*Field*\> | `False`              | Boolean field, defining whether or not to instance geoemtry for Eevee or use Point Cloud rendering for Cycles. |
| Scale Radii         | Float \<*Field*\>   | `0.8`                | Scale the resulting spheres to a proportion of the atoms Van der Waal radii.                                   |
| Eevee: Subdivisions | Integer \<*Field\>* | `2`                  | Number of subdivisions for the instanced geometry, when using Eevee atoms.                                     |
| Eevee: Shade Smooth | Boolean \<*Field*\> | `True`               | Boolean field, whether to shade the Eevee atoms smooth.                                                        |
| Material            | Material            | `MN_atomic_material` | Material to apply to the resulting geometry.                                                                   |

### Cartoon

Style to apply the traditional cartoon representation of protein structures. This style highlights alpha-helices and beta-sheets with arrows and cylinders.

![](videos/style_cartoon.mp4)

| Name           | Data Type           | Default              | Description                                                                                                    |
|------------|------------|------------|--------------------------------------|
| Atoms          | Geometry            |                      | Vertices + edges describing atomic data to apply the style to.                                                 |
| Selection      | Boolean \<*Field*\> | `True`               | Boolean field, defining whether or not to apply the style.                                                     |
| Quality        | Integer \<*Field\>* | `2`                  | Integer field controlling the number of subdivisions and thus 'quality' of the cartoon.                        |
| Smooth / Sharp | Boolean \<*Field*\> | `False`              | Boolean field, whether or not to create smooth or sharp styler cartoon geometry.                               |
| Arrows         | Boolean \<*Field*\> | `True`               | Boolean field, whether or not to create arrows for beta-sheet cartoons.                                        |
| Cylinders      | Boolean \<*Field*\> | `False`              | Boolean field, whether or not to create cylinders for alpha-helix.                                             |
| Thickness      | Float \<*Field*\>   | `0.4`                | Float field controlling the thickness of the cartoons.                                                         |
| Width          | Float \<*Field*\>   | `1.8`                | Float field controlling the width of alpha-helices and beta-sheets, along the 'flat' part of the cartoon mesh. |
| Loop Radius    | Float \<*Field*\>   | `0.18`               | Float field controlling the radius of the loop geometry between secondary structure cartoons.                  |
| BS Smooth      | Float \<*Field*\>   | `1.00`               | Factor from `0` to `1`, controlling the amount to smooth the beta-sheet geometry that is created.              |
| Shade Smooth   | Boolean \<*Field*\> | `True`               | Boolean field, whether to shade the resulting cartoon mesh smooth.                                             |
| Material       | Material            | `MN_atomic_material` | Material to apply to the resulting geometry.                                                                   |

### Ribbon

Style that creates a continuous solid ribbon or licorice tube through the alpha-carbons of the structure.

![](videos/style_ribbon.mp4)

| Name           | Data Type           | Default              | Description                                                                                                    |
|------------|------------|------------|--------------------------------------|
| Atoms           | Geometry            |                      | Vertices + edges describing atomic data to apply the style to.                                                 |
| Selection       | Boolean \<*Field*\> | `True`               | Boolean field, defining whether or not to apply the style.                                                     |
| Quality         | Integer \<*Field\>* | `2`                  | Integer field controlling the number of subdivisions and thus 'quality' of the cartoon.                        |
| Radius          | Float \<*Field*\>   | `0.18`               | Float field controlling the radius of the loop geometry between secondary structure cartoons.                  |
| BS Smooth       | Float \<*Field*\>   | `1.00`               | Factor from `0` to `1`, controlling the amount to smooth the beta-sheet geometry that is created.              |
| Shade Smoothing | Boolean \<*Field*\> | `True`               | Boolean field, whether to shade the resulting cartoon mesh smooth.                                             |
| Material        | Material            | `MN_atomic_material` | Material to apply to the resulting geometry.                                                                   |

## Selection