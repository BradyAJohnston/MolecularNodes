---
title: Cryo-EM ensembles from starfiles
author: Johannes Elferich
bibliography: references.bib
---

::: {#fig-starfile-example-render}
![](https://i.imgur.com/Hoa9TRz.png){width=80%}

The protein beta-galactosidase as reconstructed by RELION-4.0 floating above the micrograph.
The image was rendered with Blender and the MolecularNodes plugin (Brady Johnston).
The data files are available via [RELION tutorial](https://relion.readthedocs.io/en/release-4.0/SPA_tutorial/Introduction.html) Copyright CC-BY 2024 Johannes Elferich.
:::

An essential step during cryo-EM data processing is to determine the location and orientation of molecules in cryo-EM micrographs.
This information is often saved in starfiles, which can be important to use this information.
By instancing structures or densities using the location and orientation information a mesoscale model of the cryo-EM sample can be generated and visualized.

## Example data

This tutorial uses single particle analysis tutorial data from Relion 4.0. 
Download instruction are [here](https://relion.readthedocs.io/en/release-4.0/SPA_tutorial/Introduction.html). 
You will only need the precalculated results.

## Starfile panel

In the Import panel of Molecular Nodes change the type to starfile. 

![](https://i.imgur.com/aj7wCXP.png)

Then you can select a starfile, optionally assign it a name, and click the `Load` button. For the puproses of this tutorial you can choose the data file of iteration 16 of the Refine3D job 29 `Refine3D/job029/run_it016_data.star`. The initial expected result is the creation of a new object in the MolecularNodes collection, which initially will show the position and rotation of particles as RGB-colored axes. You might want to switch the shading mode to `rendered` and disable the background prop in the default scene.

![](https://i.imgur.com/9THusDv.png)



## The starfile node

