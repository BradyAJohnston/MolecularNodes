---
title: Cryo-EM particles
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

Then you can select a starfile, optionally assign it a name, and click the `Load` button. 
For the puproses of this tutorial you can choose the data file of iteration 16 of the Refine3D job 29 `Refine3D/job029/run_it016_data.star`. 
The initial expected result is the creation of a new object in the MolecularNodes collection, which initially will show the position and rotation of particles as RGB-colored axes. 
You will want to switch the shading mode to `rendered` and disable the background prop in the default scene.

![](https://i.imgur.com/9THusDv.png)



## The starfile node

To get a better view of the instances, press the `Z` in the axis gizmo of the 3D viewport. 
Further control of the instance display is possibly in the `Starfile Instances` node. 
For example, enabling `Show Micrograph` will display the micrograph the instances were derived from.
You will want to adjust the Brightness/Constrast setting under `Micrograph options`.
For this micrograph 0.2 and 0.3 work reasonably well.
Another important control is the `Image` input, which will allow you to "scroll" through the micrographs of this dataset.
And, of course, it would be way more visually appealing to see the actual structure of beta-gal instead of the axes.
For this purpose switch the `Method` of the `Import` panel to `Density` and select a mrc file. 
In this case I have chosen the masked density of PostProcessing Job 30.
Very importantly, **enable the `Center Density` option.** 
Otherwise the maps will not line up correctly with the micrograph.

![](https://i.imgur.com/c3OVFwz.png)

After you click `Load`, a new object will be created.
Select it in the outliner and in the 3D vieport `View` menu click `Frame selected` to inspect it more closely.
Adjust the `Threshold` and `Hide dust` inputs of the `Style density` node to desired values and at this point decide on a good material and color.
Note how the map is precicely centered around the origin of the blender coordinate system, which will allow for accurate placement on the micrograph.

![](https://i.imgur.com/ILA37AV.png)

After you have adjusted the display of the map to your lining, reselect the starfile object and again use `Frame selected` to readjust the display.
Now you will be able to choose the `bGal map` in teh `Molecule` input of the `Starfile Instances` node. 
Disable the display of the `bGal map` in the outliner for both the viewpoer and render to only show the starfile instances.


![](https://i.imgur.com/LiJtqdD.png)