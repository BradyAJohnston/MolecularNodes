---
title: Molecular Dynamics
author: Brady Johnston
bibliography: references.bib
---

::: {#fig-md-example-render}
![](https://imgur.com/S42CNVJ.mp4){autoplay="true" loop="true"}

The protein adenylate kinase (AdK) undergoes a structural change during its catalytic cycle between a closed and an open state that is captured in enhanced sampling simulations [@seyler2015].
The protein secondary structure is shown as a round ribbon with individual amino acid sidechains as ball-and-sticks.
The trajectory was rendered with Blender and the MolecularNodes plugin (Brady Johnston).
Trajectory files are available via [MDAnalysisData](https://www.mdanalysis.org/MDAnalysisData/adk_transitions.html#adk-dims-transitions-ensemble-dataset) Copyright CC-BY 2023 Brady Johnston.
:::

As well as importing static structures, the results from molecular dynamics simulations can be imported as models in to Blender.
This is enabled through the excellent package [`MDAnalysis`](https://www.mdanalysis.org/).
The imported structure will have an object created that will act as the topology file.
Depending on the import method, the frames of the trajectory will either be streamed from the disk, or loaded in to memory inside of the `.blend` file as their own separate objects.

## MD Trajectory Panel

Use the `MD Trajectory` panel to import trajectories.

![](https://imgur.com/laYLRZu.png)

The minimum requirements are a valid topology file, that can be read by [MDAnalysis](https://userguide.mdanalysis.org/stable/formats/index.html).

If just a toplogy file is selected, then the model will be imported without any trajectory associated with it.
If a trajectory file is additionally chosen, then a trajectory will be associated with the toplogy file.

### Import Methods

There are two methods of importing the trajectory alongside the topology file. The default will stream the trajectory file from disk, while `In Memory` will load the entire selected trajectory in to memory inside of the `.blend` file. See [Importing](#importing-a-trajectory) for more info.

### Frames

If `In Memory` is selected, then you can choose which frames from the trajectory will be imported, with options for the first frame (indexed from 0), the last frame, and how many frames to skip (if any).

### Import Filter

When importing, you can filter the atoms that are imported, to potentially not import waters or other particular selections, by specifying the `Import Filter`.
This uses the [MDAnalysis selection language](https://userguide.mdanalysis.org/1.0.0/selections.html), and is only applied on import.

### Custom Selections

If you wish to still import atoms, but create a series of custom boolean selections for custom colouring or animation, then you can create custom selections in this panel.
Create a new selection, give it a name which will be used as the attribute name inside of Geometry Nodes, and the selection string which will be used to create the selection using the [MDAnalysis selection language](https://userguide.mdanalysis.org/1.0.0/selections.html).

## Importing a Trajectory

### Streaming

The default option will associate an `MDAnlaysis` session with the read topology file.
This will stream the topology from disk, as the frame in the scene inside of Blender changes.
If the original topology or trajectory files are moved, this will break the connection to the data.
This is the most performant option, but will potentially break if changing computers.

Below is an example of importing a trajectory, by streaming the frames.
As the frame changes in the scene, the loaded frame is updated on the imported protein, based on the created MDAnalysis session.
Interpolation between frames is currently not supported with this import method.

The MDAnalysis session will be saved when the `.blend` file is saved, and should be restored when the `.blend` file is reopened.

![](https://imgur.com/nACvzzd.mp4)

### In Memory

The `In Memory` option will load all frames of the trajectory in to memory, and store them as objects inside of the `MN_data` collection in the scene.
This will ensure that all of the associated data is stored inside of the `.blend` file for portability, but will come at the cost of performance for very large trajectories.
It also breaks the connection to the underlying `MDAnalysis` session, which limits the ability to further tweak the trajectory after import.

If `In Memory` is selected, the frames are imported as individual objects and stored in a `MN_data` collection.
The interpolation between frames is then handled by nodes inside of Geometry Nodes, which aren't necessarily linked to the scene frame.

This will create a larger `.blend` file and can lead to some performance drops with large trajectories, but ensures all of the data is kept within the saved file, and can enable further creative control through Geometry Nodes.

All connection to the underlying MDAnalysis session is lost on import, and the selections and trajectory cannot changed.
To make changes you must reimport the trajectory.

![](https://imgur.com/TK8eIaK.mp4)

## Creating the Animation

To replicate the animation which we see at the top of the tutorial, we can use some of the example datasets which are provided alongside `MDAnalysis` with the `MDAnalysisData` package.
To download one of the datasets, use the code below:

``` python
# pip install MDAnalysisData
from MDAnalysisData import datasets
datasets.fetch_adk_transitions_FRODA()
```

### Loading the Trajectory

We will load the trajectory, and load all of the frame sin to memory to ensure we can make a smoother trajectory.

In the video below we have imported the trajectory, and we can adjust the number of frames in the scene, as well as the number of frames the trajectory will play back over.
We also enabled `EEVEE` atoms to display in the EEVEE render engine.

![](https://imgur.com/jKTYWp9.mp4)

#### Changing Styles

We can change the style of the imported trajectory, by adding a new style node. We can combine styles with the `Join Geometry`. For more details on adding styles, see the (importing)[01_importing.qmd] tutorial.

![](https://imgur.com/nhau0r9.mp4)

We can apply the atoms style, only to the side chains of the protein, by using the `Backbone` selection node, and using the `is_side_chain` output. This selectively applies the style to only those atoms in the selection. The combined styles now contain only the atoms for the side chains and a continuous ribbon for the protein.

![](https://imgur.com/1m3pHKM.mp4)


### Setting the Scene

We can set up the scene a bit nicer with a backdrop. In this case we create a plane using <kbd>Shift</kbd> + <Kbd>A</kbd> to add a plane, go in to [edit mode](#01-introduction-edit-mode) and extrude the backbdrop up with the <kbd>E</kbd> key. We can create a slightly curved corner by bevelling the corner. Select the two vertices of the edge and click <kbd>Ctrl</kbd> + <kbd>B</kbd>. Move the mouse and use the scroll wheel to adjust the settings, then left click to apply.

![](https://imgur.com/6LUQEnz.mp4)

### Rendering the Animation

We can change some final settings of the style, do a test `Render Image`, change the export settings for where the frames of the animation are going to be saved, then we can click `Render Animation` to render all of the frames of the animation.

![](https://imgur.com/IBKUQSr.mp4)