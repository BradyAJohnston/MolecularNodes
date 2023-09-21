---
title: "Molecular Dynamics"
---

## Importing a Trajectory

Use the `MD Trajectory` panel to import trajectories.

![](https://imgur.com/laYLRZu.png)

The minimum requirements are a valid topology file, that can be read by [MDAnalysis](https://userguide.mdanalysis.org/stable/formats/index.html).

If just a toplogy file is selected, then the model will be imported without any trajectory associated with it. If a trajectory file is additionally chosen, then a trajectory will be associated with the toplogy file. 

### Import Methods

There are two methods of importing the trajectory alongside the topology file.

#### Streaming the Trajectory
The default option will associate an `MDAnlaysis` session with the read topology file. This will stream the topology from disk, as the frame in the scene inside of Blender changes. If the original topology or trajectory files are moved, this will break the connection to the data. This is the most performant option, but will potentially break if changing computers.

#### In Memory Trajectory
The `In Memory` option will load all frames of the trajectory in to memory, and store them as objects inside of the `MN_data` collection in the scene. This will ensure that all of the associated data is stored inside of the `.blend` file for portability, but will come at the cost of performance for very large trajectories. It also breaks the connection to the underlying `MDAnalysis` session, which limits the ability to further tweak the trajectory after import.

### Frames
If `In Memory` is selected, then you can choose which frames from the trajectory will be imported, with options for the first frame (indexed from 0), the last frame, and how many frames to skip (if any).

### Import Filter

When importing, you can filter the atoms that are imported, to potentially not import waters or other particular selections, by specifying the `Import Filter`. This uses the [MDAnalysis selection language](https://userguide.mdanalysis.org/1.0.0/selections.html), and is only applied on import.

### Custom Selections

If you wish to still import atoms, but create a series of custom boolean selections for custom colouring or animation, then you can create custom selections in this panel. Create a new selection, give it a name which will be used as the attribute name inside of Geometry Nodes, and the selection string which will be used to create the selection using the [MDAnalysis selection language](https://userguide.mdanalysis.org/1.0.0/selections.html).


## Streaming Trajectory

Below is an example of importing a trajectory, by streaming the frames. As the frame changes in the scene, the loaded frame is updated on the imported protein, based on the created MDAnalysis session. Interpolation between frames is currently not supported with this import method.

![](https://imgur.com/nACvzzd.mp4)

## In Memory

If `In Memory` is selected, the frames are imported as individual objects and stored in a `MN_data` collection. The interpolation between frames is then handled by nodes inside of Geometry Nodes, which aren't necessarily linked to the scene frame.

This will create a larger `.blend` file and can lead to some performance drops with large trajectories, but ensures all of the data is kept within the saved file, and can enabled further creative control through Geometry Nodes. 

All connection to the underlying MDAnalysis session is lost on import, and the selections and trajectory cannot be updated but have to be reimported if you wish to make changes.

![](https://imgur.com/TK8eIaK.mp4)