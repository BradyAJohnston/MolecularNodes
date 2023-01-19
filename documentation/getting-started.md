---
title: "Getting Started"
editor: visual
number-sections: true
number-depth: 1
---

This is a very basic introduction to Molecular Nodes.
How to import a protein from the PDB & change the colour and style via editing the node graph.
Molecular Nodes and Geometry Nodes in general has a lot more advanced functionality, and I encourage you to watch other youtube tutorials and spend time playing around with it to see what is possible.
Everything that I have achieved so far is through playing around to see what could be done.

::: callout-ip
## YouTube Tutorials

I have also made a series of YouTube tutorials walking through some of the functionality of Molecular Nodes.
Currently this tutorial series is for the older version of Molecular Nodes and not the 2.0 version.
The basic functionality and idea remain the same, but I will be updating this series once I get enough time.

<iframe width="560" height="315" src="https://www.youtube.com/embed/CvmFaRVmZRU" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen>

</iframe>
:::

# Downloading from the PDB

To download a model from the PDB, all you need is the 4-character PDB ID.

### Navigating to Molecular Nodes panel

Navigate to the Molecular Nodes panel under the Scene Properties tab.

![](images/CleanShot%202022-05-20%20at%2015.14.12.png)

The Molecular Nodes panel should now be visible in the bottom right of the screen, you may need to scroll down the Scene Properties tab.
If there is a button saying 'Install Packages' then you can click it to install the required python packages for your Blender installation.

![](https://i.imgur.com/LHM1fmk.png)

![](https://imgur.com/pWFq1FE)

### Downloading the PDB ID

Input the 4-character code for the structure you wish to download and press the 'Download' button to download and import the structure into blender.
You can select the default starting styling, if you want to start with 'Ribbon' or 'Ball and Stick' rather than the default 'Atoms Cycles'.

::: callout-caution
## Large Structures with 'Ball and Stick'

If importing large structures, I recommend starting with 'Atoms' or 'Ribbon'.
'Ball and Stick' is computationally quite a lot for Blender to process, and can crash your Blender on particularly large structures.
:::

![](images/CleanShot%202022-05-20%20at%2015.24.38@2x-01.png)

### Viewing the Structure

There should now be atoms visibly in the viewport, and a new collection in the collections tab that contains the structure and the properties for the imported protein.

![](images/CleanShot%202022-05-20%20at%2015.25.35.png)

By default, the atoms are only visible in the Cycles render engine.
Switch over to rendered view and change to Cycles to view the atoms of your structure.

![](images/CleanShot%202022-05-20%20at%2015.28.02.png)

# Geometry Nodes

Switch over to the Geometry Nodes tab and change the viewport to rendered view, to view the atomic structure again.

![](images/CleanShot%202022-05-20%20at%2015.32.52.png)

You can now interact with and manipulate the atoms in the structure using the nodes inside of geometry nodes.

Try scaling the radius of the atoms by adjust the `Scale Radii` value on the `MOL_atomic_properties` node, or changing the base colours of the different elements.

![](images/CleanShot%202022-05-20%20at%2015.35.16.gif)

## Atomic Properties

A limitation of Molecular Nodes currently is the inability to have text attributes.
To get around this, Molecular Nodes translates the different properties into numerical versions, and provides nodes to interact with the numeric versions of the data.

### Amino Acid Names

For example, the different amino-acid names have been assigned values from 1-20 based on alphabetical order.
The properties is available from the `MOL_atomic_properties` nodes under the `AA_name` property.
The property is a green diamond, indicating that it is an *integer field* and compatible with other numeric field sockets.

![](images/CleanShot%202022-05-20%20at%2015.42.31@2x.png)

As an example, we can generate a random colour based on the *numerical representation* of the amino acid name, and colour out atoms based on that.

Create a `Random Value` node by pressing `Shift + A` inside of the Geometry Nodes window and selecting *Utilities* -\> *Random Value*, and placing the node.
Change the random type from `Float` to `Vector`.
Plug the output of the `AA_name` into the `ID` of the `Random Value` node and plug the output into the different coloured sockets for the different elements in the `MOL_style_colour` node, as shown below.

We are creating a random numerical vector, with generated from the amino acid name.
While we have created a *numeric vector* and not a *colour vector,* Blender treats an `XYZ` vector the same as an `RGB` vector and *vice-versa*, so you can use them interchangeably.

![](images/CleanShot%202022-05-20%20at%2015.46.03@2x.png)

The different atoms should now be coloured based on their amino-acid name.

![](images/CleanShot%202022-05-20%20at%2015.47.03@2x.png)
