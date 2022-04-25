# Molecular Nodes

Molecular Nodes provides a convenient method for importing structural biology files into Blender, and several nodes for working with atomic data inside of Blender's Geometry Nodes.

Blender's Geometry Nodes provides a powerful interface for procedural modelling and animation. Currently it is limited in its ability to read any kind of structured data file as input, that isn't a 3D mesh. Molecular Nodes bridges this gap by providing an interface for converting `.pdb` and other file types into meshes that are usable by Geometry Nodes.

## Installation

To install Molecular Nodes, download the latest release and install the addon through the preferences panel inside of Blender, and select the `MolecularNodes_0.3.zip` file
> Edit -> Preferences -> Install

If the addon isn't highlighted, search for it and click the tick box to enable the addon.

### Installing Atomium
> Windows: You must run Blender as Administrator to be able to install Atomium successfully.

While still in the preferences panel, click the `Install Atomium` button to download and install the [Atomium python library](https://github.com/samirelanduk/atomium). 

The addon should now be enabled and available for use. 

## Loading Structures
You can load structures using the Molecular Nodes panel, which should be available at the top-right of the 3D viewport. You can press `N` to show & hide the side panels. Click on the `Molecular Nodes` tab to reveal the panel.

Molecular Noes provides two ways to open files. You can fetch directly from the PDB by inputting the 4-character code into the top box and pressing the download button. This will download the file and open it inside of Blender.

You can also open a file saved on your local computer. Click the folder icon to navigate to your `.pdb` file, select it and press OK. Then press `Open` on the Molecular Nodes panel to read the file and load the model into Blender.