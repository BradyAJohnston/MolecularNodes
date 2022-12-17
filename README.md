# MolecularNodes 🧬🍝💻
Molecular Nodes has undergone a complete re-write from the ground-up, to fix a number of underlying bugs and brings a large number of speed improvements, but also removes a large amount of technical debt that enables easier collaboration and contribution from others wishing to work on the project.

The documenation for the [current release](https://bradyajohnston.github.io/MolecularNodes/) of Molecular Nodes is still available, while this refactored version is being worked on and improved.


## Installation

To install Molecular Nodes, DO NOT download this repo. Instead, navigate to the [release page](https://github.com/BradyAJohnston/MolecularNodes/releases) and download the most recent release. The current version of Molecular Nodes (2.0.X) requires Blender 3.4.

> The file needs to remain zipped as `MolecularNodes_2.0.X.zip`. If your browser automaticlaly unzips the file (Safari on MacOS) then re-zip the file or download with a different browser (Chrome or Firefox) and then install the `.zip` file.

Inside of Blender, navigate to `Edit` -> `Preferences` and click `Addons` then `Install` and select the `MolecularNodes_2.0.X.zip` file. Click 'Enable' which is the tick box in the top left corner. 

The MolecularNodes panel can be found under Scene Properties, where there should be a button to click to install the requried python packages.

### MolecularNodes `v0.13` Installation

The installation instructions for the previous version [can be found here](https://bradyajohnston.github.io/MolecularNodes/installation.html). They are similar but differ slightly to the curren installation instructions.

[![image](https://user-images.githubusercontent.com/36021261/205629018-a6722f88-505e-4cb6-a641-8d423aa26963.png)](https://youtu.be/CvmFaRVmZRU)

## Contributing
To contribute to the project, fork and clone the Molecular Nodes repo to your local machine. I recommend using VS Code and the [Blender VS Code](https://github.com/JacquesLucke/blender_vscode) addon which streamlines the development process. 

Once installed, you can use the `Blender: Build and Start` command with VS Code open in the addon directory, to start Blender with the addon built and installed. Any changes that are then made to the underlying addon code, can be quickly previewed inside of the running Blender by using the VS Code command `Blender: Reload Addonds`.

Once happy with your code, open a pull request to discuss and get it reviewed by others working on the project. Open a draft pull request early, or open an issue to discuss the scope and feasability of potential features.
