# MolecularNodes 🧬🍝💻



<a href="https://www.buymeacoffee.com/bradyajohnston" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/v2/default-violet.png" alt="Buy Me A Coffee" style="height: 80px !important;width: 300px !important;" >    </a><a href="https://discord.gg/hbRp2NYDqK"><img src="https://discord.com/api/guilds/940526858800336936/widget.png?style=banner1"
style="height:80px !important; width: 145px !important;"></a>

MolecularNodes is a plugin for easy visualisation of structural biology and molecular data inside of the 3D modelling and animation program Blender. It provides a method for quick 1-button import of a range of formats including downloading directly from the protein data bank, or opening molecular dynamics trajectories.

MolecularNodes is powered by two key python packages:

<img src="https://www.biotite-python.org/_static/assets/general/biotite_logo_s.png" alt="Biotite Logo" style="height: 64px !important;width: 128px !important;" ></a>

  - [Biotite](https://www.biotite-python.org) which powers the majority of the file parsing and downloading from the PDB

[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)

  - [MDAnalysis](https://www.mdanalysis.org/) which powers the reading of a variety of molecular dynamics topology and trajectory files.

## Installation

See the [installation page](https://bradyajohnston.github.io/MolecularNodes/installation.html) of the documentation, for detailed instructions on how to install the addon.

## Getting Started
These tutorials are for earlier versions of the addon. There are some differences in design, but overall the workflow is the same. Watch through the videos to get an overview of how the addon works.

[![image](https://user-images.githubusercontent.com/36021261/205629018-a6722f88-505e-4cb6-a641-8d423aa26963.png)](https://youtu.be/CvmFaRVmZRU)

## Contributing
To contribute to the project, fork and clone the Molecular Nodes repo to your local machine. I recommend using VS Code and the [Blender VS Code](https://github.com/JacquesLucke/blender_vscode) addon which streamlines the development process. 

Once installed, you can use the `Blender: Build and Start` command with VS Code open in the addon directory, to start Blender with the addon built and installed. Any changes that are then made to the underlying addon code, can be quickly previewed inside of the running Blender by using the VS Code command `Blender: Reload Addonds`.

Once happy with your code, open a pull request to discuss and get it reviewed by others working on the project. Open a draft pull request early, or open an issue to discuss the scope and feasability of potential features.

## Citation

A paper has not yet been published on the addon, but if ou use it in your academic work you can site it from Zenodo:

[![DOI](https://zenodo.org/badge/485261976.svg)](https://zenodo.org/badge/latestdoi/485261976) 
