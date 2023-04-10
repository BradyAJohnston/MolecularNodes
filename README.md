# MolecularNodes 🧬🍝💻

<a href="https://www.github.com/bradyajohnston/MolecularNodes/releases"><img src="https://img.shields.io/github/v/release/bradyajohnston/molecularnodes" style="height:30px" alt="Badge displaying license, which is MIT."></a>
<a href="https://www.github.com/bradyajohnston/MolecularNodes/releases"><img src="https://img.shields.io/github/downloads/BradyAJohnston/MolecularNodes/total.svg" style="height:30px" alt="Repo total downloads count."></a>
<a href="https://www.buymeacoffee.com/bradyajohnston"><img src="https://img.shields.io/github/license/bradyajohnston/molecularnodes" style="height:30px" alt="Badge displaying license, which is MIT."></a>
<a href="https://www.buymeacoffee.com/bradyajohnston"><img src="https://img.shields.io/github/stars/bradyajohnston/molecularnodes?style=social" style="height:30px" alt="Badge displaying count of GitHub stars."></a>

<a href="https://pypi.org/project/biotite"><img src="https://img.shields.io/badge/powered%20by-Biotite-orange.svg" style="height:30px" alt="Button linking to buymeacoffee.com to leave me a tip as a thank you."></a>
<a href="https://pypi.org/project/MDAnalysis"><img src="https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg" style="height:30px" alt="Button linking to buymeacoffee.com to leave me a tip as a thank you."></a>
<a href="https://pypi.org/project/mrcfile"><img src="https://img.shields.io/badge/powered%20by-mrcfile-orange.svg" style="height:30px" alt="Button linking to buymeacoffee.com to leave me a tip as a thank you."></a>

## About

MolecularNodes enables quick import and visualisation of structural biology data inside of Blender. Blender provides advanced industry-leading visualisation and animation technology, while MolecularNodes provides the interface that allows Blender to understand the unique data formats used in structural biology.


<img src="https://i.imgur.com/TNZLIGu.gif" style="height:300px">

MolecularNodes is a plugin for easy visualisation of structural biology and molecular data inside of the 3D modelling and animation program Blender. It provides a method for quick 1-button import of a range of formats including downloading directly from the protein data bank, or opening molecular dynamics trajectories.

MolecularNodes is powered by two key python packages:


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

## Thanks

<img src="https://download.blender.org/branding/blender_logo.png" style="height:80px" alt="Button linking to buymeacoffee.com to leave me a tip as a thank you.">
<a href="https://www.buymeacoffee.com/bradyajohnston"><img src="https://img.buymeacoffee.com/button-api/?text=Buy Me a Coffee&emoji=&slug=bradyajohnston&button_colour=eabae1&font_colour=000000&font_family=Poppins&outline_colour=000000&coffee_colour=FFDD00" style="height:80px" alt="Button linking to buymeacoffee.com to leave me a tip as a thank you."></a>


<img src="https://discord.com/api/guilds/940526858800336936/widget.png?style=banner1"
style="height:80px !important; width: 145px !important;"></a>