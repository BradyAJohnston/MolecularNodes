# Molecular Nodes 🧬🍝💻

<img src="docs/images/logo.png" align="right" style = "height:250px;"/>


![Documentation Building](https://github.com/bradyajohnston/molecularnodes/actions/workflows/docs.yml/badge.svg) ![Running Tests](https://github.com/bradyajohnston/molecularnodes/actions/workflows/tests.yml/badge.svg) [![codecov](https://codecov.io/gh/BradyAJohnston/MolecularNodes/branch/main/graph/badge.svg?token=ZB2SJFY8FU)](https://codecov.io/gh/BradyAJohnston/MolecularNodes)


<a href="https://www.github.com/bradyajohnston/MolecularNodes/releases"><img src="https://img.shields.io/github/v/release/bradyajohnston/molecularnodes" alt="Badge displaying license, which is MIT." style="height:20px"/></a> <a href="https://www.github.com/bradyajohnston/MolecularNodes/releases"><img src="https://img.shields.io/github/downloads/BradyAJohnston/MolecularNodes/total.svg" alt="Repo total downloads count." style="height:20px"/></a> <a href="https://www.buymeacoffee.com/bradyajohnston"><img src="https://img.shields.io/github/license/bradyajohnston/molecularnodes" alt="Badge displaying license, which is MIT." style="height:20px"/></a> <a href="https://www.buymeacoffee.com/bradyajohnston"><img src="https://img.shields.io/github/stars/bradyajohnston/molecularnodes?style=social" alt="Badge displaying count of GitHub stars." style="height:20px"/></a>

<a href="https://pypi.org/project/biotite"><img src="https://img.shields.io/badge/powered%20by-Biotite-orange.svg" alt="Badge showing usage of MDAnalysis as a python package powering the add-on" style="height:20px"/></a> <a href="https://pypi.org/project/MDAnalysis"><img src="https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg" alt="Badge showing usage of biotite as a python package powering the add-on" style="height:20px"/></a>

 <a href="https://patreon.com/bradyajohnston"><img src="https://img.shields.io/endpoint.svg?url=https%3A%2F%2Fshieldsio-patreon.vercel.app%2Fapi%3Fusername%3Dbradyajohnston%26type%3Dpatrons&style=for-the-badge" alt="Support me on Patreon with ongoing donations" style="height:35px"/></a>

 <a href="https://buymeacoffee.com/bradyajohnston"><img src="https://img.shields.io/badge/Buy%20Me%20a%20Coffee-ffdd00?style=for-the-badge&logo=buy-me-a-coffee&logoColor=black" alt="Support me by buying a couple of coffees as a one-off donation" style="height:35px"/></a>


 <a href="https://discord.gg/fvw6vT3vY9"><img src="https://img.shields.io/badge/Discord-%235865F2.svg?style=for-the-badge&logo=discord&logoColor=white" alt="Join the scientific visualisation in Blender discord." style="height:35px"/></a>

## About

MolecularNodes enables quick import and visualisation of structural biology data inside of Blender. Blender provides advanced industry-leading visualisation and animation technology, while MolecularNodes provides the interface that allows Blender to understand the unique data formats used in structural biology.

The add-on enables creating animations from static crystal structures, styling proteins and other molecules in a variety of highly customisable styles, importing and playing back molecular dynamics trajectories from a wide variety of sources, and even importing of EM density maps.

## Examples

See examples, tutorials and video projects that use Molecular Nodes in the [documentation page](https://bradyajohnston.github.io/MolecularNodes/examples/), or academic papers that use Molecular Nodes in the [citations page](https://bradyajohnston.github.io/MolecularNodes/citations/).

| Clockwork | Veritasium | BCON22 Talk |
| --- | --- | --- |
| ![[Clockwork](https://www.youtube.com/watch?v=lv89fSt5jBY)](https://img.youtube.com/vi/lv89fSt5jBY/0.jpg) | ![[Veritasium](https://www.youtube.com/watch?v=P_fHJIYENdI)](https://img.youtube.com/vi/P_fHJIYENdI/0.jpg) | ![[Clockwork](https://www.youtube.com/watch?v=adhTmwYwOiA)](https://img.youtube.com/vi/adhTmwYwOiA/0.jpg) |

## Installation

Molecular Nodes can be installed from within Blender in versions >=4.2, by using the _Get Extensions_ menu. More details can be found on the [installation page](https://bradyajohnston.github.io/MolecularNodes/tutorials/installation.html) of the documentation.

## Getting Started

There are video and written tutorials in the [documentation](https://bradyajohnston.github.io/MolecularNodes/tutorials/) that will walk you through the basics of using the addon.

Documentation for all of the individual nodes are also available in the [node documentation](https://bradyajohnston.github.io/MolecularNodes/nodes/style.html) page.

[![image](https://user-images.githubusercontent.com/36021261/205629018-a6722f88-505e-4cb6-a641-8d423aa26963.png)](https://youtu.be/CvmFaRVmZRU)

## Contributing

If you would like to contribute to the project, please open an issue to discuss potential new features, or comment on an existing issue if you would like to help with fixing it. I welcome any and all potential PRs.

It's recommended to clone this repository using `git clone --depth 1` as the complete commit history gets close to 1GB of data. I also recommend using VS Code with the [Blender VS Code](https://github.com/JacquesLucke/blender_vscode) addon which streamlines the development process.

Once installed, you can use the `Blender: Build and Start` command with VS Code open in the addon directory, to start Blender with the addon built and installed. Any changes that are then made to the underlying addon code, can be quickly previewed inside of the running Blender by using the VS Code command `Blender: Reload Addonds`.

Once happy with your code, open a pull request to discuss and get it reviewed by others working on the project. Open a draft pull request early, or open an issue to discuss the scope and feasability of potential features.

## Citation

A paper has not yet been published on the addon, but if you use it in your academic work you can site it from Zenodo:

[![](https://zenodo.org/badge/485261976.svg)](https://zenodo.org/badge/latestdoi/485261976)

## Thanks

A massive thanks to the [Blender Foundation](https://blender.org) which develops Blender as a free and open source program, and to the python package developers who enable the functionality of the this add-on. Primaryil Biotite and MDAnalysis teams.

<img src="https://download.blender.org/branding/blender_logo.png" alt="The Blender logo." style="height:80px"/>
