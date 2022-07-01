# Molecular Nodes 🧬🔬💻


[![DOI](https://zenodo.org/badge/485261976.svg)](https://zenodo.org/badge/latestdoi/485261976) 


### Buy Me a Coffee to Keep Development Going!
<a href="https://www.buymeacoffee.com/bradyajohnston" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/v2/default-violet.png" alt="Buy Me A Coffee" style="height: 60px !important;width: 217px !important;" ></a>

### Join a Community of Blender SciVis People!
<a href="https://discord.gg/hbRp2NYDqK"><img src="https://discord.com/api/guilds/940526858800336936/widget.png?style=banner1"></a>


### What is Molecular Nodes?
Molecular Nodes provides a convenient method for importing structural biology files into Blender, and several nodes for working with atomic data inside of Blender's Geometry Nodes.

Blender's Geometry Nodes provides a powerful interface for procedural modelling and animation. Currently it is limited in its ability to read any kind of structured data file as input, that isn't a 3D mesh. Molecular Nodes bridges this gap by providing an interface for converting `.pdb` and other file types into meshes that are usable by Geometry Nodes.

![](img/atp-animation-demo.gif)

## Acknowledgements
This addon is built crucially around the python library [`atomium`](https://github.com/samirelanduk/atomium), which provides a great lightweight and pure-python way to to parse a variety of structural biology files, that interface nicely with the built-in python interpreter that is bundled with Blender.

The molecular dynamics trajectrory import is enabled by the [MDAnalysis package](https://github.com/MDAnalysis/mdanalysis).

[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)

This addon is built using the _excellent_ plugin for creating addons [Serpens](https://joshuaknauber.notion.site/Serpens-Documentation-d44c98df6af64d7c9a7925020af11233), which made the process of compiling a clear and concise addon, significantly reducing development time an creating a better addon overall.

Importantly none of this would be possible without the [Blender](https://www.blender.org/about/foundation/) and the [Blender Foundation](https://www.blender.org/about/foundation/) and the fantastic work of all the developers, from the community and the foundation itself. 

> ### This documentation is currently not complete for all of the nodes, but there should be enough to get you started. More complete documentation is in the works.

## Installation

Please see the [installation page](https://bradyajohnston.github.io/MolecularNodes/installation.html) of the documentation on how to install Molecular Nodes in Blender, and the [examples page](https://bradyajohnston.github.io/MolecularNodes/examples.html) on how to get started. 

You can also watch the [tutorial series on youtube](https://www.youtube.com/watch?v=CvmFaRVmZRU) on how to use Molecular Nodes.

[![image](https://user-images.githubusercontent.com/36021261/176899003-b8e03b44-8918-420b-94a2-a9fac379ca96.png)](https://www.youtube.com/watch?v=CvmFaRVmZRU)


