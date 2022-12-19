# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added 

- `CHANGELOG.md` for tracking changes to the project
- Issue templates for GitHub issues
- Proper support for Fields for ribbon width, enabling 'licorice' representation amoung others
- Tooltips documentation for each of the custom node groups that can be added in Geometry Nodes
- More selection nodes for distance, XYZ slice and whole residues.

### Fixed

- Error when defaulting to `connect_via_distance()` when importing with 'Find Bonds' enabled
- Adding of a color node which was mis-labelled and couldn't be added

## [2.0.2] - 2022-12-13

### Fixed

- Error on reporting the success of improting a molecule

## [2.0.1] - 2022-12-13

### Changed
- Remove usage of [Atomium](https://github.com/samirelanduk/atomium) and switched to [Biotite](https://github.com/biotite-dev/biotite) for most internal structural file parsing
- Removed reliance upon [Serpens](https://blendermarket.com/products/serpens) to build the addon. Can now be developed by anybody without the requirements of additional Blender addons.
  - Aspects of the addon can now have changes tracked properly by git, rather than inside a binary `.blend` Serpens file.
  - Now also allows easier and clearer contributions from others who wish to contribute
  - Custom node trees still remain opaque to git unfortunately

- Attributes are now created and stored on the molecular model itself, removing the need for properties models that store the data XYZ positions. Makes the import process clearer, more robust and more easily expandable.
- Bond information is available by default as edges of the mesh along with bond types

- Improve MDAnalysis installation and usage internally
- Expose operators and functions to Blender Python console to enable user scripting with the addon

```python
import MolecularNodes as mn

# to fetch structures from the protein data bank
struc_list = ['4ozs', '1xi4', '6n2y']
for pdb in struc_list:
    mn.load.molecule_rcsb(pdb, starting_style = 1)

# to open a local structure file
mn.load.molecule_local('file_path_here.pdb', 'SomeMoleculeName', starting_style = 2)
```
