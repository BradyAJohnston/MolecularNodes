# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [UNRELEASED]

### Added
- Custom text input for creating selections of `res_id` ranges and indivual numbers thanks to @YaoYinYing ([#149](https://github.com/BradyAJohnston/MolecularNodes/pull/149)), enabling quicker creation of complex selections inside of Molecular Nodes.

### Fixed
- Refactor of package installation via pip, to help with those who require `pip` mirrors and provide more information when installation fails on ARM macs. ([#162](https://github.com/BradyAJohnston/MolecularNodes/pull/162))
- Removed redundant python submodules and general cleanup

## [2.3.1](https://github.com/BradyAJohnston/MolecularNodes/releases/tag/v2.3.1) - 2023-01-26

### Added

### Fixed

- Change translations for other languages to `bpy.app.translations.pgettext_data()` from `bpy.app.translations.pgettext()` which reduces potential conflicts with other addons. ([#147](https://github.com/BradyAJohnston/MolecularNodes/pull/147))
- Fixed `chain_id_unique` not being added when importing via MDAnalysis, stopping custom nodes relying on chain information to break. ([#156](https://github.com/BradyAJohnston/MolecularNodes/pull/156))

## [[2.3.0]](https://github.com/BradyAJohnston/MolecularNodes/releases/tag/v2.3.0) - 2023-01-26

### Added
- Panel for adding multiple selection strings, which will become boolean attributes on the imported model when importing via MDAnalysis.

### Fixed
- Ball and stick node sphers now support field input for scaling the radius
- Error with initial node setup breaking when in non-english Blender UI ([#139](https://github.com/BradyAJohnston/MolecularNodes/pull/139)) contributed by @YaoYinYing
- Problems with biological assemblies failling on larger structures. ([#143](https://github.com/BradyAJohnston/MolecularNodes/pull/145))
- Problem with Animate Frames node defaulting to wrong from range on start

## [[2.2.2]](https://github.com/BradyAJohnston/MolecularNodes/releases/tag/v2.2.2) - 2023-01-06

### Added

### Fixed
- Issue on linux and with newer versions of Numpy where `np.bool` is deprecated and was erroring on import.

## [[2.2.1]](https://github.com/BradyAJohnston/MolecularNodes/releases/tag/v2.2.1) - 2023-01-05

### Added
- multi-model `b_factor` is added when importing from `.pdb` files via biotite [#133](https://github.com/BradyAJohnston/MolecularNodes/pull/133)
- 'Invert' field option to atom_properties and other selection nodes to optionally invert the selection
- Added better detection of ligands and modifcations (such as sugars) and a separate selection node for them. Currently ligands are stored on the `res_name` attribute, starting at 100 and incrementing one for each unique ligand.

### Fixed
- `include_bonds` option was not being utilised on MD import [#132](https://github.com/BradyAJohnston/MolecularNodes/pull/132)
- `MOL_animate_res_wiggle` was wiggling the `OXT` (`res_name == 38`) oxygen when a peptide chain ended. Added additional selection to not wiggle this atom, which should only ever appear when a peptide chain terminates.
- fixed import of `vdw_radii` for elements not supported by biotite (such as Fe) by moving vdw_radii to the data dictionary rather than relying on a function from biotite which had a limited dictionary for vdw_radii lookup


## [[2.2.0]](https://github.com/BradyAJohnston/MolecularNodes/releases/tag/v2.2.0) - 2023-01-03

### Added
- `atom_name` attribute, which is a numerical representation of the atom name (C, CA, C5 etc)
  - Dicussed in [#118](https://github.com/BradyAJohnston/MolecularNodes/issues/118)
  - Allows for more precise selections for new styling and animation nodes
- Reimplemented amino acid 'wiggle' node: using the `atom_name` attribute
  - 3x faster with the improved `atom_name` attribute and refactor of the underlying nodes
- Reimplemented the amino acid to curve node
- Ribbon Style Nucleic: A ribbon representation for nucleic acids to complement the ribbon represenation of the proteins from alpha carbons. Enabled now thanks to the `atom_name` attribute.
- Capturing the `Index` field in the selection node before the selection occurs, and added an `Index` field input to the `MOL_animate_frames` node to enable selection to occur before animating between frames, if the `Pre-Sel Index` field is used in the `Index` field of the `MOL_animate_frames` node
- Added cutoff field for limiting the interpolation of atoms between frames based on a distance cutoff
- Added bonds through MDAnalysis import when a trajectory supports it [#129](https://github.com/BradyAJohnston/MolecularNodes/issues/129)
- Added `is_solvent`, `is_nucleic` and `is_peptide` boolean attributes when importing via MDAnalysis
- Added frame-specific attribute `occupancy` which is added to each frame of the trajectory when imported via MDAnalysis. [#128](https://github.com/BradyAJohnston/MolecularNodes/issues/128)

### Fixed
- Changed naming of `MOL_style_atoms` to `MOL_style_atoms_cycles` and `MOL_style_ribbon` to `MOL_style_ribbon_protein`

## [[2.1.0]](https://github.com/BradyAJohnston/MolecularNodes/releases/tag/v2.1.0) - 2023-01-01

### Added 

- `CHANGELOG.md` for tracking changes to the project
- Issue templates for GitHub issues
- Proper support for Fields for ribbon width, enabling 'licorice' representation amoung others
- Tooltips documentation for each of the custom node groups that can be added in Geometry Nodes
- More selection nodes for distance, XYZ slice and whole residues.
- Custom selections using a string when importing via MDAnalysis [#123](https://github.com/BradyAJohnston/MolecularNodes/pull/123)
- Added option to disable interpolation of atom positions between frames
- Added node for animating any field between frames of a trajectory (no fields currently added on import, but used in the new Aniamte Frames backend)
- UVs are now available for the ribbon mesh style, idea from [ErikMarklund](https://github.com/BradyAJohnston/MolecularNodes/issues/77#issuecomment-1273598042) and implemented by [quellenform](https://github.com/quellenform/blender-CurveToMeshUV)

### Fixed

- Error when defaulting to `connect_via_distance()` when importing with 'Find Bonds' enabled
- Adding of a color node which was mis-labelled and couldn't be added
- Non-`.gro` topology files were failling to add `vdw_radii` attribute [#124](https://github.com/BradyAJohnston/MolecularNodes/issues/124)
- Remove use of `np.int` which is now deprecated and was causing errors when [linking python on M1 Mac](https://github.com/BradyAJohnston/MolecularNodes/issues/108#issuecomment-1365540371)
- Attributes now available on ribbon mesh which are [sampled from backbone](https://github.com/BradyAJohnston/MolecularNodes/issues/77)
- Changed starting material to be appended instead of created, which should avoid duplication of material.

## [[2.0.2]](https://github.com/BradyAJohnston/MolecularNodes/releases/tag/v2.0.2) - 2022-12-13

### Fixed

- Error on reporting the success of improting a molecule

## [[2.0.1]](https://github.com/BradyAJohnston/MolecularNodes/releases/tag/v2.0.1) - 2022-12-13

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
