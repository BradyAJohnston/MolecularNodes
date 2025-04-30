import io
import warnings
from abc import ABCMeta
from pathlib import Path
from typing import Callable, List
import biotite.structure as struc
import bpy
import databpy
import numpy as np
from biotite import InvalidFileError
from biotite.structure import AtomArray, AtomArrayStack
from ... import blender as bl
from ... import download, utils
from ...nodes import nodes
from ...nodes.geometry import (
    GeometryNodeInterFace,
    add_style_branch,
    style_interfaces_from_tree,
)
from ..base import EntityType, MolecularEntity
from ..utilities import create_object
from . import pdb, pdbx, sdf, selections
from .reader import ReaderBase


class Molecule(MolecularEntity, metaclass=ABCMeta):
    """
    Primary Molecular Nodes class that coordinates the conversion of structural bioinformatic data
    into raw Blender data.  Most notable the conversion of atoms and bonds into a collection
    of vertices and lines.

    It associates the atomic data (the array) with the created 3D model inside of Blender
    (the object). If multiple conformations are imported, then a `frames` collection
    is also instantiated.

    The `named_attribute()` and `store_named_attribute()` methods access and set attributes on
    `object` that is in the Blender scene.

    Attributes
    ----------
    object : bpy.types.Object
        The Blender object representing the molecule.
    frames : bpy.types.Collection
        The Blender collection which holds the objects making up the frames to animate.
    array: AtomArray | AtomArrayStack:
        The numpy array which stores the atomic coordinates and associated attributes.
    select: MoleculeSelector
        A selector object that provides methods for creating atom selections based on various
        criteria such as atom name, residue type, chain ID, etc. These selections can be used
        with the `add_style` method to apply visual styles to specific parts of the molecule.
    """

    def __init__(
        self,
        array: AtomArray | AtomArrayStack,
        reader: ReaderBase | None = None,
    ):
        """
        Initialize the Molecule object.

        Parameters
        ----------
        array : AtomArray | AtomArrayStack
            The atom array or atom array stack containing the molecular data.
        reader : ReaderBase | None, optional
            The reader used to load the molecular data, by default None.
        """
        self._frames_collection: str | None = None
        self._code: str | None = None
        self._entity_type = EntityType.MOLECULE
        self._reader: ReaderBase | None = reader
        super().__init__()
        self.array = array
        self.select = MoleculeSelector(self)

    def create_object(self, name: str = "NewObject"):
        """
        Create a 3D model of the molecule, with one vertex for each atom.
        """
        self.object = create_object(
            array=self.array,
            name=name,
            collection=bl.coll.mn(),
        )
        if self._reader is not None:
            self._store_object_custom_properties(self.object, self._reader)
        self._setup_frames_collection()
        self._setup_modifiers()

    def _setup_modifiers(self):
        """
        Create the modifiers for the molecule.
        """
        self.object.modifiers.new("MolecularNodes", "NODES")
        tree = nodes.new_tree(  # type: ignore
            name=f"MN_{self.name}", input_name="Atoms", is_modifier=True
        )
        self.object.modifiers[0].node_group = tree  # type: ignore

    @classmethod
    def load(
        cls, file_path: str | Path, name: str | None = None, remove_solvent: bool = True
    ) -> "Molecule":
        """
        Load a molecule from a file.

        Parameters
        ----------
        file_path : str or Path
            The path to the file containing molecular data
        name : str or None, optional
            The name to give the molecule object. If None, uses the filename stem
        remove_solvent : bool, optional
            Whether to remove solvent molecules from the structure, default True

        Returns
        -------
        mol : Molecule
            The loaded molecule object with associated data and 3D representation

        Notes
        -----
        Supports various file formats including .cif, .bcif, .pdb, .sdf, and .mol
        """
        reader = cls._read(file_path)
        if not name:
            name = Path(file_path).stem
        mol = cls(reader.array, reader=reader)

        if remove_solvent:
            if isinstance(mol.array, AtomArrayStack):
                mol.array = mol.array[:, ~struc.filter_solvent(mol.array)]
            else:
                mol.array = mol.array[~struc.filter_solvent(mol.array)]

        mol.create_object(name=name)
        mol._reader = reader

        try:
            mol._assemblies = reader._assemblies()
        except InvalidFileError:
            mol._assemblies = ""

        return mol

    @property
    def code(self) -> str | None:
        """
        Get the code for the molecule.
        """
        return self._code

    @staticmethod
    def _read(
        file_path: str | Path | io.BytesIO,
        # remove_solvent: bool = False,
        # del_hydrogen: bool = False,
    ) -> ReaderBase:
        """
        Initially open the file, ready to extract the required data.

        Parameters
        ----------
        file_path : Union[Path, io.BytesIO]
            The file path to the file which stores the atomic coordinates.
        """
        if isinstance(file_path, io.BytesIO):
            reader = pdbx.PDBXReader(file_path)
        else:
            if isinstance(file_path, str):
                file_path = Path(file_path)

            match file_path.suffix:
                case ".cif":
                    reader = pdbx.PDBXReader(file_path)
                case ".bcif":
                    reader = pdbx.PDBXReader(file_path)
                case ".pdb":
                    reader = pdb.PDBReader(file_path)
                case ".sdf":
                    reader = sdf.SDFReader(file_path)
                case ".mol":
                    reader = sdf.SDFReader(file_path)
                case _:
                    raise InvalidFileError("The file format is not supported.")

        return reader

    @classmethod
    def fetch(
        cls,
        code: str,
        format=".bcif",
        centre: str | None = None,
        remove_solvent: bool = True,
        cache: Path | str | None = download.CACHE_DIR,
        database: str = "rcsb",
    ):
        """
        Fetch a molecule from the RCSB database.

        Parameters
        ----------
        code : str
            The PDB ID code of the molecule to fetch.
        format : str, optional
            The file format to download. Default is ".bcif".
        centre : str | None, optional
            Method to use for centering the molecule. Options are "centroid" (geometric center)
            or "mass" (center of mass). If None, no centering is performed. Default is None.
        cache : str, optional
            Path to cache directory. If None, no caching is performed.
        remove_solvent : bool, optional
            Whether to remove solvent from the molecule. Default is True.
        database : str, optional
            The database to fetch from. Default is "rcsb".

        Returns
        -------
        Molecule
            A new Molecule instance created from the downloaded data.
        """
        file_path = download.StructureDownloader(cache=cache).download(
            code=code, format=format, database=database
        )
        mol = cls.load(file_path, name=code, remove_solvent=remove_solvent)
        mol.object.mn["entity_type"] = "molecule"
        mol._code = code

        return mol

    def centre_molecule(self, method: str | None = "centroid"):
        """
        Offset positions to centre the atoms and vertices over either the geometric centroid
        or the centre of mass.
        """
        if method is None or method == "":
            return self

        adjustment = self.centroid(weight=method)
        self.position -= adjustment
        self.array.coord -= adjustment
        return self

    def centroid(self, weight: str | np.ndarray | None = None) -> np.ndarray:
        if weight == "centroid" or weight == "":
            return super().centroid()

        return super().centroid(weight)

    @property
    def tree(self) -> bpy.types.GeometryNodeTree:
        mod: bpy.types.NodesModifier = self.object.modifiers["MolecularNodes"]  # type: ignore
        if mod is None:
            raise ValueError(
                f"Unable to get MolecularNodes modifier for {self.object}, modifiers: {list(self.object.modifiers)}"
            )
        return mod.node_group  # type: ignore

    @property
    def frames(self) -> bpy.types.Collection | None:
        """
        Get the collection of frames for the molecule.

        Returns
        -------
        bpy.types.Collection
            The collection of frames for the molecule.
        """
        if self._frames_collection is None:
            return None

        return bpy.data.collections[self._frames_collection]

    @frames.setter
    def frames(self, value: bpy.types.Collection):
        """
        Set the collection of frames for the molecule.

        Parameters
        ----------
        value : bpy.types.Collection
            The collection of frames for the molecule.
        """
        if value is None:
            self._frames_collection = None
            return
        if not isinstance(value, bpy.types.Collection):
            raise TypeError("The frames must be a bpy.types.Collection.")

        self._frames_collection = value.name

    @property
    def n_models(self):
        """
        Get the number of models in the molecule.

        Returns
        -------
        int
            The number of models in the molecule.
        """
        if isinstance(self.array, struc.AtomArray):
            return 1
        else:
            return self.array.shape[0]

    def add_style(
        self,
        style: bpy.types.GeometryNodeTree | str = "spheres",
        color: str | None = "common",
        selection: "str | MoleculeSelector | None" = None,
        assembly: bool = False,
        material: bpy.types.Material | str | None = None,
    ):
        """
        Add a visual style to the molecule.

        Parameters
        ----------
        style : bpy.types.GeometryNodeTree | str, optional
            The style to apply to the molecule. Can be a GeometryNodeTree or a string
            identifying a predefined style (e.g., "spheres", "sticks", "ball_stick").
            Default is "spheres".

        color : str | None, optional
            The coloring scheme to apply. Can be "common" (element-based coloring),
            "chain", "residue", or other supported schemes. If None, no coloring
            is applied. Default is "common".

        selection : str | MoleculeSelector | None, optional
            Apply the style only to atoms matching this selection. Can be:
            - A string referring to an existing boolean attribute on the molecule
            - A MoleculeSelector object defining a selection criteria
            - None to apply to all atoms (default)

        assembly : bool, optional
            If True, set up the style to work with biological assemblies.
            Default is False.

        material : bpy.types.Material | str | None, optional
            The material to apply to the styled atoms. Can be a Blender Material object,
            a string with a material name, or None to use default materials. Default is None.

        Returns
        -------
        Molecule
            Returns self for method chaining.

        Notes
        -----
        If a MoleculeSelector is provided, it will be evaluated and stored as a new
        named attribute on the molecule with an automatically generated name (sel_N).
        """
        if style is None:
            return self

        if isinstance(selection, str) and selection not in self.list_attributes(
            drop_hidden=False
        ):
            warnings.warn(
                f"Named Attribute: '{selection}' does not exist. Style will be added but nothing will be displayed unless that attribute is created.",
                category=UserWarning,
            )

        if isinstance(selection, MoleculeSelector):
            name = "sel_0"
            i = 0
            while name in self.list_attributes():
                name = f"sel_{i}"
                i += 1

            self.store_named_attribute(
                selection.evaluate_on_array(self.array),
                name=name,
                atype=databpy.AttributeTypes.BOOLEAN,
                domain=databpy.AttributeDomains.POINT,
            )

            selection = name

        add_style_branch(
            tree=self.tree,
            style=style,
            color=color,
            selection=selection,
            material=material,
            frames=self.frames,
        )

        if assembly:
            nodes.assembly_initialise(self.object)
            nodes.assembly_insert(self.object)

        return self

    def _setup_frames_collection(self):
        if self.n_models > 1:
            self.frames = bl.coll.frames(self.name)
            for i, array in enumerate(self.array):
                databpy.create_object(
                    vertices=array.coord * self._world_scale,  # type: ignore
                    name="{}_frame_{}".format(self.name, str(i)),
                    collection=self.frames,
                )

    @property
    def styles(self) -> List[GeometryNodeInterFace]:
        """
        Get the styles in the tree.
        """
        return style_interfaces_from_tree(self.tree)

    @staticmethod
    def _store_object_custom_properties(obj, reader: ReaderBase):
        obj["entity_ids"] = reader.entity_ids()
        obj["chain_ids"] = reader.chain_ids()
        obj.mn.biological_assemblies = reader.assemblies(as_json_string=True)

    def assemblies(self, as_array=False):
        """
        Get the biological assemblies of the molecule.

        Parameters
        ----------
        as_array : bool, optional
            Whether to return the assemblies as an array of quaternions. Default is False.

        Returns
        -------
        dict or None
            The biological assemblies of the molecule, as a dictionary of
            transformation matrices, or None if no assemblies are available.
        """
        try:
            if self._reader is None:
                raise InvalidFileError
            assemblies_info = self._reader._assemblies()
        except InvalidFileError:
            return None

        if isinstance(assemblies_info, dict) and as_array:
            return utils.array_quaternions_from_dict(assemblies_info)

        return assemblies_info

    def __repr__(self) -> str:
        """
        Get the string representation of the Molecule object.

        Returns
        -------
        str
            The string representation of the Molecule object.
        """
        return f"<Molecule object: {self.name}>"


class MoleculeSelector:
    """
    A helper to create selections for Molecules and AtomArrays.

    The selection (self.mask) is not computed or returned until the `evaluate_on_array`
    method is called. Until then methods are stored for later evaluation.

    Parameters
    ----------
    mol : Molecule
        The molecule object to select from.

    Attributes
    ----------
    mol : Molecule
        The molecule object to select from.
    mask : ndarray or None
        Boolean array for the selection on the most recently evaluated array.
    pending_selections : list
        List of selection operations to be applied once `evaluate_on_array` is called.
    """

    def __init__(self, mol: Molecule | None = None):
        self.mol = mol
        self.mask: np.ndarray | None = None
        self.pending_selections: list[Callable] = []

    def _update_mask(self, operation):
        """
        Add a selection operation to the pending list.

        Parameters
        ----------
        operation : callable
            Selection operation to add

        Returns
        -------
        self : Selector
            Returns self for method chaining
        """
        self.pending_selections.append(operation)
        return self

    def reset(self):
        """
        Reset all pending selections and the mask

        Returns
        -------
        self : Selector
            Returns self for method chaining
        """
        self.pending_selections = []
        self.mask = None
        return self

    def store_selection(self, name: str) -> None:
        """
        Evaluate and store the current selection as a named attribute on the Molecule

        Parameters
        ----------
        name : str
            The name to store the selection under.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If no selection has been made.
        """
        if self.mol is None:
            raise ValueError("No Molecule associated with this selector")
        if self.mask is None or len(self.pending_selections) > 0:
            self.evaluate_on_array(self.mol.array)

        if self.mask is None:
            raise ValueError("No selection made.")

        self.mol.store_named_attribute(self.mask, name)
        return None

    def evaluate_on_array(self, array: AtomArray | AtomArrayStack) -> np.ndarray:
        """Evaluate this selection on the AtomArray.

        Parameters
        ----------
        array : AtomArray or AtomArrayStack
            The atomic structure to evaluate the selection on.

        Returns
        -------
        ndarray
            Boolean mask array indicating which atoms match the selection criteria.

        Notes
        -----
        All of the selection operations that have been staged for this Selector are
        evaluated and combined with a logical AND, using the AtomArray as input.
        """
        if isinstance(array, AtomArrayStack):
            array = array[0]  # type: ignore

        self.mask = np.ones(array.array_length(), dtype=bool)
        if not self.pending_selections:
            return self.mask

        for operation in self.pending_selections:
            self.mask &= operation(array)
        return self.mask  # type: ignore

    def atom_name(self, atom_name: str | list[str] | tuple[str, ...] | np.ndarray):
        """Select atoms by their name.

        Parameters
        ----------
        atom_name : str or list of str or tuple of str or ndarray
            The atom name(s) to select.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(
            lambda arr: selections.select_atom_names(arr, atom_name)
        )

    def is_amino_acid(self):
        """Select amino acid residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(selections.select_amino_acids)

    def is_canonical_amino_acid(self):
        """Select canonical amino acid residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(selections.select_canonical_amino_acids)

    def is_canonical_nucleotide(self):
        """Select canonical nucleotide residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(selections.select_canonical_nucleotides)

    def is_carbohydrate(self):
        """Select carbohydrate residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(selections.select_carbohydrates)

    def chain_id(self, chain_id: str | list[str] | tuple[str, ...] | np.ndarray):
        """Select atoms by chain identifier.

        Parameters
        ----------
        chain_id : list of str or tuple of str or ndarray
            The chain identifier(s) to select.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: selections.select_chain_id(arr, chain_id))

    def element(self, element: list[str] | tuple[str, ...] | np.ndarray):
        """Select atoms by element symbol.

        Parameters
        ----------
        element : list of str or tuple of str or ndarray
            The element symbol(s) to select.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: selections.select_element(arr, element))

    def is_hetero(self):
        """Select hetero atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(selections.select_hetero)

    def is_ligand(self):
        """Select ligand atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(selections.select_ligand)

    def linear_bond_continuity(self):
        """Select atoms with linear bond continuity.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(selections.select_linear_bond_continuity)

    def is_monoatomic_ion(self):
        """Select monoatomic ions.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(selections.select_monoatomic_ions)

    def is_nucleotide(self):
        """Select nucleotide residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(selections.select_nucleotides)

    def is_peptide_backbone(self):
        """Select peptide backbone atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(selections.select_peptide_backbone)

    def is_phosphate_backbone(self):
        """Select phosphate backbone atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(selections.select_phosphate_backbone)

    def is_backbone(self):
        """Select backbone atoms for peptide and nucleotide."""
        return self._update_mask(selections.select_backbone)

    def is_peptide(self):
        return self._update_mask(selections.select_peptide)

    def is_side_chain(self):
        return self._update_mask(selections.select_side_chain)

    def is_polymer(self):
        """Select polymer atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(selections.select_polymer)

    def res_id(self, num):
        """Select atoms by residue ID.

        Parameters
        ----------
        num : int or list of int
            The residue ID(s) to select.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: selections.select_res_id(arr, num))

    def res_name(self, res_name):
        """Select atoms by residue name.

        Parameters
        ----------
        res_name : str or list of str
            The residue name(s) to select.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: selections.select_res_name(arr, res_name))

    def is_solvent(self):
        """Select solvent atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(selections.select_solvent)

    def not_amino_acids(self):
        """Select non-amino acid residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~selections.select_amino_acids(arr))

    def not_atom_names(self, atomname):
        """Select atoms not matching the specified atom names.

        Parameters
        ----------
        atomname : str or list of str
            The atom name(s) to exclude.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(
            lambda arr: ~selections.select_atom_names(arr, atomname)
        )

    def not_canonical_amino_acids(self):
        """Select non-canonical amino acid residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(
            lambda arr: ~selections.select_canonical_amino_acids(arr)
        )

    def not_canonical_nucleotides(self):
        """Select non-canonical nucleotide residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(
            lambda arr: ~selections.select_canonical_nucleotides(arr)
        )

    def not_carbohydrates(self):
        """Select non-carbohydrate residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~selections.select_carbohydrates(arr))

    def not_chain_id(self, chain_id):
        """Select atoms not in the specified chains.

        Parameters
        ----------
        chain_id : str or list of str
            The chain identifier(s) to exclude.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~selections.select_chain_id(arr, chain_id))

    def not_element(self, element):
        """Select atoms not matching the specified elements.

        Parameters
        ----------
        element : str or list of str
            The element symbol(s) to exclude.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~selections.select_element(arr, element))

    def not_hetero(self):
        """Select non-hetero atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~selections.select_hetero(arr))

    def not_monoatomic_ions(self):
        """Select non-monoatomic ion atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~selections.select_monoatomic_ions(arr))

    def not_peptide(self):
        return self._update_mask(lambda arr: ~selections.select_peptide(arr))

    def not_side_chain(self):
        return self._update_mask(lambda arr: ~selections.select_side_chain(arr))

    def not_nucleotides(self):
        """Select non-nucleotide residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~selections.select_nucleotides(arr))

    def not_peptide_backbone(self):
        """Select non-peptide backbone atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~selections.select_peptide_backbone(arr))

    def not_phosphate_backbone(self):
        """Select non-phosphate backbone atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~selections.select_phosphate_backbone(arr))

    def not_polymer(self):
        """Select non-polymer atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~selections.select_polymer(arr))

    def not_res_id(self, num):
        """Select atoms not matching the specified residue IDs.

        Parameters
        ----------
        num : int or list of int
            The residue ID(s) to exclude.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~selections.select_res_id(arr, num))

    def not_res_name(self, res_name):
        """Select atoms not matching the specified residue names.

        Parameters
        ----------
        res_name : str or list of str
            The residue name(s) to exclude.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~selections.select_res_name(arr, res_name))

    def not_solvent(self):
        """Select non-solvent atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~selections.select_solvent(arr))
