import io
import warnings
from abc import ABCMeta
from pathlib import Path
from typing import TYPE_CHECKING
import biotite.structure as struc
import bpy
import databpy
import numpy as np
from biotite import InvalidFileError
from biotite.structure import AtomArray, AtomArrayStack
from nodebpy.nodes.geometry import NamedAttribute
from ... import blender as bl
from ... import download, utils
from ...nodes import geometry as g
from ...nodes.nodes import STYLE_LITERALS, STYLE_NODE_MAPPING, styles_mapping
from ..base import EntityType, MolecularEntity
from ..utilities import create_object
from . import pdb, pdbx, sdf
from .annotations import MoleculeAnnotationManager
from .reader import ReaderBase

if TYPE_CHECKING:
    from ...ui.props import MolecularNodesObjectProperties


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
        self.annotations = MoleculeAnnotationManager(self)

    def create_object(self, name: str = "NewObject"):
        """
        Create a 3D model of the molecule, with one vertex for each atom.
        """
        self.object = create_object(
            array=self.array,
            name=name,
            collection=bl.coll.mn(),
        )
        self.props.entity_type = self._entity_type.value
        if self._reader is not None:
            self._store_object_custom_properties(self.object, self._reader)
        self._setup_frames_collection()
        self._setup_modifiers()

    @property
    def props(self) -> "MolecularNodesObjectProperties":
        """The custom Blender properties for the molecule object."""
        return self.object.mn  # ty: ignore[unresolved-attribute]

    def create_data_object(self) -> bpy.types.Object:
        from ... import utils
        from ...blender import mesh

        data_obj_name = f".data_{self.name}_assemblies"
        data_obj = bpy.data.objects.get(data_obj_name)
        if not data_obj:
            transforms = utils.array_quaternions_from_dict(
                self.props.biological_assemblies
            )
            data_obj = mesh.create_data_object(array=transforms, name=data_obj_name)

        return data_obj

    @classmethod
    def load(cls, file_path: str | Path, name: str | None = None) -> "Molecule":
        """
        Load a molecule from a file.

        Parameters
        ----------
        file_path : str or Path
            The path to the file containing molecular data
        name : str or None, optional
            The name to give the molecule object. If None, uses the filename stem

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
        cache: Path | str | None = download.CACHE_DIR,
        database: str = "rcsb",
    ) -> "Molecule":
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
        mol = cls.load(file_path, name=code)
        mol._code = code

        return mol

    def centroid(self, weight: str | np.ndarray | None = None) -> np.ndarray:
        if weight == "centroid" or weight == "":
            return super().centroid()

        return super().centroid(weight)

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
        style: STYLE_LITERALS = "spheres",
        selection: "str | None" = None,
        material: bpy.types.Material | None = None,
        **kwargs,
    ) -> "Molecule":
        """
        Add a visual style to the molecule.

        Provides a simple interface for adding visual styles to the molecule. For more complex
        styling, use the manual node tree creation via the `with mol.tree:` context manager.

        Parameters
        ----------
        style : str
            The style to apply to the molecule. Can be a string identifying a predefined
            style (e.g., "spheres", "sticks", "ball_stick"). Default is "spheres".

        selection : str | None, optional
            The selection to apply the style to. Can be the name of a boolean attribute
            on the molecule, or None to apply to all atoms.

        material : bpy.types.Material | None, optional
            The material to apply to the styled atoms. Can be a Blender Material object,
            a string with a material name, or None to use default materials. Default is None.

        **kwargs : optional
            Additional keyword arguments to pass to the added style node.

        Returns
        -------
        Molecule
            Returns self for method chaining.

        Raises
        ------
        ValueError
            If an unsupported style string is passed
        """
        if isinstance(style, str) and style not in styles_mapping:
            raise ValueError(
                f"Invalid style '{style}'. Supported styles are {[key for key in styles_mapping.keys()]}"
            )

        if isinstance(selection, str) and selection not in self.list_attributes(
            drop_hidden=False
        ):
            warnings.warn(
                f"Named Attribute: '{selection}' does not exist. Style will be added but nothing will be displayed unless that attribute is created.",
                category=UserWarning,
            )

        with self.tree as tree:
            (
                tree.atoms
                >> (
                    g.AnimateFrames(
                        frames=self.frames,
                        frame=g.AnimateValue(value_min=0, value_max=self.n_models),
                    )
                    if self.n_models > 1
                    else None
                )
                >> STYLE_NODE_MAPPING[style](
                    selection=NamedAttribute(selection) if selection else None,
                    material=material,
                    **kwargs,
                )
                >> tree.join
            )

        return self

    def _setup_frames_collection(self):
        if self.n_models > 1:
            self.frames = bl.coll.frames(self.name)
            for i, array in enumerate(self.array):
                databpy.create_object(
                    vertices=array.coord * self._world_scale,
                    name="{}_frame_{}".format(self.name, str(i)),
                    collection=self.frames,
                )

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

    def __getstate__(self):
        """Custom serialization."""
        state = super().__getstate__()
        # Remove objects with circular references or PyCapsules
        if "annotations" in state:
            del state["annotations"]
        return state

    def __setstate__(self, state):
        """Custom deserialization."""
        self.__dict__.update(state)
        # Recreate objects with circular references
        if not hasattr(self, "annotations"):
            self.annotations = MoleculeAnnotationManager(self)
