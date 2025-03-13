import io
from abc import ABCMeta
from pathlib import Path

import biotite.structure as struc
import bpy
import databpy
import numpy as np
from biotite import InvalidFileError
from biotite.structure import AtomArray, AtomArrayStack
from databpy import AttributeDomains, AttributeTypes

from ... import blender as bl
from ... import download, utils
from ..base import EntityType, MolecularEntity
from ..utilities import create_object, atom_array_to_named_attributes
from . import pdb, pdbx, sdf
from .reader import ReaderBase


class Molecule(MolecularEntity, metaclass=ABCMeta):
    """
    Primary Molecular node that coordinates the conversion of structureal bioinfromatic data
    into raw Blender data.  Most notable the conversion of atoms and bonds into a collection
    of vertices and lines.

    It associates the atomic data (the array) with the created 3D model inside of Blender
    (the object). If multiple conformations are imported, then a `frames` collection
    is also instantiated.

    The `named_attribute()` and `store_named_attribute()` methods access and set attributes on
    `object` that is in the Blender scene.

    Attributes
    ----------
    file_path : str
        The file path to the file which stores the atomic coordinates.
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
        file_path : Union[str, Path, io.BytesIO]
            The file path to the file which stores the atomic coordinates.
        """
        self._frames_collection: str | None = None
        self._code: str | None = None
        self._entity_type = EntityType.MOLECULE
        self._reader: ReaderBase | None = reader
        super().__init__()
        self.array = array
        # self.object: bpy.types.Object | None = None


    def create_object(self, name: str = "NewObject"):
        """
        Create a 3D model of the molecule, with one vertex for each atom.
        """
        self.object = create_object(
            array=self.array, name=name, collection=bl.coll.mn(),
        )
        if self._reader is not None:
            self._store_object_custom_properties(self.object, self._reader)
        self._setup_frames_collection()


    @classmethod
    def load(cls, file_path: str | Path, name: str | None = None):
        reader = cls._read(file_path)
        if not name:
            name = Path(file_path).stem
        mol = cls(reader.array, reader=reader)

        # currently filtering out solvent, will make optional in another PR
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
        # del_solvent: bool = False,
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
        cache : str, optional
            Path to cache directory. If None, no caching is performed.
        database : str, optional
            The database to fetch from. Default is "rcsb".

        Returns
        -------
        Molecule
            A new Molecule instance created from the downloaded data.
        """
        file_path = download.download(
            code=code, format=format, cache=cache, database=database
        )
        mol = cls.load(file_path, name=code)
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
        return self.object.modifiers["MolecularNodes"].node_group

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
        style: str | None = "spheres",
        color: str | None = "common",
        assembly: bool = False,
    ):
        """
        Add a style to the molecule.
        """
        bl.nodes.create_starting_node_tree(
            object=self.object, coll_frames=self.frames, style=style, color=color
        )

        if assembly:
            bl.nodes.assembly_initialise(self.object)
            bl.nodes.assembly_insert(self.object)

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
