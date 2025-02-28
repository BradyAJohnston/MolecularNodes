import io
from abc import ABCMeta
from pathlib import Path
import json

import biotite.structure as struc
from biotite.structure import AtomArray, AtomArrayStack
import bpy
import numpy as np
from biotite import InvalidFileError
from . import pdbx, pdb, sdf


from ... import blender as bl
from ... import utils, download
from databpy import Domains, AttributeTypes
import databpy
from ..base import MolecularEntity, EntityType


class Molecule(MolecularEntity, metaclass=ABCMeta):
    """
    Abstract base class for representing a molecule.

    It associates the atomic data (the array) with the created 3D model inside of Blender
    (the object). If multiple conformations are imported, then a `frames` collection
    is also instantiated.

    The `named_attribute()` and `store_named_attribute()` methods access and set attributes on
    `object` that is in the Blender scene.

    Attributes
    ----------
    file_path : str
        The file path to the file which stores the atomic coordinates.
    file : Any
        The opened file.
    object : bpy.types.Object
        The Blender object representing the molecule.
    frames : bpy.types.Collection
        The Blender collection which holds the objects making up the frames to animate.
    array: np.ndarray:
        The numpy array which stores the atomic coordinates and associated attributes.
    entity_ids : np.ndarray
        The entity IDs of the molecule.
    chain_ids : np.ndarray
        The chain IDs of the molecule.
    """

    def __init__(self, file_path: str | Path | io.BytesIO):
        """
        Initialize the Molecule object.

        Parameters
        ----------
        file_path : Union[str, Path, io.BytesIO]
            The file path to the file which stores the atomic coordinates.
        """
        self._frames_collection: str | None = None
        self._entity_type = EntityType.MOLECULE
        super().__init__()
        self._read(file_path)
        self.atom_array: AtomArray | AtomArrayStack
        self._create_object()
        self._code: str | None = None

    @property
    def code(self) -> str | None:
        """
        Get the code for the molecule.
        """
        return self._code

    def _read(self, file_path: str | Path | io.BytesIO) -> None:
        """
        Initially open the file, ready to extract the required data.

        Parameters
        ----------
        file_path : Union[Path, io.BytesIO]
            The file path to the file which stores the atomic coordinates.
        """
        if isinstance(file_path, io.BytesIO):
            self._reader = pdbx.PDBXReader(file_path)
        else:
            if isinstance(file_path, str):
                file_path = Path(file_path)

            match file_path.suffix:
                case ".cif":
                    self._reader = pdbx.PDBXReader(file_path)
                case ".bcif":
                    self._reader = pdbx.PDBXReader(file_path)
                case ".pdb":
                    self._reader = pdb.PDBReader(file_path)
                case ".sdf":
                    self._reader = sdf.SDFReader(file_path)
                case _:
                    raise InvalidFileError("The file format is not supported.")

        self.atom_array = self._reader.array.copy()
        self._assemblies = self._reader._assemblies()
        del self._reader

    @classmethod
    def fetch(
        cls, code: str, format=".bcif", cache=download.CACHE_DIR, database: str = "rcsb"
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
        mol = cls(file_path)
        mol.object.mn["entity_type"] = "molecule"
        mol._code = code
        return mol

    @property
    def frames(self) -> bpy.types.Collection:
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

    def __len__(self) -> int:
        """
        Get the number of atoms in the molecule.

        Returns
        -------
        int
            The number of atoms in the molecule.
        """
        if hasattr(self, "object"):
            if self.object:
                return len(self.object.data.vertices)
        if self.atom_array:
            return len(self.atom_array)
        else:
            return 0

    @property
    def n_models(self):
        """
        Get the number of models in the molecule.

        Returns
        -------
        int
            The number of models in the molecule.
        """
        if isinstance(self.atom_array, struc.AtomArray):
            return 1
        else:
            return self.atom_array.shape[0]

    def add_style(self, style: str | None = "spheres", color: str | None = None):
        """
        Add a style to the molecule.
        """
        bl.nodes.create_starting_node_tree(
            object=self.object, coll_frames=self.frames, style=style, color=color
        )

    def _store_all_attributes(self):
        """We store all potential attributes form the AtomArray onto the mesh in Blender."""

        annotations_to_skip = ["coord"]

        for attr in self.atom_array.get_annotation_categories():
            if attr in annotations_to_skip:
                continue

            data = self.atom_array.get_annotation(attr)  # type: ignore

            if attr == "vdw_radii":
                data *= self._world_scale

            # geometry nodes doesn't support strings at the moment so we can only store the
            # numeric and boolean attributes on the mesh. All string attributes should have
            # already been converted to a corresponding numeric attribute during the
            # reader process

            if not (
                np.issubdtype(data.dtype, np.number)
                or np.issubdtype(data.dtype, np.bool_)
            ):
                continue

            # the integer versions of strings have been added as annotations that just
            # append `_int` onto the name so we don't overrite the original data
            self.store_named_attribute(
                data=data,
                name=attr.replace("_int", ""),
            )

    def _setup_frames_collection(self):
        if self.n_models > 1:
            self.frames = bl.coll.frames(self.name)
            for i, array in enumerate(self.atom_array):
                databpy.create_object(
                    vertices=array.coord * self._world_scale,  # type: ignore
                    name="{}_frame_{}".format(self.name, str(i)),
                    collection=self.frames,
                )

    def _store_object_properties(self):
        try:
            self.object["entity_ids"] = np.unique(self.atom_array.entity_id).tolist()  # type: ignore
        except AttributeError:
            self.object["entity_ids"] = None

        try:
            self.object["chain_ids"] = np.unique(self.atom_array.chain_id).tolist()  # type: ignore
        except AttributeError:
            self.object["chain_ids"] = None

        try:
            self.object.mn.biological_assemblies = json.dumps(self.assemblies())  # type: ignore
        except InvalidFileError:
            self.object.mn.biological_assemblies = ""  # type: ignore

    def _create_object(
        self,
        name="NewObject",
    ) -> None:
        """
        Create a 3D model of the molecule, with one vertex for each atom.

        Parameters
        ----------
        verbose : bool, optional
            Whether to print verbose output. Default is False.

        Returns
        -------
        bpy.types.Object
            The created 3D model, as an object in the 3D scene.
        """

        if isinstance(self.atom_array, AtomArrayStack):
            is_stack = True
            array = self.atom_array[0]
        else:
            is_stack = False
            array = self.atom_array

        self.object = databpy.create_object(
            vertices=array.coord * self._world_scale,  # type: ignore
            edges=array.bonds.as_array()[:, :2] if array.bonds is not None else None,
            name=name,
        )
        # Add information about the bond types to the model on the edge domain
        # Bond types: 'ANY' = 0, 'SINGLE' = 1, 'DOUBLE' = 2, 'TRIPLE' = 3, 'QUADRUPLE' = 4
        # 'AROMATIC_SINGLE' = 5, 'AROMATIC_DOUBLE' = 6, 'AROMATIC_TRIPLE' = 7
        # https://www.biotite-python.org/apidoc/biotite.structure.BondType.html#biotite.structure.BondType
        if array.bonds:
            self.store_named_attribute(
                array.bonds.as_array()[:, 2], "bond_type", domain=Domains.EDGE
            )
        if is_stack:
            self._setup_frames_collection()

        self._store_all_attributes()
        self._store_object_properties()

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
            assemblies_info = self._assemblies
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
