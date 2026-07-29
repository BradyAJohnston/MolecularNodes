from abc import ABCMeta, abstractmethod
from pathlib import Path
from typing import Self
import bpy
from ... import blender as bl
from ..base import EntityType, MolecularEntity


class Ensemble(MolecularEntity, metaclass=ABCMeta):
    def __init__(self, file_path: str | Path):
        """
        Initialize an Ensemble object.

        Parameters
        ----------
        file_path : Union[str, Path]
            The path to the file.
        """
        super().__init__()
        self._entity_type = EntityType.ENSEMBLE
        self.type: str = "ensemble"
        self.file_path: Path = bl.path_resolve(file_path)

    @classmethod
    def load(
        cls, file_path: str | Path, name: str | None = None, node_setup: bool = True
    ) -> Self:
        """
        Load an ensemble from a file and create its 3D object.

        Parameters
        ----------
        file_path : str or Path
            The path to the ensemble file, such as a CellPack ``.cif``/``.bcif`` model
            or a RELION/cisTEM ``.star`` file.
        name : str or None, optional
            The name to give the created object. If None, the file name is used.
        node_setup : bool, optional
            Whether to set up the geometry nodes for the ensemble, by default True.

        Returns
        -------
        Ensemble
            The loaded ensemble, with its data and 3D representation created.
        """
        self = cls(file_path)
        self.create_object(name=name or Path(file_path).name, node_setup=node_setup)
        # record the source so the entity can be reloaded into a fresh session
        self.object.mn.filepath = str(file_path)
        return self

    @property
    def instance_collection(self) -> bpy.types.Collection:
        """
        The instances of the ensemble.

        Returns
        -------
        bpy.types.Collection
            The collection containing the ensemble instances.
        """
        return bpy.data.collections[self._instance_collection_name]

    @instance_collection.setter
    def instance_collection(self, value: bpy.types.Collection) -> None:
        """
        Set the instance collection.

        Parameters
        ----------
        value : bpy.types.Collection
            The collection to set as the instance collection.

        Raises
        ------
        ValueError
            If the value is not a bpy.types.Collection.
        """
        if not isinstance(value, bpy.types.Collection):
            raise ValueError("The instances must be a bpy.types.Collection.")
        self._instance_collection_name = value.name

    @abstractmethod
    def create_object(
        self,
        name: str = "NewEnsemble",
        node_setup: bool = True,
    ) -> bpy.types.Object:
        """
        Create a 3D object for the ensemble.

        Parameters
        ----------
        name : str, optional
            The name of the model, by default "NewEnsemble"
        node_setup : bool, optional
            Whether to setup nodes for the data and instancing objects, by default True

        Notes
        -----
        Creates a data object which stores all of the required instancing information. If
        there are molecules to be instanced, they are also created in their own data collection.
        """
        ...
