from abc import ABCMeta
from pathlib import Path
from typing import Union
import bpy
from ... import blender as bl
from ..base import EntityType, MolecularEntity


class Ensemble(MolecularEntity, metaclass=ABCMeta):
    def __init__(self, file_path: Union[str, Path]):
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

    def create_object(
        self,
        name: str = "NewEnsemble",
        node_setup: bool = True,
        world_scale: float = 0.01,
        fraction: float = 1.0,
        simplify=False,
    ) -> bpy.types.Object:
        """
        Create a 3D object for the ensemble.

        Parameters
        ----------
        name : str, optional
            The name of the model, by default "NewEnsemble"
        node_setup : bool, optional
            Whether to setup nodes for the data and instancing objects, by default True
        world_scale : float, optional
            Scaling transform for the coordinates before loading in to Blender, by default 0.01
        fraction : float, optional
            The fraction of the instances to display on loading. Reducing can help with performance, by default 1.0
        simplify : bool, optional
            Whether to instance the given models or simplify them for debugging and performance, by default False

        Notes
        -----
        Creates a data object which stores all of the required instancing information. If
        there are molecules to be instanced, they are also created in their own data collection.
        """
        pass
