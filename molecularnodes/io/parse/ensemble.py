import bpy
from abc import ABCMeta
import numpy as np
from ... import blender as bl
import warnings


class Ensemble(metaclass=ABCMeta):
    def __init__(self, file_path):
        """
        Initialize an Ensemble object.

        Parameters
        ----------
        file_path : str
            The path to the file.

        """
        self.type: str = "ensemble"
        self.file_path: str = file_path
        self.object: bpy.types.Object = None
        self.instances: bpy.types.Collection = None
        self.frames: bpy.types.Collection = None

    @classmethod
    def create_model(cls, name: str = "NewEnsemble", node_setup: bool = True, world_scale: float = 0.01, fraction: float = 1.0, simplify=False):
        """
        Create a 3D model in the of the ensemble.

        Parameters
        ----------
        name : str, optional
            The name of the model (default is "NewEnsemble").
        node_setup : bool, optional
            Whether to setup nodes for the data and instancing objects. (default is True).
        world_scale : float, optional
            Scaling transform for the coordinates before loading in to Blender. (default is 0.01).
        fraction : float, optional
            The fraction of the instances to display on loading. Reducing can help with performance. (default is 1.0).
        simplify : bool, optional
            Whether to isntance the given models or simplify them for debugging and performance. (default is False).

        Creates a data object which stores all of the required instancing information. If 
        there are molecules to be instanced, they are also created in their own data collection.

        Parameters:
        - name (str): The name of the model. Default is "NewEnsemble".
        - node_setup (bool): Whether to set up nodes. Default is True.
        - world_scale (float): The scale of the world. Default is 0.01.
        - fraction (float): The fraction of molecules to be instanced. Default is 1.0.
        - simplify (bool): Whether to simplify the model. Default is False.

        """
        pass

    def get_attribute(self, name='position', evaluate=False) -> np.ndarray | None:
        """
        Get the value of an object for the data molecule.

        Parameters
        ----------
        name : str, optional
            The name of the attribute. Default is 'position'.
        evaluate : bool, optional
            Whether to first evaluate all node trees before getting the requsted attribute. 
            False (default) will sample the underlying atomic geometry, while True will 
            sample the geometry that is created through the Geometry Nodes tree.

        Returns
        -------
        np.ndarray
            The value of the attribute.
        """
        if not self.object:
            warnings.warn(
                'No object yet created. Use `create_model()` to create a corresponding object.'
            )
            return None
        return bl.obj.get_attribute(self.object, name=name, evaluate=evaluate)
