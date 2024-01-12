from abc import ABCMeta


class Ensemble(metaclass=ABCMeta):
    def __init__(self, file_path):
        """
        Initialize an Ensemble object.

        Parameters
        ----------
        file_path : str
            The path to the file.

        """
        self.type = "ensemble"
        self.file_path = file_path
        self.data_object = None
        self.data_collection = None

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
