from abc import ABCMeta
import numpy as np

from ...blender import nodes
from .. import load

class Molecule(metaclass=ABCMeta):
    """
    Abstract base class for representing a molecule.
    """
    def entity_ids(self, as_int = False):
        if as_int:
            return np.unique(self.structure.entity_id, return_inverse = True)[1]
        return np.unique(self.structure.entity_id)
    
    def chain_ids(self, as_int = False):
        if as_int:
            return np.unique(self.structure.chain_id, return_inverse=True)[1]
        return np.unique(self.structure.chain_id)
    
    def create_model(
        self, 
        name: str = 'NewMolecule', 
        style: str = 'spheres', 
        build_assembly = False,
        centre: bool = False, 
        del_solvent: bool = True, 
        collection=None, 
        verbose: bool = False
        ) -> None:
        """
        Create a model in the 3D scene.

        Parameters
        ----------
        name : str, optional
            The name of the model. Default is 'NewMolecule'.
        centre : bool, optional
            Whether to center the model in the scene. Default is False.
        style : str, optional
            The style of the model. Default is 'spheres'.
        del_solvent : bool, optional
            Whether to delete solvent molecules. Default is True.
        collection : str, optional
            The collection to add the model to. Default is None.
        verbose : bool, optional
            Whether to print verbose output. Default is False.

        Returns
        -------
        None
        
        """
        from biotite import InvalidFileError
        
        mol, coll_frames = load.create_model(
            array=self.structure,
            name=name,
            centre=centre,
            del_solvent=del_solvent,
            style=style,
            collection=collection,
            verbose=verbose
        )
    
        if style:
            nodes.create_starting_node_tree(
                object = mol, 
                coll_frames=coll_frames, 
                style = style
                )
        
        try:
            mol['biological_assemblies'] = self.assemblies()
        except InvalidFileError:
            pass
        
        if build_assembly and style:
            nodes.assembly_insert(mol)
        
        return mol