from abc import ABCMeta

from ...blender import nodes
from .. import loading

class Molecule(metaclass=ABCMeta):
    """
    Abstract base class for representing a molecule.
    """

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
        
        mol, coll_frames = loading.create_model(
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
            mol['biological_assemblies'] = self.assemblies
        except InvalidFileError:
            pass
        
        if build_assembly:
            nodes.assembly_insert(mol)
        
        return mol