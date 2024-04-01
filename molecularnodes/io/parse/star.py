import numpy as np
import bpy
from .ensemble import Ensemble
from ... import blender as bl

@bpy.app.handlers.persistent
def _rehydrate_ensembles(scene):
    for obj in bpy.data.objects:
        if hasattr(obj, 'mn') and 'molecule_type' in obj.mn.keys():
            if obj.mn['molecule_type'] == 'star':
                ensemble = StarFile.from_blender_object(obj)
                if not hasattr(bpy.types.Scene, 'MN_starfile_ensembles'):
                    bpy.types.Scene.MN_starfile_ensembles = []
                bpy.types.Scene.MN_starfile_ensembles.append(ensemble)

class StarFile(Ensemble):
    def __init__(self, file_path):
        super().__init__(file_path)
        
        
    @classmethod
    def from_starfile(cls, file_path):
        self = cls(file_path)
        self.data = self._read()
        self.star_type = None
        self.positions = None
        self.current_image = -1
        self._create_mn_columns()
        self.n_images = self._n_images()
        return self
    
    @classmethod
    def from_blender_object(cls, blender_object):
        import bpy
        self = cls(blender_object["starfile_path"])
        self.object = blender_object
        self.star_node = bl.nodes.get_star_node(self.object)
        self.micrograph_material = bl.nodes.MN_micrograph_material()
        self.data = self._read()
        self.star_type = None
        self.positions = None
        self.current_image = -1
        self._create_mn_columns()
        self.n_images = self._n_images()
        bpy.app.handlers.depsgraph_update_post.append(self._update_micrograph_texture)
        return self
        

    def _read(self):
        import starfile
        star = starfile.read(self.file_path)
        return star

    def _n_images(self):
        if isinstance(self.data, dict):
            return len(self.data)
        return 1

    def _create_mn_columns(self):
        # only RELION 3.1 and cisTEM STAR files are currently supported, fail gracefully
        if isinstance(self.data, dict) and 'particles' in self.data and 'optics' in self.data:
            self.star_type = 'relion'
        elif "cisTEMAnglePsi" in self.data:
            self.star_type = 'cistem'
        else:
            raise ValueError(
                'File is not a valid RELION>=3.1 or cisTEM STAR file, other formats are not currently supported.'
            )

        # Get absolute position and orientations
        if self.star_type == 'relion':
            df = self.data['particles'].merge(self.data['optics'], on='rlnOpticsGroup')

            # get necessary info from dataframes
            # Standard cryoEM starfile don't have rlnCoordinateZ. If this column is not present
            # Set it to "0"
            if "rlnCoordinateZ" not in df:
                if "rlnDefocusU" in df:
                    df['rlnCoordinateZ'] = (df['rlnDefocusU'] + df['rlnDefocusV']) / 2
                    df['rlnCoordinateZ'] = df['rlnCoordinateZ'] - df['rlnCoordinateZ'].median()
                else:
                    df['rlnCoordinateZ'] = 0

            self.positions = df[['rlnCoordinateX', 'rlnCoordinateY',
                      'rlnCoordinateZ']].to_numpy()
            if 'rlnMicrographOriginalPixelSize' in df:
                df['MNPixelSize'] = df['rlnMicrographOriginalPixelSize']
            else:
                df['MNPixelSize'] = df['rlnImagePixelSize']
            pixel_size = df['MNPixelSize'].to_numpy().reshape((-1, 1))
            self.positions = self.positions * pixel_size
            shift_column_names = ['rlnOriginXAngst',
                                  'rlnOriginYAngst', 'rlnOriginZAngst']
            if all([col in df.columns for col in shift_column_names]):
                shifts_ang = df[shift_column_names].to_numpy()
                self.positions -= shifts_ang
            df['MNAnglePhi'] = df['rlnAngleRot']
            df['MNAngleTheta'] = df['rlnAngleTilt']
            df['MNAnglePsi'] = df['rlnAnglePsi']
            
            try:
                df['MNImageId'] = df['rlnMicrographName'].astype(
                    'category').cat.codes.to_numpy()
            except KeyError:
                try:
                    df['MNImageId'] = df['rlnTomoName'].astype(
                        'category').cat.codes.to_numpy()
                except KeyError:
                    df['MNImageId'] = 0.0
            
            self.data = df

        elif self.star_type == 'cistem':
            df = self.data
            df['cisTEMZFromDefocus'] = (
                df['cisTEMDefocus1'] + df['cisTEMDefocus2']) / 2
            df['cisTEMZFromDefocus'] = df['cisTEMZFromDefocus'] - \
                df['cisTEMZFromDefocus'].median()
            self.positions = df[['cisTEMOriginalXPosition',
                      'cisTEMOriginalYPosition', 'cisTEMZFromDefocus']].to_numpy()
            df['MNAnglePhi'] = df['cisTEMAnglePhi']
            df['MNAngleTheta'] = df['cisTEMAngleTheta']
            df['MNAnglePsi'] = df['cisTEMAnglePsi']
            df['MNPixelSize'] = df['cisTEMPixelSize']
            df['MNImageId'] = df['cisTEMOriginalImageFilename'].astype(
                'category').cat.codes.to_numpy()
    
    def _convert_mrc_to_tiff(self):
        import mrcfile
        from pathlib import Path
        if self.star_type == 'relion':
            micrograph_path = self.object['rlnMicrographName_categories'][self.star_node.inputs['Image'].default_value - 1]
        elif self.star_type == 'cistem':
            micrograph_path = self.object['cisTEMOriginalImageFilename_categories'][self.star_node.inputs['Image'].default_value - 1].strip("'")
        else:
            return False
        
        # This could be more elegant
        if not Path(micrograph_path).exists():
            pot_micrograph_path = Path(self.file_path).parent / micrograph_path
            if not pot_micrograph_path.exists():
                if self.star_type == 'relion':
                    pot_micrograph_path = Path(self.file_path).parent.parent.parent / micrograph_path
                    if not pot_micrograph_path.exists():
                        raise FileNotFoundError(f"Micrograph file {micrograph_path} not found")
                else:
                    raise FileNotFoundError(f"Micrograph file {micrograph_path} not found")
            micrograph_path = pot_micrograph_path

        tiff_path = Path(micrograph_path).with_suffix('.tiff')
        if not tiff_path.exists():
            with mrcfile.open(micrograph_path) as mrc:
                micrograph_data = mrc.data.copy()

            # For 3D data sum over the z axis. Probalby would be nicer to load the data as a volume
            if micrograph_data.ndim == 3:
                micrograph_data = np.sum(micrograph_data, axis=0)
            # Normalize the data to 0-1
            micrograph_data = (micrograph_data - micrograph_data.min()) / (micrograph_data.max() - micrograph_data.min())
            
            if micrograph_data.dtype != np.float32:
                micrograph_data = micrograph_data.astype(np.float32)
            from PIL import Image
            # Need to invert in Y to generate the correct tiff
            Image.fromarray(micrograph_data[::-1,:]).save(tiff_path)
        return tiff_path
    
    def _update_micrograph_texture(self, *_):
        try:
            show_micrograph = self.star_node.inputs['Show Micrograph']
            _ = self.object['mn']
        except ReferenceError:
            bpy.app.handlers.depsgraph_update_post.remove(self._update_micrograph_texture)
            return
        if self.star_node.inputs['Image'].default_value == self.current_image:
            return
        else:
            self.current_image = self.star_node.inputs['Image'].default_value
        if not show_micrograph:
            return
        tiff_path = self._convert_mrc_to_tiff()
        if tiff_path:
            try:
                image_obj = bpy.data.images[tiff_path.name]
            except KeyError:
                image_obj = bpy.data.images.load(str(tiff_path))
            image_obj.colorspace_settings.name = 'Non-Color'
            self.micrograph_material.node_tree.nodes['Image Texture'].image = image_obj
            self.star_node.inputs['Micrograph'].default_value = image_obj

                 

    def create_model(self, name='StarFileObject', node_setup=True, world_scale=0.01):
        from molecularnodes.blender.nodes import get_star_node, MN_micrograph_material
        blender_object = bl.obj.create_object(
            self.positions * world_scale, collection=bl.coll.mn(), name=name)

        blender_object.mn['molecule_type'] = 'star'
        
        # create attribute for every column in the STAR file
        for col in self.data.columns:
            col_type = self.data[col].dtype
            # If col_type is numeric directly add
            if np.issubdtype(col_type, np.number):
                bl.obj.set_attribute(
                    blender_object, col, self.data[col].to_numpy().reshape(-1), 'FLOAT', 'POINT')

            # If col_type is object, convert to category and add integer values
            elif col_type == object:
                codes = self.data[col].astype(
                    'category').cat.codes.to_numpy().reshape(-1)
                bl.obj.set_attribute(blender_object, col, codes, 'INT', 'POINT')
                # Add the category names as a property to the blender object
                blender_object[f'{col}_categories'] = list(
                    self.data[col].astype('category').cat.categories)

        if node_setup:
            bl.nodes.create_starting_nodes_starfile(
                blender_object, n_images=self.n_images)
            self.node_group = blender_object.modifiers['MolecularNodes'].node_group

        blender_object["starfile_path"] = str(self.file_path)
        self.object = blender_object
        self.star_node = get_star_node(self.object)
        self.micrograph_material = MN_micrograph_material()
        bpy.app.handlers.depsgraph_update_post.append(self._update_micrograph_texture)
        return blender_object
