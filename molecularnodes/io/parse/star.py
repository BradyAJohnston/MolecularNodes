import numpy as np
from .ensemble import Ensemble
from ... import blender as bl


class StarFile(Ensemble):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.data = self._read()
        self.star_type = None
        self.positions = None
        self._create_mn_columns()
        self.n_images = self._n_images()
        

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
                df['rlnCoordinateZ'] = 0

            self.positions = df[['rlnCoordinateX', 'rlnCoordinateY',
                      'rlnCoordinateZ']].to_numpy()
            pixel_size = df['rlnImagePixelSize'].to_numpy().reshape((-1, 1))
            self.positions = self.positions * pixel_size
            shift_column_names = ['rlnOriginXAngst',
                                  'rlnOriginYAngst', 'rlnOriginZAngst']
            if all([col in df.columns for col in shift_column_names]):
                shifts_ang = df[shift_column_names].to_numpy()
                self.positions = self.positions - shifts_ang
            df['MNAnglePhi'] = df['rlnAngleRot']
            df['MNAngleTheta'] = df['rlnAngleTilt']
            df['MNAnglePsi'] = df['rlnAnglePsi']
            df['MNImageId'] = df['rlnMicrographName'].astype(
                'category').cat.codes.to_numpy()
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
            df['MNImageId'] = df['cisTEMOriginalImageFilename'].astype(
                'category').cat.codes.to_numpy()


    def create_model(self, name='StarFileObject', node_setup=True, world_scale=0.01):

        object = bl.obj.create_object(
            self.positions * world_scale, collection=bl.coll.mn(), name=name)

        object.mn['molecule_type'] = 'star'
        object['mn_object'] = self

        # create attribute for every column in the STAR file
        for col in self.data.columns:
            col_type = self.data[col].dtype
            # If col_type is numeric directly add
            if np.issubdtype(col_type, np.number):
                bl.obj.set_attribute(
                    object, col, self.data[col].to_numpy().reshape(-1), 'FLOAT', 'POINT')

            # If col_type is object, convert to category and add integer values
            elif col_type == np.object:
                codes = self.data[col].astype(
                    'category').cat.codes.to_numpy().reshape(-1)
                bl.obj.set_attribute(object, col, codes, 'INT', 'POINT')
                # Add the category names as a property to the blender object
                object[f'{col}_categories'] = list(
                    self.data[col].astype('category').cat.categories)

        if node_setup:
            bl.nodes.create_starting_nodes_starfile(
                object, n_images=self.n_images)
            self.node_group = object.modifiers['MolecularNodes'].node_group

        self.object = object

        return object
