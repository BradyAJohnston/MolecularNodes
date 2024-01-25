import numpy as np
from .ensemble import Ensemble
from ... import blender as bl


class StarFile(Ensemble):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.data = self._read()
        self.n_images = self._n_images()

    def _read(self):
        import starfile
        star = starfile.read(self.file_path)
        return star

    def _n_images(self):
        if isinstance(self.data, dict):
            return len(self.data)
        return 1

    def create_model(self, name='StarFileObject', node_setup=True, world_scale=0.01):

        star = self.data

        # only RELION 3.1 and cisTEM STAR files are currently supported, fail gracefully
        if isinstance(star, dict) and 'particles' in star and 'optics' in star:
            star_type = 'relion'
        elif "cisTEMAnglePsi" in star:
            star_type = 'cistem'
        else:
            raise ValueError(
                'File is not a valid RELION>=3.1 or cisTEM STAR file, other formats are not currently supported.'
            )

        # Get absolute position and orientations
        if star_type == 'relion':
            df = star['particles'].merge(star['optics'], on='rlnOpticsGroup')

            # get necessary info from dataframes
            # Standard cryoEM starfile don't have rlnCoordinateZ. If this column is not present
            # Set it to "0"
            if "rlnCoordinateZ" not in df:
                df['rlnCoordinateZ'] = 0

            xyz = df[['rlnCoordinateX', 'rlnCoordinateY',
                      'rlnCoordinateZ']].to_numpy()
            pixel_size = df['rlnImagePixelSize'].to_numpy().reshape((-1, 1))
            xyz = xyz * pixel_size
            shift_column_names = ['rlnOriginXAngst',
                                  'rlnOriginYAngst', 'rlnOriginZAngst']
            if all([col in df.columns for col in shift_column_names]):
                shifts_ang = df[shift_column_names].to_numpy()
                xyz = xyz - shifts_ang
            df['MNAnglePhi'] = df['rlnAngleRot']
            df['MNAngleTheta'] = df['rlnAngleTilt']
            df['MNAnglePsi'] = df['rlnAnglePsi']
            image_id = df['rlnMicrographName'].astype(
                'category').cat.codes.to_numpy()

        elif star_type == 'cistem':
            df = star
            df['cisTEMZFromDefocus'] = (
                df['cisTEMDefocus1'] + df['cisTEMDefocus2']) / 2
            df['cisTEMZFromDefocus'] = df['cisTEMZFromDefocus'] - \
                df['cisTEMZFromDefocus'].median()
            xyz = df[['cisTEMOriginalXPosition',
                      'cisTEMOriginalYPosition', 'cisTEMZFromDefocus']].to_numpy()
            df['MNAnglePhi'] = df['cisTEMAnglePhi']
            df['MNAngleTheta'] = df['cisTEMAngleTheta']
            df['MNAnglePsi'] = df['cisTEMAnglePsi']
            image_id = df['cisTEMOriginalImageFilename'].astype(
                'category').cat.codes.to_numpy()

        object = bl.obj.create_object(
            xyz * world_scale, collection=bl.coll.mn(), name=name)

        object.mn['molecule_type'] = 'star'
        object.mn['star_type'] = star_type

        # create the attribute and add the data for the image id
        bl.obj.set_attribute(object, 'MNImageId', image_id, 'INT', 'POINT')

        # create attribute for every column in the STAR file
        for col in df.columns:
            col_type = df[col].dtype
            # If col_type is numeric directly add
            if np.issubdtype(col_type, np.number):
                bl.obj.set_attribute(
                    object, col, df[col].to_numpy().reshape(-1), 'FLOAT', 'POINT')

            # If col_type is object, convert to category and add integer values
            elif col_type == object:
                codes = df[col].astype(
                    'category').cat.codes.to_numpy().reshape(-1)
                bl.obj.set_attribute(object, col, codes, 'INT', 'POINT')
                # Add the category names as a property to the blender object
                object[f'{col}_categories'] = list(
                    df[col].astype('category').cat.categories)

        if node_setup:
            bl.nodes.create_starting_nodes_starfile(
                object, n_images=self.n_images)
            self.node_group = object.modifiers['MolecularNodes'].node_group

        self.object = object

        return object
