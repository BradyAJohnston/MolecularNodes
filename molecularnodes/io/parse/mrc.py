from .density import Density

from ...blender import coll, obj, nodes
import bpy
import numpy as np
import os


class MRC(Density):
    """
    A class for parsing EM density files in the format `.map` or `.map.gz`.

    It utilises `mrcfile` for file parsing, which is then converted into `pyopevdb` grids, 
    that can be written as `.vdb` files and the imported into Blender as volumetric objects.
    """

    def __init__(self, file_path, center=False, invert=False):
        super().__init__(self)
        self.file_path = file_path
        self.grid = self.map_to_grid(self.file_path, center=center)
        self.file_vdb = self.map_to_vdb(
            self.file_path, center=center, invert=invert)

    def create_model(
        self,
        name='NewDensity',
        style='density_surface',
        setup_nodes=True,
        invert: bool = False,
        center: bool = False,
        world_scale: float = 0.01
    ) -> bpy.types.Object:
        """
        Loads an MRC file into Blender as a volumetric object.

        Parameters
        ----------
        file : str
            Path to the MRC file.
        name : str, optional
            If not None, renames the object with the new name.
        invert : bool, optional
            Whether to invert the data from the grid, defaulting to False. Some file types
            such as EM tomograms have inverted values, where a high value == low density.
        world_scale : float, optional
            Scale of the object in the world. Defaults to 0.01.
        center : bool, optional
            Whether to center the volume on the origin. Defaults to False.

        Returns
        -------
        bpy.types.Object
            The loaded volumetric object.
        """
        object = obj.import_vdb(self.file_vdb, collection=coll.mn())
        self.object = object
        object.mn['molecule_type'] = 'density'

        if name and name != "":
            # Rename object to specified name
            object.name = name

        if setup_nodes:
            nodes.create_starting_nodes_density(
                object, style=style, threshold=self.threshold)

        return object

    def map_to_vdb(
        self,
        file: str,
        invert: bool = False,
        world_scale=0.01,
        center: bool = False,
        overwrite=False
    ) -> (str, float):
        """
        Converts an MRC file to a .vdb file using pyopenvdb.

        Parameters
        ----------
        file : str
            The path to the input MRC file.
        invert : bool, optional
            Whether to invert the data from the grid, defaulting to False. Some file types
            such as EM tomograms have inverted values, where a high value == low density.
        world_scale : float, optional
            The scaling factor to apply to the voxel size of the input file. Defaults to 0.01.
        center : bool, optional
            Whether to center the volume on the origin. Defaults to False.
        overwrite : bool, optional
            If True, the .vdb file will be overwritten if it already exists. Defaults to False.

        Returns
        -------
        str
            The path to the converted .vdb file.
        """
        import pyopenvdb as vdb

        file_path = self.path_to_vdb(file)

        # If the map has already been converted to a .vdb and overwrite is False, return that instead
        if os.path.exists(file_path) and not overwrite:
            # Also check that the file has the same invert and center settings
            grid = vdb.readAllGridMetadata(file_path)[0]
            if 'MN_invert' in grid and grid['MN_invert'] == invert and 'MN_center' in grid and grid['MN_center'] == center:
                self.threshold = grid['MN_initial_threshold']
                return file_path

        # Read in the MRC file and convert it to a pyopenvdb grid
        grid = self.map_to_grid(file, invert=invert, center=center)

        grid.transform.scale(np.array((1, 1, 1)) *
                             world_scale * grid['MN_voxel_size'])
        if center:
            grid.transform.translate(-np.array(grid['MN_box_size'])
                                     * 0.5 * world_scale * grid['MN_voxel_size'])

        # Write the grid to a .vdb file
        vdb.write(file_path, grids=[grid])
        self.threshold = grid['MN_initial_threshold']

        # Return the path to the output file
        return file_path

    def map_to_grid(self, file: str, invert: bool = False, center: bool = False):
        """
        Reads an MRC file and converts it into a pyopenvdb FloatGrid object.

        This function reads a file in MRC format, and converts it into a pyopenvdb FloatGrid object,
        which can be used to represent volumetric data in Blender.

        Parameters
        ----------
        file : str
            The path to the MRC file.
        invert : bool, optional
            Whether to invert the data from the grid, defaulting to False. Some file types
            such as EM tomograms have inverted values, where a high value == low density.

        Returns
        -------
        pyopenvdb.FloatGrid
            A pyopenvdb FloatGrid object containing the density data.
        """
        import mrcfile
        import pyopenvdb as vdb

        volume = mrcfile.read(file)

        dataType = volume.dtype

        # enables different grid types

        if dataType == "float32" or dataType == "float64":
            grid = vdb.FloatGrid()
        elif dataType == "int8" or dataType == "int16" or dataType == "int32":
            volume = volume.astype('int32')
            grid = vdb.Int32Grid()
        elif dataType == "int64":
            grid = vdb.Int64Grid()

        if invert:
            volume = np.max(volume) - volume

        initial_threshold = np.quantile(volume, 0.995)

        # The np.copy is needed to force numpy to actually rewrite the data in memory
        # since openvdb seems to read is straight from memory without checking the striding
        # The np.transpose is needed to convert the data from zyx to xyz
        volume = np.copy(np.transpose(volume, (2, 1, 0)), order='C')
        try:
            grid.copyFromArray(volume)
        except ValueError:
            print(f"Grid data type '{volume.dtype}' is an unsupported type.")

        grid.gridClass = vdb.GridClass.FOG_VOLUME
        grid.name = 'density'

        # Set some metadata for the vdb file, so we can check if it's already been converted
        # correctly
        grid['MN_invert'] = invert
        grid['MN_initial_threshold'] = initial_threshold
        grid['MN_center'] = center
        with mrcfile.open(file) as mrc:
            grid['MN_voxel_size'] = float(mrc.voxel_size.x)
            grid['MN_box_size'] = (int(mrc.header.nx), int(
                mrc.header.ny), int(mrc.header.nz))

        return grid
