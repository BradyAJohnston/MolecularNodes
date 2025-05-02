import os
from pathlib import Path
import bpy
import databpy
import mrcfile
import numpy as np
from ...blender import coll
from ...nodes import nodes
from .base import Density


class MRC(Density):
    """
    A class for parsing EM density files in the format `.map` or `.map.gz`.

    It utilises `mrcfile` for file parsing, which is then converted into `pyopevdb` grids,
    that can be written as `.vdb` files and the imported into Blender as volumetric objects.
    """

    def __init__(
        self,
        file_path: str | Path,
        center: bool = False,
        invert: bool = False,
        overwrite: bool = False,
    ):
        super().__init__(file_path=file_path)
        self.grid = self.map_to_grid(str(self.file_path), center=center)
        self.file_vdb = self.map_to_vdb(
            str(self.file_path), center=center, invert=invert, overwrite=overwrite
        )

    def create_object(
        self, name="NewDensity", style="density_surface", setup_nodes=True
    ) -> bpy.types.Object:
        """
        Loads an MRC file into Blender as a volumetric object.

        Parameters
        ----------
        name : str, optional
            If not empty, renames the object with the new name. Default is "NewDensity".
        style : str, optional
            The style of the density object. Default is "density_surface".
        setup_nodes : bool, optional
            Whether to create starting node tree. Default is True.

        Returns
        -------
        bpy.types.Object
            The loaded volumetric object.
        """
        # import and ensure object is at world origin to get corect alignment with
        # structures
        self.object = databpy.import_vdb(self.file_vdb, collection=coll.mn())
        self.object.location = (0, 0, 0)
        self.object.mn.entity_type = self._entity_type.value

        if name and name != "":
            self.name = name

        self.create_starting_node_tree(style=style)

        return self.object

    def create_starting_node_tree(self, style="density_surface"):
        """
        Creates a starting node tree for the density object.

        Parameters
        ----------
        style : str, optional
            The style of the density object, defaulting to 'density_surface'.
        """
        nodes.create_starting_nodes_density(
            object=self.object, style=style, threshold=self.threshold
        )

    def map_to_vdb(
        self,
        file: str,
        invert: bool = False,
        world_scale=0.01,
        center: bool = False,
        overwrite=False,
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
        import openvdb as vdb

        file_path = self.path_to_vdb(file, center=center, invert=invert)

        # If the map has already been converted to a .vdb and overwrite is False, return that instead
        if os.path.exists(file_path) and not overwrite:
            # Also check that the file has the same invert and center settings
            grid = vdb.readAllGridMetadata(file_path)[0]
            if (
                "MN_invert" in grid
                and grid["MN_invert"] == invert
                and "MN_center" in grid
                and grid["MN_center"] == center
            ):
                self.threshold = grid["MN_initial_threshold"]
                return file_path

        print("Reading new file")
        # Read in the MRC file and convert it to a pyopenvdb grid
        grid = self.map_to_grid(file=file, invert=invert, center=center)

        print(f"{dir(grid.transform)}")
        grid.transform.preScale(
            np.array((1, 1, 1)) * world_scale * grid["MN_voxel_size"]
        )

        if center:
            offset = -np.array(grid["MN_box_size"]) * 0.5
            offset *= grid["MN_voxel_size"] * world_scale
            print("transforming")
            grid.transform.postTranslate(offset)

        if os.path.exists(file_path):
            os.remove(file_path)

        # Write the grid to a .vdb file
        print("writing new file")
        vdb.write(file_path, grids=[grid])
        self.threshold = grid["MN_initial_threshold"]
        del grid

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

        import openvdb as vdb

        volume = mrcfile.read(file)

        dataType = volume.dtype

        # enables different grid types

        if dataType == "float32" or dataType == "float64":
            grid = vdb.FloatGrid()
        elif dataType == "int8" or dataType == "int16" or dataType == "int32":
            volume = volume.astype("int32")
            grid = vdb.Int32Grid()
        elif dataType == "int64":
            grid = vdb.Int64Grid()
        else:
            grid = vdb.FloatGrid()

        if invert:
            volume = np.max(volume) - volume

        initial_threshold = np.quantile(volume, 0.995)

        # The np.copy is needed to force numpy to actually rewrite the data in memory
        # since openvdb seems to read is straight from memory without checking the striding
        # The np.transpose is needed to convert the data from zyx to xyz
        volume = np.copy(np.transpose(volume, (2, 1, 0)), order="C")
        try:
            grid.copyFromArray(volume.astype(float))
        except Exception as e:
            print(
                f"Grid data type '{volume.dtype}' is an unsupported type.\nError: {e}"
            )

        grid.gridClass = vdb.GridClass.FOG_VOLUME
        grid.name = "density"

        # Set some metadata for the vdb file, so we can check if it's already been converted
        # correctly
        grid["MN_invert"] = invert
        grid["MN_initial_threshold"] = initial_threshold
        grid["MN_center"] = center
        with mrcfile.open(file) as mrc:
            grid["MN_voxel_size"] = float(mrc.voxel_size.x)
            grid["MN_box_size"] = (
                int(mrc.header.nx),
                int(mrc.header.ny),
                int(mrc.header.nz),
            )

        return grid
