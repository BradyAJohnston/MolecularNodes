import os
from pathlib import Path
import bpy
import databpy
import numpy as np
from gridData import Grid
from ...blender import coll
from ...nodes import nodes
from .base import Density


class Grids(Density):
    """
    A class for parsing grid files in the various formats:
    See: https://www.mdanalysis.org/GridDataFormats/gridData/formats.html#supported-file-formats
    Plain: `.dx`, `.plt`, `.ccp4`, `.mrc`, `.map`, `.pickle`
    Compressed: `.dx.gz`, `.ccp4.gz`, `.ccp4.bz2`, `.mrc.gz`, `.mrc.bz2`, `.map.gz`, `.map.bz2`

    It utilises `GridDataFormats` for file parsing, which is then converted into `openvdb` grids,
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
        self.grid: Grid = None  # GridDataFormats Grid object
        self.file_vdb = self.grid_to_vdb(
            str(self.file_path), center=center, invert=invert, overwrite=overwrite
        )

    def _parse_grid_with_fallback(
        self,
        file: str,
        file_format: str | None,
        metadata: dict,
    ):
        """
        Try parsing the grid using GridDataFormats first; on failure, fall back to
        reading with the mrcfile package and constructing a minimal grid-like object.

        Returns an object with attributes: grid (np.ndarray), delta (np.ndarray),
        origin (np.ndarray), metadata (dict).
        """
        # First attempt: GridDataFormats
        try:
            return Grid(file, file_format=file_format, metadata=metadata)
        except Exception:
            pass

        # Fallback: try using mrcfile for MRC/CCP4/MAP formats
        try:
            import mrcfile  # type: ignore

            with mrcfile.open(file, permissive=True) as mrc:
                # Determine voxel size (delta)
                try:
                    vx, vy, vz = (
                        float(mrc.voxel_size.x),
                        float(mrc.voxel_size.y),
                        float(mrc.voxel_size.z),
                    )
                except Exception:
                    # Compute from header cell dimensions if voxel_size unavailable
                    h = mrc.header
                    nx, ny, nz = int(h.nx), int(h.ny), int(h.nz)
                    # cella stored in Angstroms
                    ax, ay, az = float(h.cella.x), float(h.cella.y), float(h.cella.z)
                    vx = ax / nx if nx else 1.0
                    vy = ay / ny if ny else 1.0
                    vz = az / nz if nz else 1.0

                # Origin (Angstroms). If not set, default to (0,0,0)
                try:
                    ox, oy, oz = (
                        float(mrc.header.origin.x),
                        float(mrc.header.origin.y),
                        float(mrc.header.origin.z),
                    )
                except Exception:
                    ox = oy = oz = 0.0

                # The OpenVDB Python API expects array in Z, Y, X order by default.
                # mrcfile returns data with shape (nz, ny, nx), which matches that order.
                delta = np.array([vx, vy, vz], dtype=float)
                origin = np.array([ox, oy, oz], dtype=float)

                return Grid(mrc.data, origin=origin, delta=delta, metadata=metadata)
        except Exception as e:
            # Re-raise the original parsing error with context from fallback
            raise RuntimeError(
                f"Failed to parse grid file '{file}' with GridDataFormats and mrcfile fallback: {e}"
            )

    def create_object(
        self, name="NewDensity", style="density_surface", setup_nodes=True
    ) -> bpy.types.Object:
        """
        Loads a grid into Blender as a volumetric object.

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

        # TODO: setup_nodes is unused. Using it causes pytest failures
        node_density = self.create_starting_node_tree(style=style)
        # set the active index for UI to the added style
        self.object.mn.styles_active_index = self.tree.nodes.find(node_density.name)

        return self.object

    def create_starting_node_tree(self, style="density_surface"):
        """
        Creates a starting node tree for the density object.

        Parameters
        ----------
        style : str, optional
            The style of the density object, defaulting to 'density_surface'.
        """

        gobj = self.grid
        grid = gobj.grid
        threshold = np.quantile(grid, 0.995)
        threshold_range = (np.min(grid), np.max(grid))
        threshold_type = None
        if np.issubdtype(grid.dtype, np.floating):
            threshold_type = "NodeSocketFloat"
        elif np.issubdtype(grid.dtype, np.floating):
            threshold_type = "NodeSocketInt"

        x_range = y_range = z_range = None
        if gobj.origin.size == 3:
            origin = gobj.origin.copy()
            if gobj.metadata["center"]:
                origin = -np.array(grid.shape) * 0.5 * gobj.delta
            origin *= self._world_scale
            length = grid.shape * gobj.delta * self._world_scale
            x_range = (origin[0], origin[0] + length[0])
            y_range = (origin[1], origin[1] + length[1])
            z_range = (origin[2], origin[2] + length[2])

        return nodes.create_starting_nodes_density(
            object=self.object,
            style=style,
            threshold=threshold,
            threshold_range=threshold_range,
            threshold_type=threshold_type,
            x_range=x_range,
            y_range=y_range,
            z_range=z_range,
        )

    def grid_to_vdb(
        self,
        file: str,
        invert: bool = False,
        world_scale=0.01,
        center: bool = False,
        overwrite=False,
    ) -> str:
        """
        Converts an grid file to a .vdb file using openvdb.

        Parameters
        ----------
        file : str
            The path to the input grid file.
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

        bpy.utils.expose_bundled_modules()

        is_pyopenvdb = False
        if bpy.app.version >= (4, 4, 0):
            import openvdb as vdb  # type: ignore
        else:
            import pyopenvdb as vdb  # type: ignore

            is_pyopenvdb = True

        file_path = self.path_to_vdb(file, center=center, invert=invert)

        file_format = None
        if file.endswith((".map", ".map.gz", ".map.bz2")):
            file_format = "mrc"

        metadata = {
            "filepath": file,
            "invert": invert,
            "center": center,
        }
        gobj = self._parse_grid_with_fallback(
            file, file_format=file_format, metadata=metadata
        )
        if invert:
            gobj.grid = np.max(gobj.grid) - gobj.grid
        self.grid = gobj

        # If the grid has already been converted to a .vdb and overwrite is False, return that instead
        if os.path.exists(file_path) and not overwrite:
            # Also check that the file has the same invert and center settings
            vdb_grid = vdb.readAllGridMetadata(file_path)[0]
            if (
                "MN_invert" in vdb_grid
                and vdb_grid["MN_invert"] == invert
                and "MN_center" in vdb_grid
                and vdb_grid["MN_center"] == center
            ):
                return file_path

        grid = gobj.grid

        dataType = grid.dtype
        # enables different grid types
        if dataType == "int8" or dataType == "int16" or dataType == "int32":
            grid = grid.astype("int32")
            vdb_grid = vdb.Int32Grid()
        elif dataType == "int64":
            vdb_grid = vdb.Int64Grid()
        elif dataType == "float32":
            vdb_grid = vdb.FloatGrid()
        elif dataType == "float64":
            vdb_grid = vdb.DoubleGrid()
        else:
            grid = grid.astype(float)
            vdb_grid = vdb.FloatGrid()

        # The np.copy is needed to force numpy to actually rewrite the data in memory
        # since openvdb seems to read is straight from memory without checking the striding
        vdb_grid.copyFromArray(np.copy(grid, "C"))
        vdb_grid.gridClass = vdb.GridClass.FOG_VOLUME
        vdb_grid.name = "density"

        if center:
            offset = -np.array(grid.shape) * 0.5 * gobj.delta
        else:
            offset = np.array(gobj.origin)

        # apply transformations
        if is_pyopenvdb:
            vdb_grid.transform.scale(np.array(gobj.delta) * world_scale)
            vdb_grid.transform.translate(offset * world_scale)
        else:
            vdb_grid.transform.preScale(np.array(gobj.delta) * world_scale)
            vdb_grid.transform.postTranslate(offset * world_scale)

        # Set some metadata for the vdb file, so we can check if it's already
        # been converted correctly
        vdb_grid["MN_invert"] = invert
        vdb_grid["MN_center"] = center

        # Write the grid to a .vdb file
        if os.path.exists(file_path):
            os.remove(file_path)
        vdb.write(file_path, grids=[vdb_grid])
        # delete the vdb grid to avoid nanobind leak errors
        del vdb_grid
        # Return the path to the output vdb file
        return file_path
