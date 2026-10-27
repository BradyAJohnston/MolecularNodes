import gzip
import os
from pathlib import Path
from typing import Sequence
import bpy
import databpy
import numpy as np
from gridData import Grid
from ... import blender as bl
from ...blender import coll
from ...nodes.geometry import (
    DensityStyleISOSurface,
    DensityStyleSurface,
    DensityStyleWire,
    ImportVolume,
)
from .base import Density
from .cube import count_datasets, is_cube, read_cube


def cell_vectors(lengths, angles) -> np.ndarray:
    """
    Unit cell vectors from cell lengths and angles.

    Uses the PDB / CCP4 convention: ``a`` along x, ``b`` in the xy plane.

    Parameters
    ----------
    lengths : array-like
        Cell lengths ``(a, b, c)``.
    angles : array-like
        Cell angles ``(alpha, beta, gamma)`` in degrees.

    Returns
    -------
    np.ndarray
        3x3 array whose rows are the ``a``, ``b`` and ``c`` vectors.
    """
    a, b, c = (float(x) for x in lengths)
    alpha, beta, gamma = np.radians(np.asarray(angles, dtype=float))
    cy = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(max(1.0 - np.cos(beta) ** 2 - cy**2, 0.0))
    vectors = np.array(
        [
            [a, 0.0, 0.0],
            [b * np.cos(gamma), b * np.sin(gamma), 0.0],
            [c * np.cos(beta), c * cy, c * cz],
        ]
    )
    # cos(90 deg) is ~6e-17, not 0; zero it so right angles stay exact
    vectors[np.abs(vectors) < 1e-12 * np.abs(vectors).max()] = 0.0
    return vectors


def _escape_format(path: str) -> str:
    """Escape braces in a path, for the `Format String` node in `Import Volume`."""
    return path.replace("{", "{{").replace("}", "}}")


def _dx_delta(file: str) -> np.ndarray | None:
    """
    Read the full ``delta`` matrix from the header of an OpenDX file.

    GridDataFormats keeps only the diagonal of ``delta``, which drops the shear of
    non-orthogonal grids. Returns ``None`` if the header does not hold three
    ``delta`` lines.
    """
    opener = gzip.open if file.endswith(".gz") else open
    deltas = []
    with opener(file, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith("delta"):
                deltas.append([float(x) for x in line.split()[1:]])
            elif "data follows" in line:
                break
    if len(deltas) != 3 or any(len(d) != 3 for d in deltas):
        return None
    return np.array(deltas, dtype=float)


class Grids(Density):
    """
    A class for parsing grid files in the various formats:
    See: https://www.mdanalysis.org/GridDataFormats/gridData/formats.html#supported-file-formats
    Plain: `.dx`, `.plt`, `.ccp4`, `.mrc`, `.map`, `.pickle`
    Compressed: `.dx.gz`, `.ccp4.gz`, `.ccp4.bz2`, `.mrc.gz`, `.mrc.bz2`, `.map.gz`, `.map.bz2`

    Gaussian cube files (`.cube`, `.cub`, optionally gzipped) are read by
    :func:`~molecularnodes.entities.density.cube.read_cube`, with Bohr converted to Å.
    Each dataset of a cube file with several orbitals becomes one numbered `.vdb`
    file, which the `Index` of `Import Volume` picks between.

    It utilises `GridDataFormats` for file parsing, which is then converted into `openvdb` grids,
    that are written as `.vdb` files. The density's object is an empty mesh whose node tree
    imports them with `Import VDB`, so several files can load into one density and play
    as an animation.

    Float64 data, which includes every `.dx` file, is written as a 32-bit float grid:
    Blender loads 64-bit grids as `DOUBLE` grids, which its grid nodes (Get Named Grid,
    Sample Grid) cannot read. Values therefore carry float32 precision in Blender.
    """

    def __init__(
        self,
        file_path: str | Path | Sequence[str | Path],
        center: bool = False,
        invert: bool = False,
        overwrite: bool = False,
    ):
        paths = [file_path] if isinstance(file_path, (str, Path)) else list(file_path)
        if not paths:
            raise ValueError("No grid files given")
        super().__init__(file_path=paths[0])
        self.file_paths: list[Path] = [bl.path_resolve(p) for p in paths]
        self.grid: Grid = None  # GridDataFormats Grid object, of the first file
        # one grid per file, or per dataset of a cube file with several orbitals
        jobs: list[tuple[Path, int | None]] = []
        for path in self.file_paths:
            n = count_datasets(path) if is_cube(path) else 1
            jobs += [(path, i if n > 1 else None) for i in range(n)]
        if len(jobs) == 1:
            targets = [None]
        else:
            # several grids are written as numbered .vdb files next to the first
            # file, which `Import Volume` picks between by index
            base = self.path_to_vdb(
                str(self.file_paths[0]), center=center, invert=invert
            ).removesuffix(".vdb")
            targets = [f"{base}_{i}.vdb" for i in range(len(jobs))]
        self._cubes = {}  # parsed cube files, reused for each of their datasets
        self.files_vdb: list[str] = []
        for (path, dataset), target in reversed(list(zip(jobs, targets))):
            # converted last to first, so `grid` and `axes` describe the first grid
            self.files_vdb.insert(
                0,
                self.grid_to_vdb(
                    str(path),
                    center=center,
                    invert=invert,
                    overwrite=overwrite,
                    file_vdb=target,
                    dataset=dataset,
                ),
            )
        del self._cubes
        self.file_vdb = self.files_vdb[0]
        #: path for `Import Volume`, with `{index}` in place of the number when
        #: there are several files
        self.vdb_path = (
            _escape_format(self.file_vdb)
            if len(self.files_vdb) == 1
            else _escape_format(base) + "_{index}.vdb"
        )

    @classmethod
    def load(
        cls,
        file_path: str | Path | Sequence[str | Path],
        name: str | None = None,
        style: str | None = "density_surface",
        invert: bool = False,
        center: bool = False,
        overwrite: bool = False,
        follow_frame: bool | None = None,
        frame_start: int | None = None,
    ) -> "Grids":
        """Load density grid files and create their density object.

        Parameters
        ----------
        file_path : str | Path | Sequence[str | Path]
            Path to the grid file (see the class docstring for supported formats), or
            a list of paths to load as one density, such as the frames of an animation.
            The `Index` input of the tree's `Import Density` node picks the grid shown.
        name : str | None, optional
            Name for the created object. If None, the name of the first file is used.
        style : str | None, optional
            The density style to apply, by default "density_surface". If None, no
            style is added.
        invert : bool, optional
            Whether to invert the density values, by default False.
        center : bool, optional
            Whether to center the density at the origin, by default False.
        overwrite : bool, optional
            Whether to overwrite an existing converted ``.vdb`` file, by default False.
        follow_frame : bool | None, optional
            Whether the scene frame picks the grid, playing the files as an animation
            with one file per frame. Defaults to True when several files are given.
        frame_start : int | None, optional
            Scene frame that shows the first file when following the frame. Defaults
            to the scene's start frame.

        Returns
        -------
        Grids
            The loaded density entity.

        Notes
        -----
        Float64 grids, including every ``.dx`` file, are stored as 32-bit floats,
        because Blender's grid nodes (Get Named Grid, Sample Grid) return 0 for
        64-bit ``DOUBLE`` grids.
        """
        density = cls(
            file_path=file_path, center=center, invert=invert, overwrite=overwrite
        )
        density.create_object(
            name=name or density.file_paths[0].name,
            style=style,
            # several files play as an animation; the datasets of one file don't
            follow_frame=len(density.file_paths) > 1
            if follow_frame is None
            else follow_frame,
            frame_start=frame_start,
        )
        # record the source so the entity can be reloaded into a fresh session
        density.object.mn.filepath = str(density.file_paths[0])
        return density

    def _parse_grid_with_fallback(
        self,
        file: str,
        file_format: str | None,
        metadata: dict,
        dataset: int | None = None,
    ):
        """
        Try parsing the grid using GridDataFormats first; on failure, fall back to
        reading with the mrcfile package and constructing a minimal grid-like object.

        Returns an object with attributes: grid (np.ndarray), delta (np.ndarray),
        origin (np.ndarray), metadata (dict).
        """
        if is_cube(file):
            return self._parse_cube(file, metadata, dataset=dataset or 0)

        # First attempt: GridDataFormats
        try:
            return Grid(file, file_format=file_format, metadata=metadata)
        except Exception:
            pass

        # Fallback: try using mrcfile for MRC/CCP4/MAP formats
        try:
            import mrcfile

            with mrcfile.open(file, permissive=True) as mrc:
                h = mrc.header
                angles = [
                    float(h.cellb.alpha),
                    float(h.cellb.beta),
                    float(h.cellb.gamma),
                ]
                if not np.allclose(angles, 90.0):
                    return self._parse_mrc_non_orthogonal(mrc, angles, metadata)

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
                # We also swap the X and Z axes of the numpy array when constructing the Grid below.
                delta = np.array([vx, vy, vz], dtype=float)
                origin = np.array([ox, oy, oz], dtype=float)

                return Grid(
                    grid=np.swapaxes(mrc.data, 0, 2),
                    origin=origin,
                    delta=delta,
                    metadata=metadata,
                )
        except Exception as e:
            # Re-raise the original parsing error with context from fallback
            raise RuntimeError(
                f"Failed to parse grid file '{file}' with GridDataFormats and mrcfile fallback: {e}"
            )

    @staticmethod
    def _parse_mrc_non_orthogonal(mrc, angles, metadata: dict) -> Grid:
        """
        Build a grid from an MRC/CCP4 map whose unit cell is not orthogonal.

        GridDataFormats rejects these maps. The voxel axes are the cell vectors
        divided by the sampling, and are stored on the returned grid as
        ``mn_axes`` (rows are the step in Å along each array axis).
        """
        h = mrc.header
        # mrcfile data is indexed [section, row, column]; reorder it to crystal
        # x, y, z using mapc / mapr / maps, as GridDataFormats does
        axes_order = np.hstack([h.mapc, h.mapr, h.maps])
        axes_c_order = np.argsort(axes_order)
        transpose_order = np.argsort(axes_order[::-1])
        data = np.transpose(mrc.data, axes=transpose_order)

        sampling = np.array([h.mx, h.my, h.mz], dtype=float)
        sampling[sampling == 0] = np.array(data.shape, dtype=float)[sampling == 0]
        axes = (
            cell_vectors([h.cella.x, h.cella.y, h.cella.z], angles) / sampling[:, None]
        )
        starts = np.hstack([h.nxstart, h.nystart, h.nzstart])[axes_c_order]
        origin = (
            np.array([h.origin.x, h.origin.y, h.origin.z], dtype=float)
            + starts.astype(float) @ axes
        )

        grid = Grid(
            grid=data,
            origin=origin,
            delta=np.linalg.norm(axes, axis=1),
            metadata=metadata,
        )
        grid.mn_axes = axes
        return grid

    def _parse_cube(self, file: str, metadata: dict, dataset: int = 0) -> Grid:
        """
        Build a grid from one dataset of a Gaussian cube file.

        The voxel axes, which may be sheared, are stored on the returned grid as
        ``mn_axes`` (rows are the step in Å along each array axis).
        """
        cubes = getattr(self, "_cubes", {})
        if file not in cubes:
            cubes[file] = read_cube(file)
        cube = cubes[file]
        grid = Grid(
            grid=cube.values[dataset],
            origin=cube.origin,
            delta=np.linalg.norm(cube.axes, axis=1),
            metadata=metadata,
        )
        grid.mn_axes = cube.axes
        return grid

    def _voxel_axes(self, file: str, gobj) -> np.ndarray:
        """
        The step in Å along each array axis of the grid, as rows of a 3x3 array.

        Diagonal for orthogonal grids. Sheared or rotated grids come from the full
        ``delta`` matrix of a ``.dx`` file or the cell angles of an MRC/CCP4 map.
        """
        axes = getattr(gobj, "mn_axes", None)
        if axes is not None:
            return np.asarray(axes, dtype=float)
        if file.endswith((".dx", ".dx.gz")):
            delta = _dx_delta(file)
            if delta is not None:
                return delta
        return np.diag(np.asarray(gobj.delta, dtype=float))

    def create_object(
        self,
        name="NewDensity",
        style="density_surface",
        follow_frame: bool = False,
        frame_start: int | None = None,
    ) -> bpy.types.Object:
        """
        Create the density's object: an empty mesh whose node tree imports the grids.

        The tree starts with an `Import Volume` node that reads the converted
        ``.vdb`` files, followed by the style.

        Parameters
        ----------
        name : str, optional
            If not empty, renames the object with the new name. Default is "NewDensity".
        style : str, optional
            The style of the density object. Default is "density_surface". If None,
            the imported volume is shown unstyled.
        follow_frame : bool, optional
            Whether the scene frame picks the grid, by default False.
        frame_start : int | None, optional
            Scene frame that shows the first grid when following the frame. Defaults
            to the scene's start frame.

        Returns
        -------
        bpy.types.Object
            The density's object.
        """
        # an empty mesh at the world origin, so the grids line up with structures
        self.object = databpy.create_object(
            name=name or "NewDensity", collection=coll.mn()
        )
        self.object.mn.entity_type = self._entity_type.value

        if name and name != "":
            self.name = name

        threshold = float(np.quantile(self.grid.grid, 0.995))
        with self.tree.reset() as (_, join):
            (
                self.import_volume(follow_frame=follow_frame, frame_start=frame_start)
                >> (
                    None
                    if style is None
                    else {
                        "density_iso_surface": DensityStyleISOSurface,
                        "density_surface": DensityStyleSurface,
                        "density_wire": DensityStyleWire,
                    }[style](threshold=threshold)
                )
                >> join
            )

        return self.object

    def import_volume(
        self,
        name: str = "density",
        index: int = 0,
        follow_frame: bool = False,
        frame_start: int | None = None,
    ) -> ImportVolume:
        """
        An `Import Volume` node that reads this density's ``.vdb`` files.

        Start each branch of the density's tree with one, inside the tree:

        >>> with density.tree.reset() as (_, join):
        ...     density.import_volume() >> DensityStyleSurface() >> join

        Parameters
        ----------
        name : str, optional
            Name the grid is kept under for `Color Sample Grid`, by default "density".
            Give each branch that imports a different grid its own name.
        index : int, optional
            Which of several files to import, when not following the frame.
        follow_frame : bool, optional
            Whether the scene frame picks the file, to play the files as an animation.
        frame_start : int | None, optional
            Scene frame that shows the first file when following the frame. Defaults
            to the scene's start frame.

        Returns
        -------
        ImportVolume
            The node, added to the tree being built.
        """
        return ImportVolume(
            path=self.vdb_path,
            index=index,
            follow_frame=follow_frame,
            frame_start=bpy.context.scene.frame_start
            if frame_start is None
            else frame_start,
            count=len(self.files_vdb),
            name=name,
        )

    def grid_to_vdb(
        self,
        file: str,
        invert: bool = False,
        world_scale=0.1,
        center: bool = False,
        overwrite=False,
        file_vdb: str | None = None,
        dataset: int | None = None,
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
        file_vdb : str | None, optional
            Path to write the .vdb file to. Defaults to one next to `file`, named after it.
        dataset : int | None, optional
            For a cube file with several datasets (orbitals), the one to convert.

        Returns
        -------
        str
            The path to the converted .vdb file.
        """

        bpy.utils.expose_bundled_modules()
        import openvdb as vdb  # type: ignore

        file_path = file_vdb or self.path_to_vdb(file, center=center, invert=invert)

        file_format = None
        if file.endswith((".map", ".map.gz", ".map.bz2")):
            file_format = "mrc"

        metadata = {
            "filepath": file,
            "invert": invert,
            "center": center,
        }
        gobj = self._parse_grid_with_fallback(
            file, file_format=file_format, metadata=metadata, dataset=dataset
        )
        source = file if dataset is None else f"{file}#{dataset}"
        if invert:
            gobj.grid = np.max(gobj.grid) - gobj.grid
        self.grid = gobj
        axes = self._voxel_axes(file, gobj)
        self.axes = axes
        orthogonal = bool(np.array_equal(axes, np.diag(np.diag(axes))))

        # If the grid has already been converted to a .vdb and overwrite is False, return that instead
        if os.path.exists(file_path) and not overwrite:
            # Also check that the file has the same invert and center settings
            vdb_grid = vdb.readAllGridMetadata(file_path)[0]
            if (
                "MN_invert" in vdb_grid
                and vdb_grid["MN_invert"] == invert
                and "MN_center" in vdb_grid
                and vdb_grid["MN_center"] == center
                # sheared grids converted before the affine transform was written
                # have a wrong transform, so only reuse ones flagged as affine
                and (orthogonal or ("MN_affine" in vdb_grid and vdb_grid["MN_affine"]))
                and vdb_grid.valueTypeName != "double"
                # numbered files of a sequence are reused only for the same source
                and ("MN_source" not in vdb_grid or vdb_grid["MN_source"] == source)
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
            # Blender loads a float64 VDB grid as a DOUBLE grid, and Get Named Grid /
            # Sample Grid return 0 for it, so store float64 data (every .dx file) as
            # float32 to keep the grid usable in geometry nodes
            grid = grid.astype("float32")
            vdb_grid = vdb.FloatGrid()
        else:
            grid = grid.astype(float)
            vdb_grid = vdb.FloatGrid()

        # The np.copy is needed to force numpy to actually rewrite the data in memory
        # since openvdb seems to read is straight from memory without checking the striding
        vdb_grid.copyFromArray(np.copy(grid, "C"))
        vdb_grid.gridClass = vdb.GridClass.FOG_VOLUME
        vdb_grid.name = "density"

        if orthogonal:
            # orthogonal grid: scale and translate
            if center:
                offset = -np.array(grid.shape) * 0.5 * gobj.delta
            else:
                offset = np.array(gobj.origin)

            vdb_grid.transform.preScale(np.array(gobj.delta) * world_scale)
            vdb_grid.transform.postTranslate(offset * world_scale)
        else:
            # sheared or rotated grid: full affine transform. openvdb takes a
            # row-vector matrix, rows 0-2 are the voxel axes and row 3 the offset
            if center:
                offset = -(np.array(grid.shape) * 0.5) @ axes
            else:
                offset = np.array(gobj.origin)
            matrix = np.eye(4)
            matrix[:3, :3] = axes * world_scale
            matrix[3, :3] = offset * world_scale
            vdb_grid.transform = vdb.createLinearTransform(matrix.tolist())

        # Set some metadata for the vdb file, so we can check if it's already
        # been converted correctly
        vdb_grid["MN_invert"] = invert
        vdb_grid["MN_center"] = center
        vdb_grid["MN_source"] = source
        if not orthogonal:
            vdb_grid["MN_affine"] = True

        # Write the grid to a .vdb file
        if os.path.exists(file_path):
            os.remove(file_path)
        vdb.write(file_path, grids=[vdb_grid])
        # delete the vdb grid to avoid nanobind leak errors
        del vdb_grid
        # Return the path to the output vdb file
        return file_path
