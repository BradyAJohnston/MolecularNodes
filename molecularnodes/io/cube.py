import numpy as np
from molecularnodes.io.parse.density import Density
from pathlib import Path

fl = "tests/data/cube/frame_0.cube"


class CubeFile:
    def __init__(self, file: str) -> None:
        self.file = file

    def _read(self):
        self


class CubeDensity(Density):
    def __init__(self, file_path: str) -> None:
        super().__init__(file_path)
        self.array = self._cube_to_numpy(self.file_path)
        self.vdb_path = self._numpy_to_vdb(self.array)

    def _cube_to_numpy(self, path: Path) -> np.ndarray:
        with open(path) as f:
            lines = f.readlines()[23:]
        array = np.genfromtxt("\n".join(lines).split()).reshape((80, 80, 80))

        return array

    def _numpy_to_vdb(self, array: np.ndarray) -> str:
        import pyopenvdb as vdb

        grid = vdb.FloatGrid()

        grid.copyFromArray(array)
        self.threshold = np.quantile(array, 0.995)
        grid.gridClass = vdb.GridClass.FOG_VOLUME
        grid.name = "density"

        vdb_file_path = self.path_to_vdb(self.file_path, center=False, invert=False)

        world_scale = 0.01
        bohr_to_angstrom = 0.529177249
        voxel_size = np.array((0.226767, 0.309443, 0.181886))
        flip_array = np.array((1, 1, 1))
        origin_point = np.array((85.415621, 82.203086, 87.210861))

        # scale the grid
        grid.transform.scale(flip_array * world_scale * voxel_size * bohr_to_angstrom)
        grid.transform.translate(world_scale * origin_point * bohr_to_angstrom)

        vdb.write(vdb_file_path, grids=[grid])
        return vdb_file_path


if __name__ == "__main__":
    cube = CubeDensity(fl)
    print(cube.array)
    print(cube.vdb_path)
