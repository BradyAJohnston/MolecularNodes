import numpy as np

from .density import Density
from ...blender import obj, nodes


class OpenDX(Density):
    def __init__(self, file_path):
        self.file_path: str = file_path
        self.data = self._read()
        self.file_vdb = self._grid_to_vdb(self.data)

    def _read(self):
        from gridData import Grid

        g = Grid(self.file_path)

        return g

    def _grid_to_vdb(self, data, world_scale=0.01):
        import pyopenvdb as vdb

        volume_grid = vdb.FloatGrid()
        volume_grid.copyFromArray(data.grid)
        volume_grid.gridClass = vdb.GridClass.FOG_VOLUME
        volume_grid.name = 'density'
        volume_grid.transform.translate(data.origin * world_scale)
        volume_grid.transform.scale(np.array(data.delta) * world_scale)

        self.file_vdb = self.path_to_vdb(self.file_path)

        vdb.write(self.file_vdb, grids=[volume_grid])

        return self.file_vdb

    def create_model(
        self,
        name='NewDensity',
        style='density_surface',
        setup_nodes=True,
        world_scale=0.01,
    ):
        self.object = obj.import_vdb(self.file_vdb)

        if setup_nodes:
            nodes.create_starting_nodes_density(
                self.object, style=style, threshold=1)

        return self.object
