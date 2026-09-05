from pathlib import Path
from .cellpack import CellPack
from .star import StarFile


def load_starfile(file_path, node_setup=True, world_scale=0.01):
    ensemble = StarFile.from_starfile(file_path)
    ensemble.create_object(
        name=Path(file_path).name, node_setup=node_setup, world_scale=world_scale
    )

    return ensemble


def load_cellpack(
    file_path,
    name="NewCellPackModel",
    node_setup=True,
    world_scale=0.01,
    fraction: float = 1,
):
    ensemble = CellPack(file_path)
    ensemble.create_object(
        name=name, node_setup=node_setup, world_scale=world_scale, fraction=fraction
    )

    return ensemble
