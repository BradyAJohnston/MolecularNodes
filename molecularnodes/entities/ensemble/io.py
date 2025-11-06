from pathlib import Path
from .cellpack import CellPack
from .star import StarFile


def load_starfile(file_path, node_setup=True):
    ensemble = StarFile.from_starfile(file_path)
    ensemble.create_object(
        name=Path(file_path).name, node_setup=node_setup
    )

    return ensemble


def load_cellpack(
    file_path,
    name="NewCellPackModel",
    node_setup=True,
    fraction: float = 1,
):
    ensemble = CellPack(file_path)
    ensemble.create_object(
        name=name, node_setup=node_setup, fraction=fraction
    )

    return ensemble
