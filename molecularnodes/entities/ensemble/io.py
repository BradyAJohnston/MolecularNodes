from pathlib import Path
from .cellpack import CellPack
from .star import StarFile


def load_starfile(
    file_path: str | Path, node_setup: bool = True, world_scale: float = 0.1
) -> StarFile:
    ensemble = StarFile.from_starfile(file_path)
    ensemble.create_object(
        name=Path(file_path).name, node_setup=node_setup, world_scale=world_scale
    )

    return ensemble


def load_cellpack(
    file_path: str | Path,
    name: str | None = None,
    node_setup: bool = True,
    world_scale: float = 0.01,
    fraction: float = 1.0,
) -> CellPack:
    ensemble = CellPack(file_path)
    ensemble.create_object(
        name=name or Path(file_path).name,
        node_setup=node_setup,
        world_scale=world_scale,
        fraction=fraction,
    )

    return ensemble
