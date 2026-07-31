import re
from pathlib import Path
from .cellpack import CellPack
from .cryosparc import CryoSPARC
from .star import StarFile


def load_starfile(file_path, node_setup=True, world_scale=0.01):
    ensemble = StarFile.from_starfile(file_path)
    ensemble.create_object(
        name=Path(file_path).name, node_setup=node_setup, world_scale=world_scale
    )

    return ensemble


def load_cryosparc(
    file_path: Path,
    node_setup=True,
    world_scale=0.01,
):
    ensemble = CryoSPARC(file_path=file_path)
    if search_results := re.search(r"(J\d+)_(.+)_exported.cs", str(file_path)):
        name = f"CryoSPARC {search_results.group(1)} {search_results.group(2)}"
    else:
        name = CryoSPARC.DEFAULT_NAME
    ensemble.create_object(name=name, node_setup=node_setup, world_scale=world_scale)

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
