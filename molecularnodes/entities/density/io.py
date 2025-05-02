from pathlib import Path
from .mrc import MRC


def load(
    file_path: str | Path,
    name: str = "NewDensity",
    invert: bool = False,
    setup_nodes: bool = True,
    style: str = "density_surface",
    center: bool = False,
    overwrite: bool = False,
):
    density = MRC(
        file_path=file_path, center=center, invert=invert, overwrite=overwrite
    )
    density.create_object(
        name=Path(file_path).name, setup_nodes=setup_nodes, style=style
    )
    return density
