from pathlib import Path
from .grids import Grids


def load(
    file_path: str | Path,
    name: str = "NewDensity",
    invert: bool = False,
    style: str = "density_surface",
    center: bool = False,
    overwrite: bool = False,
):
    density = Grids(
        file_path=file_path, center=center, invert=invert, overwrite=overwrite
    )
    density.create_object(name=Path(file_path).name, style=style)
    # record the source so the entity can be reloaded into a fresh session
    density.object.mn.filepath = str(file_path)
    return density
