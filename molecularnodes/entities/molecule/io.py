from pathlib import Path
import io
from pathlib import Path
from ...blender import path_resolve

from ...download import download
from .base import Molecule
from .pdb import PDB
from .pdbx import BCIF, CIF
from .sdf import SDF


def parse(filepath) -> Molecule:
    """Parse a molecular structure file into a Molecule object.

    Parameters
    ----------
    filepath : str or io.BytesIO
        Path to the molecular structure file or BytesIO object containing file data

    Returns
    -------
    Molecule
        The parsed molecular structure object

    Raises
    ------
    ValueError
        If the file format is not supported
    InvalidFileError
        If the file cannot be parsed with the standard parser
    """
    # TODO: I don't like that we might be dealing with bytes or a filepath here,
    # I need to work out a nicer way to have it be cleanly one or the other

    if isinstance(filepath, io.BytesIO):
        suffix = ".bcif"
    else:
        filepath = path_resolve(filepath)
        suffix = Path(filepath).suffix

    parser = {
        ".pdb": PDB,
        ".pdbx": CIF,
        ".cif": CIF,
        ".bcif": BCIF,
        ".mol": SDF,
        ".sdf": SDF,
    }

    if suffix not in parser:
        raise ValueError(f"Unable to open local file. Format '{suffix}' not supported.")

    return parser[suffix](filepath)


def fetch(
    code: str,
    style: str | None = "spheres",
    centre: str | None = None,
    del_solvent: bool = True,
    del_hydrogen: bool = False,
    cache_dir: str | Path | None = None,
    build_assembly: bool = False,
    database: str = "rcsb",
    format: str = "bcif",
    color: str = "common",
) -> Molecule:
    """Fetch and create a molecular structure from online databases.

    Parameters
    ----------
    code : str
        The PDB code of the structure to fetch
    style : str or None, optional
        The visualization style to apply, by default "spheres"
    centre : str, optional
        Method for centering the structure, by default ""
    del_solvent : bool, optional
        Whether to delete solvent molecules, by default True
    del_hydrogen : bool, optional
        Whether to delete hydrogen atoms, by default False
    cache_dir : str or None, optional
        Directory to cache downloaded files, by default None
    build_assembly : bool, optional
        Whether to build the biological assembly, by default False
    database : str, optional
        Database to fetch from ("rcsb", "alphafold"), by default "rcsb"
    format : str, optional
        File format to download ("bcif", "pdb", etc), by default "bcif"
    color : str, optional
        Coloring scheme to apply, by default "common"

    Returns
    -------
    Molecule
        The created molecular structure object

    """
    if build_assembly:
        centre = ""

    file_path = download(code=code, format=format, cache=cache_dir, database=database)

    mol = parse(file_path)

    obj = mol.create_object(
        name=code,
        centre=centre,
        style=style,
        del_solvent=del_solvent,
        del_hydrogen=del_hydrogen,
        build_assembly=build_assembly,
        color=color,
    )

    obj.mn["code"] = code
    obj.mn["entity_type"] = format

    return mol


def load_local(
    file_path,
    centre: str | None = "",
    style: str | None = "spheres",
    del_solvent=True,
    del_hydrogen=False,
    build_assembly=False,
):
    mol = parse(file_path)
    mol.create_object(
        name=Path(file_path).stem,
        style=style,
        build_assembly=build_assembly,
        centre=centre,
        del_solvent=del_solvent,
        del_hydrogen=del_hydrogen,
    )
    return mol
