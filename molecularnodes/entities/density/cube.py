"""
Reader for Gaussian cube files (``.cube`` / ``.cub``).

The format is described at https://paulbourke.net/dataformats/cube/. A cube file
holds a volumetric grid (orbitals, electron density, electrostatic potential) on
voxel axes that can be sheared, together with the atoms of the structure.
"""

import gzip
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np

BOHR_TO_ANGSTROM = 0.529177210903

CUBE_EXTENSIONS = (".cube", ".cub", ".cube.gz", ".cub.gz")


def is_cube(file: str | Path) -> bool:
    """Whether a path has a cube file extension."""
    return str(file).lower().endswith(CUBE_EXTENSIONS)


@dataclass
class CubeFile:
    """
    The contents of a cube file, in Å.

    Attributes
    ----------
    comments : tuple[str, str]
        The two comment lines of the header.
    origin : np.ndarray
        Position of the centre of voxel ``[0, 0, 0]``.
    axes : np.ndarray
        3x3 array whose rows are the step along each array axis of ``values``.
    values : np.ndarray
        The grid values, with shape ``(n_datasets, n1, n2, n3)``.
    atomic_numbers : np.ndarray
        Atomic number of each atom.
    charges : np.ndarray
        The nuclear charge column of each atom.
    positions : np.ndarray
        Positions of the atoms, shape ``(n_atoms, 3)``.
    dataset_ids : list[int]
        The orbital numbers of the datasets, when the file lists them (the number
        of atoms in the header is negative), otherwise empty.
    """

    comments: tuple[str, str]
    origin: np.ndarray
    axes: np.ndarray
    values: np.ndarray
    atomic_numbers: np.ndarray
    charges: np.ndarray
    positions: np.ndarray
    dataset_ids: list[int] = field(default_factory=list)


def count_datasets(file: str | Path) -> int:
    """
    The number of datasets (such as orbitals) in a cube file, read from its header.
    """
    opener = gzip.open if str(file).endswith(".gz") else open
    with opener(file, "rt") as f:
        f.readline()
        f.readline()
        n_atoms = int(f.readline().split()[0])
        if n_atoms >= 0:
            return 1
        for _ in range(3 + abs(n_atoms)):
            f.readline()
        return int(f.readline().split()[0])


def read_cube(file: str | Path) -> CubeFile:
    """
    Read a cube file.

    Bohr values are converted to Å. Following the format, a negative voxel count
    on the first axis means the file is already in Å.

    Parameters
    ----------
    file : str | Path
        Path to a ``.cube`` or ``.cub`` file, optionally gzipped.

    Returns
    -------
    CubeFile
        The grid and atoms of the file.
    """
    opener = gzip.open if str(file).endswith(".gz") else open
    with opener(file, "rt") as f:
        comments = (f.readline().rstrip("\n"), f.readline().rstrip("\n"))
        # Gaussian may add a fifth column (values per voxel) to this line
        header = f.readline().split()
        n_atoms = int(header[0])
        origin = np.array(header[1:4], dtype=float)
        counts = []
        axes = []
        for _ in range(3):
            line = f.readline().split()
            counts.append(int(line[0]))
            axes.append([float(x) for x in line[1:4]])
        atoms = np.array(
            [f.readline().split()[:5] for _ in range(abs(n_atoms))], dtype=float
        ).reshape(-1, 5)
        tokens = f.read().split()

    dataset_ids = []
    if n_atoms < 0:
        # orbital cube: the count of orbitals, then their numbers, before the data
        n_ids = int(tokens[0])
        dataset_ids = [int(x) for x in tokens[1 : 1 + n_ids]]
        tokens = tokens[1 + n_ids :]

    in_angstrom = counts[0] < 0
    scale = 1.0 if in_angstrom else BOHR_TO_ANGSTROM
    shape = tuple(abs(n) for n in counts)
    values = np.array(tokens, dtype=float)
    n_datasets = max(len(dataset_ids), 1)
    if values.size != np.prod(shape) * n_datasets:
        raise ValueError(
            f"Cube file '{file}' holds {values.size} values, expected "
            f"{np.prod(shape) * n_datasets} for a {shape} grid with "
            f"{n_datasets} dataset(s)"
        )
    # the last axis varies fastest, and datasets are interleaved per voxel
    values = np.moveaxis(values.reshape(*shape, n_datasets), -1, 0)

    return CubeFile(
        comments=comments,
        origin=origin * scale,
        axes=np.array(axes, dtype=float) * scale,
        values=values,
        atomic_numbers=atoms[:, 0].astype(int),
        charges=atoms[:, 1],
        positions=atoms[:, 2:5] * scale,
        dataset_ids=dataset_ids,
    )
