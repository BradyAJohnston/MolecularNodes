"""Minimal reader for small-molecule crystallographic ("core CIF") files.

Files from crystallographic databases (ICSD, COD) use the core CIF dictionary:
flat ``_atom_site_fract_*`` tags with a unit cell and symmetry operators,
rather than the mmCIF ``atom_site`` category with Cartesian coordinates that
biotite parses. This module implements just enough of the format to visualise
the crystal, as an MDAnalysis topology parser plus single-frame coordinate
reader: the asymmetric unit is expanded by the symmetry operators to fill one
unit cell, and fractional coordinates are converted to Cartesian angstroms.

Deliberately out of scope: disorder and partial-occupancy handling beyond
storing the occupancy value, anisotropic displacement parameters, and multiple
data blocks (tags from later blocks simply override earlier ones).
"""

import re
import shlex
from pathlib import Path
from typing import NamedTuple
import numpy as np
from MDAnalysis.coordinates.base import SingleFrameReaderBase
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import (
    Atomids,
    Atomnames,
    Elements,
    Occupancies,
    Resids,
    Resnames,
    Segids,
)
from MDAnalysis.lib.mdamath import triclinic_vectors
from MDAnalysis.topology.base import TopologyReaderBase

# positions are considered the same symmetry-equivalent atom at this many
# decimal places of the fractional coordinate
_FRACT_DECIMALS = 4

_CELL_TAGS = (
    "_cell_length_a",
    "_cell_length_b",
    "_cell_length_c",
    "_cell_angle_alpha",
    "_cell_angle_beta",
    "_cell_angle_gamma",
)
_SYMOP_TAGS = ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz")


def is_core_cif(file_path: str | Path) -> bool:
    """Cheaply sniff whether a .cif file is a core (small-molecule) CIF rather
    than an mmCIF, by which style of atom_site tag appears first."""
    with open(file_path) as file:
        for line in file:
            tag = line.strip().lower()
            if tag.startswith("_atom_site."):
                return False
            if tag.startswith("_atom_site_fract_"):
                return True
    return False


class Crystal(NamedTuple):
    names: np.ndarray
    elements: np.ndarray
    occupancies: np.ndarray
    positions: np.ndarray  # Cartesian angstroms, one expanded unit cell
    cell: np.ndarray  # [a, b, c, alpha, beta, gamma]


def _tokenize(text: str) -> list[str]:
    """CIF tokens: whitespace-separated, '...'/"..." quoted, '#' comments and
    semicolon-delimited multiline text fields."""
    tokens = []
    lines = iter(text.splitlines())
    for line in lines:
        if line.startswith(";"):
            block = [line[1:]]
            for continuation in lines:
                if continuation.startswith(";"):
                    break
                block.append(continuation)
            tokens.append("\n".join(block).strip())
            continue
        try:
            tokens.extend(shlex.split(line, comments=True))
        except ValueError:
            # unbalanced quote; take the line whole rather than failing
            tokens.extend(line.split("#")[0].split())
    return tokens


def _parse_tags(text: str) -> dict[str, list[str]]:
    "Map each CIF tag to its list of values (single values become one-item lists)."
    tokens = _tokenize(text)
    tags: dict[str, list[str]] = {}
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token.lower() == "loop_":
            columns = []
            i += 1
            while i < len(tokens) and tokens[i].startswith("_"):
                columns.append(tokens[i])
                i += 1
            values = []
            while i < len(tokens):
                value = tokens[i]
                if (
                    value.startswith("_")
                    or value.lower() in ("loop_",)
                    or value.lower().startswith("data_")
                ):
                    break
                values.append(value)
                i += 1
            for c, column in enumerate(columns):
                tags[column] = values[c :: len(columns)]
        elif token.startswith("_") and i + 1 < len(tokens):
            tags[token] = [tokens[i + 1]]
            i += 2
        else:
            i += 1
    return tags


def _float(value: str) -> float:
    "A CIF number, dropping the '(su)' uncertainty suffix; '.'/'?' become nan."
    value = re.sub(r"\(\d+\)$", "", value.strip())
    if value in (".", "?", ""):
        return np.nan
    return float(value)


def _element(symbol: str) -> str:
    "Element from a type symbol or site label, e.g. 'Ru4+' -> 'Ru', 'O3' -> 'O'."
    match = re.match(r"([A-Za-z]{1,2})", symbol.strip())
    if match is None:
        return ""
    return match.group(1).capitalize()


def _apply_symop(op: str, fract: np.ndarray) -> np.ndarray:
    "Apply one 'x+1/2, -y, z' style operator to an (n, 3) fractional array."
    x, y, z = fract[:, 0], fract[:, 1], fract[:, 2]
    components = op.lower().split(",")
    if len(components) != 3 or not set(op.lower()) <= set("xyz0123456789+-*/., "):
        raise ValueError(f"unsupported symmetry operator: {op!r}")
    return np.stack(
        [eval(c, {"__builtins__": {}}, {"x": x, "y": y, "z": z}) for c in components],
        axis=1,
    )


def parse_core_cif(file_path: str | Path) -> Crystal:
    "Parse a core CIF file and symmetry-expand its sites to one unit cell."
    tags = _parse_tags(Path(file_path).read_text())

    try:
        cell = np.array([_float(tags[tag][0]) for tag in _CELL_TAGS], float)
        fract = np.array(
            [[_float(v) for v in tags[f"_atom_site_fract_{axis}"]] for axis in "xyz"],
            float,
        ).T
    except KeyError as e:
        raise ValueError(f"core CIF file is missing the {e} tag") from None

    labels = tags.get("_atom_site_label", [f"X{i}" for i in range(len(fract))])
    symbols = tags.get("_atom_site_type_symbol", labels)
    occupancies = [_float(v) for v in tags.get("_atom_site_occupancy", [])]
    occupancies += [1.0] * (len(fract) - len(occupancies))

    symops = ["x, y, z"]
    for tag in _SYMOP_TAGS:
        if tag in tags:
            symops = tags[tag]
            break

    # expand each site by every operator, fold into the unit cell and keep the
    # unique positions per site so special positions are not duplicated
    names, elements, occs, fractions = [], [], [], []
    for site, label, symbol, occupancy in zip(fract, labels, symbols, occupancies):
        images = np.concatenate([_apply_symop(op, site[None, :]) for op in symops])
        images = np.round(images % 1.0, _FRACT_DECIMALS) % 1.0
        images = np.unique(images, axis=0)
        fractions.append(images)
        names.extend([label] * len(images))
        elements.extend([_element(symbol)] * len(images))
        occs.extend([1.0 if np.isnan(occupancy) else occupancy] * len(images))

    box = triclinic_vectors(cell.astype(np.float32))
    positions = (np.concatenate(fractions) @ box).astype(np.float32)

    return Crystal(
        names=np.array(names),
        elements=np.array(elements),
        occupancies=np.array(occs, float),
        positions=positions,
        cell=cell,
    )


class CoreCIFParser(TopologyReaderBase):
    "MDAnalysis topology parser for small-molecule crystallographic CIF files."

    format = "CORECIF"

    def parse(self, **kwargs) -> Topology:
        crystal = parse_core_cif(self.filename)
        n_atoms = len(crystal.positions)
        attrs = [
            Atomids(np.arange(n_atoms) + 1),
            Atomnames(crystal.names),
            Elements(crystal.elements),
            Occupancies(crystal.occupancies),
            Resids(np.array([1])),
            Resnames(np.array(["XTL"])),
            Segids(np.array(["SYSTEM"])),
        ]
        return Topology(
            n_atoms=n_atoms,
            n_res=1,
            n_seg=1,
            attrs=attrs,
            atom_resindex=np.zeros(n_atoms, int),
            residue_segindex=np.zeros(1, int),
        )


class CoreCIFReader(SingleFrameReaderBase):
    "MDAnalysis coordinate reader for small-molecule crystallographic CIF files."

    format = "CORECIF"
    units = {"time": None, "length": "Angstrom"}

    def _read_first_frame(self) -> None:
        crystal = parse_core_cif(self.filename)
        self.n_atoms = len(crystal.positions)
        self.ts = self._Timestep.from_coordinates(crystal.positions, **self._ts_kwargs)
        self.ts.dimensions = crystal.cell
        self.ts.frame = 0
