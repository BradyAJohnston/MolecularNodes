"""
Benchmarks for parsing steps that the load benchmarks don't isolate: small-molecule
files, and biological assemblies on a structure with many copies (1f2n is an
icosahedral capsid), where `Molecule.load` would otherwise hide the cost.
"""

from pathlib import Path
import biotite.structure.io.pdbx as pdbx
import pytest
from molecularnodes.entities.molecule.pdb import PDBReader
from molecularnodes.entities.molecule.pdbx import CIFAssemblyParser
from molecularnodes.entities.molecule.reader import read_structure

DATA_DIR = Path(__file__).parent.parent / "tests" / "data"


@pytest.mark.benchmark(group="parse")
def test_read_sdf(benchmark):
    reader = benchmark(read_structure, DATA_DIR / "caffeine.sdf")
    assert reader.array.array_length() > 0


@pytest.mark.benchmark(group="assemblies")
def test_cif_assemblies(benchmark):
    cif = pdbx.CIFFile.read(DATA_DIR / "1f2n.cif")
    assemblies = benchmark(lambda: CIFAssemblyParser(cif).get_assemblies())
    assert len(assemblies) > 0


@pytest.mark.benchmark(group="assemblies")
def test_pdb_assemblies(benchmark):
    reader = PDBReader(DATA_DIR / "1f2n.pdb")
    assemblies = benchmark(reader._assemblies)
    assert len(assemblies) > 0
