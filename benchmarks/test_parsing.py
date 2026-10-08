"""
Benchmarks for the pure-Python parsing layer that sits underneath
`mn.Molecule.load`: reading files into annotated biotite arrays, computing the
standard annotations, chain colors, biological assemblies and selections.
These isolate the parsing cost from building the Blender object and node tree.
"""

from pathlib import Path
import biotite.structure.io.pdbx as pdbx
import pytest
import molecularnodes as mn
from molecularnodes import color
from molecularnodes.entities.molecule.pdb import PDBReader
from molecularnodes.entities.molecule.pdbx import CIFAssemblyParser
from molecularnodes.entities.molecule.reader import ReaderBase, read_structure

DATA_DIR = Path(__file__).parent.parent / "tests" / "data"


@pytest.mark.benchmark(group="parse")
@pytest.mark.parametrize(
    "filename", ["1f2n.bcif", "1f2n.cif", "1f2n.pdb", "8U8W.bcif", "caffeine.sdf"]
)
def test_read_structure(benchmark, filename):
    reader = benchmark(read_structure, DATA_DIR / filename)
    assert reader.array.array_length() > 0


@pytest.mark.benchmark(group="parse")
def test_set_standard_annotations(benchmark):
    array = read_structure(DATA_DIR / "1f2n.bcif").get_structure()
    result = benchmark(ReaderBase.set_standard_annotations, array)
    assert "is_backbone" in result.get_annotation_categories()


@pytest.mark.benchmark(group="parse")
def test_color_chains(benchmark):
    array = read_structure(DATA_DIR / "1f2n.bcif").array[0]
    colors = benchmark(color.color_chains, array.atomic_number, array.chain_id)
    assert len(colors) == array.array_length()


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


@pytest.mark.benchmark(group="selection")
def test_selection_from_string(benchmark):
    mol = mn.Molecule.load(DATA_DIR / "1f2n.bcif")
    sel = benchmark(mol.selections.from_string, "around 5.0 resname CA")
    assert sel is not None
