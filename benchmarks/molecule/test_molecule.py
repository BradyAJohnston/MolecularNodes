import pytest
import molecularnodes as mn
from tests.conftest import DATA_DIR


# @pytest.mark.benchmark
@pytest.mark.parametrize("format", ["bcif", "cif", "pdb"])
def test_molecule_load(benchmark, format):
    """
    Benchmark the molecule loading process for different formats.
    """

    # Load the molecule using the specified format
    def load_molecule():
        return mn.Molecule.fetch("4ozs", format=format)

    benchmark(load_molecule)
