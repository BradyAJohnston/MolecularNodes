import databpy as db
import molecularnodes as mn
from .utils import GeometrySet


def test_get_set(snapshot):
    code = "1BNA"
    mol = mn.Molecule.fetch(code).add_style("ribbon")
    geom = GeometrySet(mol.object)
    assert snapshot == geom
