import pytest
import molecularnodes as mn
from molecularnodes.nodes import nodes as nodes
from molecularnodes.nodes.geometry import StyleBallAndStick, StyleSpheres, StyleSurface
from .constants import data_dir
from .utils import GeometrySet

formats = ["mol", "sdf"]


@pytest.mark.parametrize("format", formats)
def test_open(format):
    molecule = mn.Molecule.load(data_dir / f"caffeine.{format}")

    assert molecule.array


@pytest.mark.parametrize("format", formats)
@pytest.mark.parametrize(
    "style",
    ["ball_and_stick", "spheres", "surface"],
)
def test_load(snapshot, format, style):
    mol = mn.Molecule.load(data_dir / f"caffeine.{format}")
    assert mol.object
    with mol.tree.reset() as tree:
        (
            tree.atoms
            >> {
                "ball_and_stick": StyleBallAndStick(sphere_geometry="Mesh"),
                "spheres": StyleSpheres(geometry="Mesh"),
                "surface": StyleSurface(),
            }[style]
            >> tree.join
        )

    assert snapshot == GeometrySet(mol.object).summary()
