import pytest
import molecularnodes as mn
from molecularnodes.nodes import nodes as nodes
from molecularnodes.nodes.assets import StyleBallAndStick, StyleSpheres, StyleSurface
from .constants import attributes, data_dir
from .utils import NumpySnapshotExtension

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
def test_load(snapshot_custom: NumpySnapshotExtension, format, style):
    mol = mn.Molecule.load(data_dir / f"caffeine.{format}")
    assert mol.object
    with mol.tree as tree:
        atoms, join = tree.reset()

        (
            atoms
            >> {
                "ball_and_stick": StyleBallAndStick(sphere_geometry="Mesh"),
                "spheres": StyleSpheres(geometry="Mesh"),
                "surface": StyleSurface(),
            }[style]
            >> join
        )

    for attribute in attributes:
        try:
            print(f"{attribute=}")
            assert snapshot_custom == mol.named_attribute(attribute, evaluate=True)
        except AttributeError as e:
            assert e
