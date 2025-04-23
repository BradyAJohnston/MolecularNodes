import pytest
import molecularnodes as mn
from molecularnodes.nodes import nodes as nodes
from .constants import attributes, data_dir
from .utils import NumpySnapshotExtension

formats = ["mol", "sdf"]


@pytest.mark.parametrize("format", formats)
def test_open(snapshot_custom, format):
    molecule = mn.Molecule.load(data_dir / f"caffeine.{format}")

    assert molecule.array


@pytest.mark.parametrize("format", formats)
@pytest.mark.parametrize("style", ["ball_and_stick", "spheres", "surface"])
def test_load(snapshot_custom: NumpySnapshotExtension, format, style):
    mol = mn.Molecule.load(data_dir / f"caffeine.{format}").add_style(style=style)
    assert mol.object

    if style != "surface":
        mol.styles[0].sphere_geometry = "Mesh"

    for attribute in attributes:
        try:
            print(f"{attribute=}")
            assert snapshot_custom == mol.named_attribute(attribute, evaluate=True)
        except AttributeError as e:
            assert e
