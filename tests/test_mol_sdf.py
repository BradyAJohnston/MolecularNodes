import molecularnodes as mn
import molecularnodes.blender as bl
import pytest

from .constants import data_dir, attributes
from .utils import NumpySnapshotExtension


formats = ["mol", "sdf"]


@pytest.mark.parametrize("format", formats)
def test_open(snapshot_custom, format):
    molecule = mn.Molecule(data_dir / f"caffeine.{format}")

    assert molecule.atom_array


@pytest.mark.parametrize("format", formats)
@pytest.mark.parametrize("style", ["ball_and_stick", "spheres", "surface"])
def test_load(snapshot_custom: NumpySnapshotExtension, format, style):
    mol = mn.Molecule(data_dir / f"caffeine.{format}")
    mol.add_style(style=style)
    assert mol.object

    if style == "spheres":
        bl.nodes.get_style_node(mol.object).inputs[
            "Sphere As Mesh"
        ].default_value = True
    mn.blender.nodes.realize_instances(mol.object)

    for attribute in attributes:
        try:
            assert snapshot_custom == mol.named_attribute(attribute, evaluate=True)
        except AttributeError as e:
            assert e
