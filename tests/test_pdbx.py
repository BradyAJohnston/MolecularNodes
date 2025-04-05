import random
import molecularnodes as mn
from .constants import data_dir
from .utils import NumpySnapshotExtension


def test_ss_label_to_int():
    examples = ["TURN_TY1_P68", "BEND64", "HELX_LH_PP_P9", "STRN44"]
    assert [3, 3, 1, 2] == [
        mn.entities.molecule.pdbx._ss_label_to_int(x) for x in examples
    ]


def test_entity_parsing():
    mn.Molecule.fetch("6VBU", format="bcif")
    assert True


def test_get_ss_from_mmcif(snapshot_custom: NumpySnapshotExtension):
    mol = mn.Molecule.load(data_dir / "1cd3.cif")
    random.seed(6)
    random_idx = random.sample(range(len(mol)), 100)
    assert snapshot_custom == mol.array.sec_struct[random_idx]


def test_secondary_structure_no_helix(snapshot_custom):
    m = mn.Molecule.fetch("7ZL4", cache=data_dir)
    assert snapshot_custom == m.named_attribute("sec_struct")
