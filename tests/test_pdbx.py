import molecularnodes as mn
import numpy as np

import random
from .constants import data_dir
from .utils import NumpySnapshotExtension, sample_attribute


def test_ss_label_to_int():
    examples = ["TURN_TY1_P68", "BEND64", "HELX_LH_PP_P9", "STRN44"]
    assert [3, 3, 1, 2] == [mn.io.parse.cif._ss_label_to_int(x) for x in examples]


def test_get_ss_from_mmcif(snapshot_custom: NumpySnapshotExtension):
    mol = mn.io.load(data_dir / "1cd3.cif")

    # mol2, fil2 = mn.io.fetch('1cd3')

    random.seed(6)
    random_idx = random.sample(range(len(mol)), 100)

    # assert (mol.sec_struct == mol2.sec_struct)[random_idx].all()

    assert snapshot_custom == mol.array.sec_struct[random_idx]


def test_secondary_structure_no_helix(snapshot_custom):
    m = mn.io.fetch("7ZL4", cache_dir=data_dir)

    assert snapshot_custom == sample_attribute(
        m.object, "sec_struct", n=500, evaluate=False
    )
