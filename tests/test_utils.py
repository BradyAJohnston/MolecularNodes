import numpy as np
import molecularnodes as mn
from molecularnodes.utils import frame_mapper


def test_correct_1d():
    assert np.allclose(
        mn.utils.correct_periodic_1d(np.array((0.9, 0.1)), np.array((0.1, 0.9)), 1.0),
        np.array((1.1, -0.1)),
    )


def test_frame_mapper_basic():
    assert frame_mapper(10) == 10
    assert frame_mapper(0) == 0
    assert frame_mapper(-2) == 0


def test_frame_mapper_with_offset():
    assert frame_mapper(10, offset=2) == 8
    assert frame_mapper(2, offset=2) == 0
    assert frame_mapper(0, offset=5) == 0


def test_frame_mapper_with_subframes():
    assert frame_mapper(10, subframes=1) == 5
    assert frame_mapper(9, subframes=2) == 3


def test_frame_mapper_with_mapping():
    mapping = np.array([0, 0, 0, 1, 2])
    assert frame_mapper(1, mapping=mapping) == 0
    assert frame_mapper(3, mapping=mapping) == 1
    assert frame_mapper(3, mapping=list(mapping)) == 1


def test_frame_mapper_with_mapping_and_subframes():
    mapping = np.array([0, 0, 0, 1, 2])
    assert frame_mapper(5, subframes=1, mapping=mapping) == 0
