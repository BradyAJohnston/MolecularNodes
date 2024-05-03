import molecularnodes as mn
import numpy as np


def test_random_rgb(snapshot):
    n = 100
    colors = np.array(list(map(
        lambda x: mn.color.random_rgb(x),
        range(n)
    )))
    assert colors.tolist() == snapshot


def test_colos_from_atomic_numbers(snapshot):
    length = len(mn.color.iupac_colors_rgb)
    colors = mn.color.colors_from_elements(np.array(list(range(length))))
    assert colors.tolist() == snapshot
