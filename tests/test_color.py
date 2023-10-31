import molecularnodes as mn
import numpy as np

def test_random_rgb(snapshot):
    n = 100
    threshold = n * 4
    colors = np.array(list(map(
        lambda x: mn.color.random_rgb(x), 
        range(n)
    )))
    snapshot.assert_match(
        np.array2string(colors, precision=3, threshold=threshold), 
        "color_values.txt"
    )

def test_colos_from_atomic_numbers(snapshot):
    length = len(mn.color.iupac_colors_rgb)
    threshold = length * 4
    colors = mn.color.colors_from_elements(np.array(list(range(length))))
    snapshot.assert_match(
        np.array2string(colors, precision=3, threshold=threshold), 
        "atomic_number_colors.txt"
    )