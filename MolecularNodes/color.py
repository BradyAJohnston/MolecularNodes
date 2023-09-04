import random
import colorsys
import numpy as np

def random_rgb():
    """Random Pastel RGB values
    """
    r, g, b = colorsys.hls_to_rgb(random.random(), 0.6, 0.6)
    return np.array((r, g, b, 1))