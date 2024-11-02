import numpy as np


def centre(position: np.ndarray, weight: np.ndarray | None = None):
    "Calculate the weighted centroid of the vectors"
    if weight is None:
        return np.mean(position, axis=0)
    return np.sum(position * weight.reshape((-1, 1)), axis=0) / np.sum(weight)
