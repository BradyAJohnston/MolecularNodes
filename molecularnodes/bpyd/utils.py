import numpy as np


def centre(position: np.ndarray, weight: np.ndarray | None = None):
    "Calculate the weighted centroid of the vectors"
    if weight is None:
        return np.mean(position, axis=0)
    return np.sum(position * weight.reshape((-1, 1)), axis=0) / np.sum(weight)


def lerp(a: np.ndarray, b: np.ndarray, t: float = 0.5) -> np.ndarray:
    """
    Linearly interpolate between two values.

    Parameters
    ----------
    a : array_like
        The starting value.
    b : array_like
        The ending value.
    t : float, optional
        The interpolation parameter. Default is 0.5.

    Returns
    -------
    array_like
        The interpolated value(s).

    Notes
    -----
    This function performs linear interpolation between `a` and `b` using the
    interpolation parameter `t` such that the result lies between `a` and `b`.

    Examples
    --------
    >>> lerp(1, 2, 0.5)
    1.5

    >>> lerp(3, 7, 0.2)
    3.8

    >>> lerp([1, 2, 3], [4, 5, 6], 0.5)
    array([2.5, 3.5, 4.5])

    """
    return np.add(a, np.multiply(np.subtract(b, a), t))
