"""Color manipulation and palette generation utilities for molecular visualization.

This module provides tools for color space conversions, palette generation,
and molecular structure coloring using various schemes including IUPAC colors,
equidistant hue distributions, and LAB color space operations.
"""

import colorsys
import math
import random
from typing import Dict, List, Optional, Sequence, Tuple, Union
import numpy as np
import numpy.typing as npt

# Color constants
RGB_MAX_VALUE = 255
ALPHA_OPAQUE = 1.0
PASTEL_LIGHTNESS = 0.6
PASTEL_SATURATION = 0.6
CARBON_ATOMIC_NUMBER = 6

# Palette generation constants
MAX_DISTINCT_COLORS = 20  # Reasonable limit before looping with offset
HUE_OFFSET_FACTOR = 0.618034  # Golden ratio for pleasing offset distribution


def clamp(value: float, min_value: float, max_value: float) -> float:
    """Clamp a value between minimum and maximum bounds.

    Parameters
    ----------
    value : float
        The value to clamp.
    min_value : float
        The minimum allowed value.
    max_value : float
        The maximum allowed value.

    Returns
    -------
    float
        The clamped value.
    """
    return max(min_value, min(value, max_value))


class Lab:
    """LAB color space representation and conversion utilities.

    The LAB color space (also known as CIELAB) is a perceptually uniform color space
    that allows for better color manipulation and distance calculations than RGB.

    Attributes
    ----------
    KN_CONSTANT : float
        Lightness scaling constant.
    ILLUMINANT_X : float
        X component of D65 illuminant.
    ILLUMINANT_Y : float
        Y component of D65 illuminant.
    ILLUMINANT_Z : float
        Z component of D65 illuminant.
    THRESHOLD_* : float
        Various thresholds for color space conversions.
    """

    # LAB color space constants
    KN_CONSTANT = 18
    ILLUMINANT_X = 0.950470  # D65 illuminant X
    ILLUMINANT_Y = 1.0  # D65 illuminant Y
    ILLUMINANT_Z = 1.088830  # D65 illuminant Z
    THRESHOLD_T0 = 0.137931034
    THRESHOLD_T1 = 0.206896552
    THRESHOLD_T2 = 0.12841855
    THRESHOLD_T3 = 0.008856452

    def __init__(
        self, lightness: float = 0.1, green_red: float = 0.0, blue_yellow: float = 0.0
    ) -> None:
        """Initialize a LAB color.

        Parameters
        ----------
        lightness : float, default=0.1
            L* component (lightness) from 0 to 100.
        green_red : float, default=0.0
            a* component (green-red axis).
        blue_yellow : float, default=0.0
            b* component (blue-yellow axis).
        """
        self.lab = [lightness, green_red, blue_yellow]
        self.l = lightness
        self.a = green_red
        self.b = blue_yellow

    @staticmethod
    def zero() -> "Lab":
        """Create a LAB color representing black.

        Returns
        -------
        Lab
            A LAB color with all components set to zero.
        """
        return Lab(0, 0, 0)

    @staticmethod
    def distance(color_a: "Lab", color_b: "Lab") -> float:
        """Calculate Euclidean distance between two LAB colors.

        Parameters
        ----------
        color_a : Lab
            First color.
        color_b : Lab
            Second color.

        Returns
        -------
        float
            Euclidean distance between the colors.
        """
        lightness_diff = color_b.l - color_a.l
        green_red_diff = color_b.a - color_a.a
        blue_yellow_diff = color_b.b - color_a.b
        return math.sqrt(lightness_diff**2 + green_red_diff**2 + blue_yellow_diff**2)

    @staticmethod
    def darken(output_color: "Lab", input_color: "Lab", amount: float) -> "Lab":
        """Darken a LAB color by reducing its lightness.

        Parameters
        ----------
        output_color : Lab
            The LAB color to store the result in.
        input_color : Lab
            The input color to darken.
        amount : float
            Amount to darken (positive values darken).

        Returns
        -------
        Lab
            The darkened color (same as output_color).
        """
        output_color.l = input_color.l - Lab.KN_CONSTANT * amount
        output_color.a = input_color.a
        output_color.b = input_color.b
        return output_color

    @staticmethod
    def lighten(output_color: "Lab", input_color: "Lab", amount: float) -> "Lab":
        """Lighten a LAB color by increasing its lightness.

        Parameters
        ----------
        output_color : Lab
            The LAB color to store the result in.
        input_color : Lab
            The input color to lighten.
        amount : float
            Amount to lighten (positive values lighten).

        Returns
        -------
        Lab
            The lightened color (same as output_color).
        """
        return Lab.darken(output_color, input_color, -amount)

    @staticmethod
    def darken_color(color: npt.NDArray[np.float32], amount: float) -> List[float]:
        """Darken an RGB color using LAB color space.

        Parameters
        ----------
        color : npt.NDArray[np.float32]
            RGBA color array with values in [0, 1] range.
        amount : float
            Amount to darken (positive values darken).

        Returns
        -------
        List[float]
            Darkened RGBA color as list.
        """
        temp_lab_color = Lab.from_color(color)
        return Lab.to_color(Lab.darken(temp_lab_color, temp_lab_color, amount))

    @staticmethod
    def lighten_color(color: npt.NDArray[np.float32], amount: float) -> List[float]:
        """Lighten an RGB color using LAB color space.

        Parameters
        ----------
        color : npt.NDArray[np.float32]
            RGBA color array with values in [0, 1] range.
        amount : float
            Amount to lighten (positive values lighten).

        Returns
        -------
        List[float]
            Lightened RGBA color as list.
        """
        return Lab.darken_color(color, -amount)

    @staticmethod
    def from_color(color: npt.NDArray[np.float32]) -> "Lab":
        """Convert RGBA color to LAB color space.

        Parameters
        ----------
        color : npt.NDArray[np.float32]
            RGBA color array with values in [0, 1] range.

        Returns
        -------
        Lab
            Color in LAB space.
        """
        red, green, blue, alpha = color * RGB_MAX_VALUE
        xyz_x, xyz_y, xyz_z = Lab.rgb_to_xyz(red, green, blue)
        lightness = 116 * xyz_y - 16
        return Lab(
            lightness if lightness >= 0 else 0,
            500 * (xyz_x - xyz_y),
            200 * (xyz_y - xyz_z),
        )

    @staticmethod
    def to_color(lab_color: "Lab") -> List[float]:
        """Convert LAB color to RGBA color space.

        Parameters
        ----------
        lab_color : Lab
            Color in LAB space.

        Returns
        -------
        List[float]
            RGBA color as list with values in [0, 1] range.
        """
        # Convert LAB to XYZ
        xyz_y_norm = (lab_color.l + 16) / 116
        xyz_x_norm = (
            xyz_y_norm if math.isnan(lab_color.a) else xyz_y_norm + lab_color.a / 500
        )
        xyz_z_norm = (
            xyz_y_norm if math.isnan(lab_color.b) else xyz_y_norm - lab_color.b / 200
        )

        # Apply illuminant scaling
        xyz_y = Lab.ILLUMINANT_Y * Lab.lab_xyz(xyz_y_norm)
        xyz_x = Lab.ILLUMINANT_X * Lab.lab_xyz(xyz_x_norm)
        xyz_z = Lab.ILLUMINANT_Z * Lab.lab_xyz(xyz_z_norm)

        # XYZ to RGB transformation matrix (sRGB)
        red = Lab.xyz_rgb(3.2404542 * xyz_x - 1.5371385 * xyz_y - 0.4985314 * xyz_z)
        green = Lab.xyz_rgb(-0.9692660 * xyz_x + 1.8760108 * xyz_y + 0.0415560 * xyz_z)
        blue = Lab.xyz_rgb(0.0556434 * xyz_x - 0.2040259 * xyz_y + 1.0572252 * xyz_z)

        return [
            round(clamp(red, 0, RGB_MAX_VALUE)) / RGB_MAX_VALUE,
            round(clamp(green, 0, RGB_MAX_VALUE)) / RGB_MAX_VALUE,
            round(clamp(blue, 0, RGB_MAX_VALUE)) / RGB_MAX_VALUE,
            ALPHA_OPAQUE,
        ]

    @staticmethod
    def xyz_rgb(xyz_component: float) -> float:
        """Convert XYZ component to RGB component using sRGB gamma correction.

        Parameters
        ----------
        xyz_component : float
            XYZ color component value.

        Returns
        -------
        float
            RGB component in [0, 255] range.
        """
        return RGB_MAX_VALUE * (
            12.92 * xyz_component
            if xyz_component <= 0.00304
            else 1.055 * math.pow(xyz_component, 1 / 2.4) - 0.055
        )

    @staticmethod
    def lab_xyz(normalized_component: float) -> float:
        """Convert normalized LAB component to XYZ.

        Parameters
        ----------
        normalized_component : float
            Normalized LAB component.

        Returns
        -------
        float
            XYZ component value.
        """
        return (
            normalized_component**3
            if normalized_component > Lab.THRESHOLD_T1
            else Lab.THRESHOLD_T2 * (normalized_component - Lab.THRESHOLD_T0)
        )

    @staticmethod
    def rgb_xyz(rgb_component: float) -> float:
        """Convert RGB component to XYZ using inverse sRGB gamma correction.

        Parameters
        ----------
        rgb_component : float
            RGB component in [0, 255] range.

        Returns
        -------
        float
            XYZ component value.
        """
        normalized = rgb_component / RGB_MAX_VALUE
        return (
            normalized / 12.92
            if normalized <= 0.04045
            else math.pow((normalized + 0.055) / 1.055, 2.4)
        )

    @staticmethod
    def xyz_lab(xyz_component: float) -> float:
        """Convert XYZ component to normalized LAB component.

        Parameters
        ----------
        xyz_component : float
            XYZ component value.

        Returns
        -------
        float
            Normalized LAB component.
        """
        return (
            math.pow(xyz_component, 1 / 3)
            if xyz_component > Lab.THRESHOLD_T3
            else xyz_component / Lab.THRESHOLD_T2 + Lab.THRESHOLD_T0
        )

    @staticmethod
    def rgb_to_xyz(red: float, green: float, blue: float) -> List[float]:
        """Convert RGB color to XYZ color space.

        Parameters
        ----------
        red : float
            Red component in [0, 255] range.
        green : float
            Green component in [0, 255] range.
        blue : float
            Blue component in [0, 255] range.

        Returns
        -------
        List[float]
            XYZ components as [X, Y, Z].
        """
        # Convert RGB to linear RGB
        linear_red = Lab.rgb_xyz(red)
        linear_green = Lab.rgb_xyz(green)
        linear_blue = Lab.rgb_xyz(blue)

        # Apply sRGB to XYZ transformation matrix
        xyz_x = Lab.xyz_lab(
            (
                0.4124564 * linear_red
                + 0.3575761 * linear_green
                + 0.1804375 * linear_blue
            )
            / Lab.ILLUMINANT_X
        )
        xyz_y = Lab.xyz_lab(
            (
                0.2126729 * linear_red
                + 0.7151522 * linear_green
                + 0.0721750 * linear_blue
            )
            / Lab.ILLUMINANT_Y
        )
        xyz_z = Lab.xyz_lab(
            (
                0.0193339 * linear_red
                + 0.1191920 * linear_green
                + 0.9503041 * linear_blue
            )
            / Lab.ILLUMINANT_Z
        )

        return [xyz_x, xyz_y, xyz_z]


def random_rgb(seed: Optional[int] = None) -> npt.NDArray[np.float32]:
    """Generate a random pastel RGB color.

    Parameters
    ----------
    seed : Optional[int], default=None
        Random seed for reproducible color generation.
        If None, uses system time.

    Returns
    -------
    npt.NDArray[np.float32]
        RGBA color array with values in [0, 1] range.

    Notes
    -----
    Uses HLS color space with fixed lightness and saturation to generate
    pleasing pastel colors with random hues.
    """
    if seed is not None:
        random.seed(seed)

    hue = random.random()  # Random hue from 0 to 1
    red, green, blue = colorsys.hls_to_rgb(hue, PASTEL_LIGHTNESS, PASTEL_SATURATION)
    return np.array([red, green, blue, ALPHA_OPAQUE], dtype=np.float32)


def plddt(confidence_scores: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
    """Generate colors based on AlphaFold2 pLDDT confidence scores.

    Parameters
    ----------
    confidence_scores : npt.NDArray[np.float32]
        Array of pLDDT confidence scores (0-100 scale).

    Returns
    -------
    npt.NDArray[np.float32]
        Array of RGBA colors corresponding to confidence levels.
        Shape: (N, 4) where N is the number of input scores.

    Notes
    -----
    Color scheme follows AlphaFold2 convention:
    - Blue (>90): Very high confidence
    - Light blue (70-90): High confidence
    - Yellow (50-70): Low confidence
    - Red (<50): Very low confidence
    """
    num_scores = len(confidence_scores)
    colors = np.zeros((num_scores, 4), dtype=np.float32)

    # Define confidence color thresholds and corresponding RGBA values
    very_high_confidence = np.array([0.000000, 0.086496, 0.672395, 1.000000])  # Blue
    high_confidence = np.array([0.130157, 0.597176, 0.896205, 1.000000])  # Light blue
    low_confidence = np.array([1.000169, 0.708345, 0.006512, 1.000000])  # Yellow
    very_low_confidence = np.array([1.000169, 0.205070, 0.059507, 1.000000])  # Red

    for i, score in enumerate(confidence_scores):
        if score > 90:
            colors[i, :] = very_high_confidence
        elif score > 70:
            colors[i, :] = high_confidence
        elif score > 50:
            colors[i, :] = low_confidence
        else:
            colors[i, :] = very_low_confidence

    return colors


def color_from_atomic_number(atomic_number: int) -> npt.NDArray[np.int32]:
    """Get IUPAC standard color for an element by atomic number.

    Parameters
    ----------
    atomic_number : int
        Atomic number of the element (1-based indexing).

    Returns
    -------
    npt.NDArray[np.int32]
        RGBA color array with values in [0, 255] range.

    Raises
    ------
    IndexError
        If atomic_number is out of range for available elements.
    """
    element_colors = list(iupac_colors_rgb.values())
    red, green, blue = element_colors[atomic_number - 1]  # Convert to 0-based indexing
    return np.array([red, green, blue, RGB_MAX_VALUE], dtype=np.int32)


def color_from_element(element_symbol: str) -> npt.NDArray[np.int32]:
    """Get IUPAC standard color for an element by symbol.

    Parameters
    ----------
    element_symbol : str
        Chemical element symbol (e.g., 'C', 'N', 'O').

    Returns
    -------
    npt.NDArray[np.int32]
        RGBA color array with values in [0, 255] range.

    Raises
    ------
    KeyError
        If element_symbol is not found in IUPAC color dictionary.
    """
    red, green, blue = iupac_colors_rgb[element_symbol]
    return np.array([red, green, blue, RGB_MAX_VALUE], dtype=np.int32)


def colors_from_elements(atomic_numbers: Sequence[int]) -> npt.NDArray[np.int32]:
    """Get IUPAC colors for multiple elements by atomic numbers.

    Parameters
    ----------
    atomic_numbers : Sequence[int]
        Sequence of atomic numbers.

    Returns
    -------
    npt.NDArray[np.int32]
        Array of RGBA colors with shape (N, 4) and values in [0, 255] range.
    """
    color_array = np.array([color_from_atomic_number(num) for num in atomic_numbers])
    return color_array


def generate_equidistant_color_palette(
    unique_identifiers: Sequence[Union[str, int]],
    max_colors: int = MAX_DISTINCT_COLORS,
    lightness: float = PASTEL_LIGHTNESS,
    saturation: float = PASTEL_SATURATION,
    return_format: str = "rgba_int",
) -> Dict[
    Union[str, int], Union[Tuple[int, int, int, int], Tuple[float, float, float, float]]
]:
    """Generate a color palette with equidistant hues for unique identifiers.

    This function creates visually distinct colors by distributing hues
    evenly around the color wheel. For large numbers of unique values,
    it intelligently loops back with golden ratio offsets to maintain
    visual distinction.

    Parameters
    ----------
    unique_identifiers : Sequence[Union[str, int]]
        Sequence of unique identifiers to assign colors to.
        Duplicates will be automatically removed.
    max_colors : int, default=20
        Maximum number of distinct colors before looping with offset.
    lightness : float, default=0.6
        HLS lightness value (0.0 to 1.0).
    saturation : float, default=0.6
        HLS saturation value (0.0 to 1.0).
    return_format : str, default="rgba_int"
        Return format: "rgba_int" for [0,255] or "rgba_float" for [0,1].

    Returns
    -------
    Dict[Union[str, int], Union[Tuple[int, int, int, int], Tuple[float, float, float, float]]]
        Dictionary mapping each unique identifier to an RGBA color tuple.

    Notes
    -----
    Uses HLS color space with configurable lightness and saturation.
    For more than max_colors unique values, cycles through hues with
    golden ratio offsets to maintain visual distinction.

    Examples
    --------
    >>> palette = generate_equidistant_color_palette(['A', 'B', 'C'])
    >>> len(palette)
    3
    >>> # Test with many values - should loop with offset
    >>> many_ids = [f"chain_{i}" for i in range(25)]
    >>> big_palette = generate_equidistant_color_palette(many_ids)
    >>> len(big_palette)
    25
    """
    unique_ids = np.unique(unique_identifiers)
    num_unique = len(unique_ids)

    if num_unique == 0:
        return {}

    colors = []

    for i in range(num_unique):
        if num_unique <= max_colors:
            # Simple equidistant distribution
            hue = (i + 1) / (max_colors + 1)
        else:
            # For large numbers, use cycling with golden ratio offset
            cycle = i // max_colors
            position_in_cycle = i % max_colors
            base_hue = (position_in_cycle + 1) / (max_colors + 1)
            # Add golden ratio offset for each cycle to maintain distinction
            hue_offset = (cycle * HUE_OFFSET_FACTOR) % 1.0
            hue = (base_hue + hue_offset) % 1.0

        # Convert HLS to RGB
        red, green, blue = colorsys.hls_to_rgb(hue, lightness, saturation)

        if return_format == "rgba_float":
            colors.append((red, green, blue, ALPHA_OPAQUE))
        else:  # rgba_int
            colors.append(
                (
                    int(red * RGB_MAX_VALUE),
                    int(green * RGB_MAX_VALUE),
                    int(blue * RGB_MAX_VALUE),
                    RGB_MAX_VALUE,
                )
            )

    return dict(zip(unique_ids, colors))


def create_color_array(
    identifiers: Sequence[Union[str, int]],
    return_format: str = "rgba_int",
    **palette_kwargs,
) -> npt.NDArray[Union[np.int32, np.float32]]:
    """Create a color array for a sequence of identifiers.

    This is a consolidated function that handles the common pattern of:
    1. Generate a color palette for unique identifiers
    2. Map each identifier to its corresponding color
    3. Return as a numpy array

    Parameters
    ----------
    identifiers : Sequence[Union[str, int]]
        Sequence of identifiers to color.
    return_format : str, default="rgba_int"
        Return format: "rgba_int" for [0,255] or "rgba_float" for [0,1].
    **palette_kwargs
        Additional arguments passed to generate_equidistant_color_palette.

    Returns
    -------
    npt.NDArray[Union[np.int32, np.float32]]
        Array of RGBA colors with shape (N, 4).

    Examples
    --------
    >>> colors = create_color_array(['A', 'B', 'A', 'C'])
    >>> colors.shape
    (4, 4)
    >>> np.array_equal(colors[0], colors[2])  # Same ID, same color
    True
    """
    color_palette = generate_equidistant_color_palette(
        identifiers, return_format=return_format, **palette_kwargs
    )

    colors = np.array(
        [color_palette[identifier] for identifier in identifiers],
        dtype=np.int32 if return_format == "rgba_int" else np.float32,
    )

    return colors


def color_chains_equidistant(
    chain_identifiers: Sequence[Union[str, int]],
) -> npt.NDArray[np.int32]:
    """Assign equidistant colors to protein chain identifiers.

    This is a convenience wrapper around create_color_array.

    Parameters
    ----------
    chain_identifiers : Sequence[Union[str, int]]
        Sequence of chain identifiers (e.g., ['A', 'B', 'A', 'C']).

    Returns
    -------
    npt.NDArray[np.int32]
        Array of RGBA colors with shape (N, 4) where N is the length
        of chain_identifiers. Values are in [0, 255] range.
    """
    return create_color_array(chain_identifiers, return_format="rgba_int")


def color_chains_with_atoms(
    atomic_numbers: Sequence[int], chain_identifiers: Sequence[Union[str, int]]
) -> npt.NDArray[np.float32]:
    """Color atoms using element colors, but override carbon with chain colors.

    This creates a hybrid coloring scheme where non-carbon atoms use their
    standard IUPAC element colors, while carbon atoms (which are abundant
    in protein structures) are colored by chain to improve visualization.

    Parameters
    ----------
    atomic_numbers : Sequence[int]
        Atomic numbers for each atom.
    chain_identifiers : Sequence[Union[str, int]]
        Chain identifiers for each atom.

    Returns
    -------
    npt.NDArray[np.float32]
        Array of RGBA colors with shape (N, 4) and values in [0, 1] range.

    Notes
    -----
    Carbon atoms (atomic number 6) are colored by chain using equidistant
    hues, while all other atoms use their standard IUPAC element colors.
    This approach provides clear chain distinction while maintaining
    chemical intuition for heteroatoms.

    Examples
    --------
    >>> atomic_nums = [6, 7, 6, 8]  # C, N, C, O
    >>> chains = ['A', 'A', 'B', 'B']
    >>> colors = color_chains_with_atoms(atomic_nums, chains)
    >>> colors.shape
    (4, 4)
    """
    # Get element colors for all atoms
    element_colors = colors_from_elements(atomic_numbers)

    # Generate chain colors using the consolidated function
    chain_colors = create_color_array(chain_identifiers, return_format="rgba_int")

    # Create mask for carbon atoms (atomic number 6)
    carbon_atom_mask = np.array(atomic_numbers) == CARBON_ATOMIC_NUMBER

    # Override carbon colors with chain colors
    element_colors[carbon_atom_mask] = chain_colors[carbon_atom_mask]

    # Convert to [0, 1] range for consistency with other color functions
    return element_colors.astype(np.float32) / RGB_MAX_VALUE


# ============================================================================
# UNIFIED COLOR PALETTE API
# ============================================================================


def create_palette(
    identifiers: Sequence[Union[str, int]], palette_type: str = "equidistant", **kwargs
) -> Union[Dict, npt.NDArray]:
    """Unified color palette generation function.

    This is the main function for generating color palettes in MolecularNodes.
    It provides a simple, consistent interface for all color generation needs.

    Parameters
    ----------
    identifiers : Sequence[Union[str, int]]
        Sequence of identifiers to assign colors to.
    palette_type : str, default="equidistant"
        Type of palette: "equidistant", "random", or "element".
    **kwargs
        Additional arguments passed to specific palette functions.
        Common options:
        - return_format: "rgba_int" (default) or "rgba_float"
        - as_array: bool, if True returns numpy array instead of dict
        - max_colors: int, maximum distinct colors before cycling
        - lightness: float, HLS lightness (0-1)
        - saturation: float, HLS saturation (0-1)

    Returns
    -------
    Union[Dict, npt.NDArray]
        Color palette as dictionary or array based on parameters.

    Examples
    --------
    >>> # Simple equidistant palette
    >>> colors = create_palette(['A', 'B', 'C'])
    >>>
    >>> # Large dataset with cycling
    >>> many_chains = [f"chain_{i}" for i in range(50)]
    >>> colors = create_palette(many_chains, as_array=True)
    >>>
    >>> # Custom lightness and saturation
    >>> colors = create_palette(['X', 'Y'], lightness=0.8, saturation=0.4)
    """
    as_array = kwargs.pop("as_array", False)

    if palette_type == "equidistant":
        if as_array:
            return create_color_array(identifiers, **kwargs)
        else:
            return generate_equidistant_color_palette(identifiers, **kwargs)
    elif palette_type == "random":
        # For random, generate one color per unique identifier
        unique_ids = np.unique(identifiers)
        seed = kwargs.get("seed")
        palette = {}
        for i, uid in enumerate(unique_ids):
            color_seed = None if seed is None else seed + i
            rgb_color = random_rgb(color_seed)
            if kwargs.get("return_format") == "rgba_int":
                palette[uid] = tuple((rgb_color * RGB_MAX_VALUE).astype(int))
            else:
                palette[uid] = tuple(rgb_color)

        if as_array:
            colors = [palette[uid] for uid in identifiers]
            dtype = (
                np.int32 if kwargs.get("return_format") == "rgba_int" else np.float32
            )
            return np.array(colors, dtype=dtype)
        return palette
    elif palette_type == "element":
        if as_array:
            return colors_from_elements(identifiers)
        else:
            unique_atomic_nums = np.unique(identifiers)
            return {
                num: tuple(color_from_atomic_number(num)) for num in unique_atomic_nums
            }
    else:
        raise ValueError(f"Unknown palette_type: {palette_type}")


# ============================================================================
# BACKWARD COMPATIBILITY ALIASES
# ============================================================================

# These aliases preserve the original API while using the new implementations
equidistant_colors = generate_equidistant_color_palette
color_chains = color_chains_with_atoms


# Additional convenience aliases
def get_chain_colors(chain_ids: Sequence[Union[str, int]]) -> npt.NDArray[np.int32]:
    """Simple alias for chain coloring - most common use case."""
    return create_color_array(chain_ids)


def get_element_colors(atomic_numbers: Sequence[int]) -> npt.NDArray[np.int32]:
    """Simple alias for element coloring."""
    return colors_from_elements(atomic_numbers)


# IUPAC standard colors for chemical elements (RGB values in 0-255 range)
iupac_colors_rgb = {
    "H": (255, 255, 255),  # Hydrogen
    "He": (217, 255, 255),  # Helium
    "Li": (204, 128, 255),  # Lithium
    "Be": (194, 255, 0),  # Beryllium
    "B": (255, 181, 181),  # Boron
    "C": (144, 144, 144),  # Carbon
    "N": (48, 80, 248),  # Nitrogen
    "O": (255, 13, 13),  # Oxygen
    "F": (144, 224, 80),  # Fluorine
    "Ne": (179, 227, 245),  # Neon
    "Na": (171, 92, 242),  # Sodium
    "Mg": (138, 255, 0),  # Magnesium
    "Al": (191, 166, 166),  # Aluminum
    "Si": (240, 200, 160),  # Silicon
    "P": (255, 128, 0),  # Phosphorus
    "S": (255, 255, 48),  # Sulfur
    "Cl": (31, 240, 31),  # Chlorine
    "K": (143, 64, 212),  # Potassium
    "Ar": (128, 209, 227),  # Argon
    "Ca": (61, 255, 0),  # Calcium
    "Sc": (230, 230, 230),  # Scandium
    "Ti": (191, 194, 199),  # Titanium
    "V": (166, 166, 171),  # Vanadium
    "Cr": (138, 153, 199),  # Chromium
    "Mn": (156, 122, 199),  # Manganese
    "Fe": (224, 102, 51),  # Iron
    "Ni": (199, 138, 138),  # Nickel
    "Co": (255, 217, 143),  # Cobalt
    "Cu": (200, 128, 51),  # Copper
    "Zn": (125, 128, 176),  # Zinc
    "Ga": (194, 143, 143),  # Gallium
    "Ge": (102, 143, 143),  # Germanium
    "As": (189, 128, 227),  # Arsenic
    "Se": (255, 161, 0),  # Selenium
    "Br": (166, 41, 41),  # Bromine
    "Kr": (92, 184, 209),  # Krypton
    "Rb": (112, 46, 176),  # Rubidium
    "Sr": (0, 255, 0),  # Strontium
    "Y": (148, 255, 255),  # Yttrium
    "Zr": (148, 224, 224),  # Zirconium
    "Nb": (115, 194, 201),  # Niobium
    "Mo": (84, 181, 181),  # Molybdenum
    "Tc": (59, 158, 158),  # Technetium
    "Ru": (36, 125, 125),  # Ruthenium
    "Rh": (10, 125, 140),  # Rhodium
    "Pd": (0, 105, 133),  # Palladium
    "Ag": (192, 192, 192),  # Silver
    "Cd": (255, 217, 143),  # Cadmium
    "In": (166, 117, 115),  # Indium
    "Sn": (102, 128, 128),  # Tin
    "Sb": (158, 99, 181),  # Antimony
    "Te": (212, 122, 0),  # Tellurium
    "I": (148, 0, 148),  # Iodine
    "Xe": (66, 158, 176),  # Xenon
    "Cs": (87, 23, 143),  # Cesium
    "Ba": (0, 201, 0),  # Barium
    "La": (112, 212, 255),  # Lanthanum
    "Ce": (255, 255, 199),  # Cerium
    "Pr": (217, 255, 199),  # Praseodymium
    "Nd": (199, 255, 199),  # Neodymium
    "Pm": (163, 255, 199),  # Promethium
    "Sm": (143, 255, 199),  # Samarium
    "Eu": (97, 255, 199),  # Europium
    "Gd": (69, 255, 199),  # Gadolinium
    "Tb": (48, 255, 199),  # Terbium
    "Dy": (31, 255, 199),  # Dysprosium
    "Ho": (0, 255, 156),  # Holmium
    "Er": (0, 230, 117),  # Erbium
    "Tm": (0, 212, 82),  # Thulium
    "Yb": (0, 191, 56),  # Ytterbium
    "Lu": (0, 171, 36),  # Lutetium
    "Hf": (77, 194, 255),  # Hafnium
    "Ta": (77, 166, 255),  # Tantalum
    "W": (33, 148, 214),  # Tungsten
    "Re": (38, 125, 171),  # Rhenium
    "Os": (38, 102, 150),  # Osmium
    "Ir": (23, 84, 135),  # Iridium
    "Pt": (208, 208, 224),  # Platinum
    "Au": (255, 209, 35),  # Gold
    "Hg": (184, 184, 208),  # Mercury
    "Tl": (166, 84, 77),  # Thallium
    "Pb": (87, 89, 97),  # Lead
    "Bi": (158, 79, 181),  # Bismuth
    "Th": (255, 161, 0),  # Thorium
    "Pa": (255, 161, 0),  # Protactinium
    "U": (255, 161, 0),  # Uranium
    "Np": (255, 161, 0),  # Neptunium
    "Pu": (255, 161, 0),  # Plutonium
    "Am": (255, 161, 0),  # Americium
    "Cm": (255, 161, 0),  # Curium
    "Bk": (255, 161, 0),  # Berkelium
    "Cf": (255, 161, 0),  # Californium
    "Es": (255, 161, 0),  # Einsteinium
    "Fm": (255, 161, 0),  # Fermium
    "Md": (255, 161, 0),  # Mendelevium
    "No": (255, 161, 0),  # Nobelium
    "Lr": (255, 161, 0),  # Lawrencium
    "Rf": (204, 0, 89),  # Rutherfordium
    "Db": (209, 0, 79),  # Dubnium
    "Sg": (217, 0, 69),  # Seaborgium
    "Bh": (224, 0, 56),  # Bohrium
    "Hs": (230, 0, 46),  # Hassium
    "Mt": (235, 0, 38),  # Meitnerium
    "Ds": (240, 0, 33),  # Darmstadtium
    "Rg": (241, 0, 30),  # Roentgenium
    "Cn": (242, 0, 26),  # Copernicium
    "Nh": (242, 0, 26),  # Nihonium
    "Fl": (242, 0, 26),  # Flerovium
    "Mc": (242, 0, 26),  # Moscovium
    "Lv": (242, 0, 26),  # Livermorium
    "Ts": (242, 0, 26),  # Tennessine
    "Og": (242, 0, 26),  # Oganesson
}
