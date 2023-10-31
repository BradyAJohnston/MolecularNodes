import random
import colorsys
import numpy as np

def random_rgb(seed=None):
    """Random Pastel RGB values
    """
    random.seed(seed)
    r, g, b = colorsys.hls_to_rgb(random.random(), 0.6, 0.6)
    return np.array((r, g, b, 1))


def color_from_atomic_number(atomic_number: int):
    r, g, b = list(iupac_colors_rgb.values())[int(atomic_number - 1)]
    return np.array((r, g, b, 1))


def colors_from_elements(atomic_numbers):
    colors = np.array(list(map(color_from_atomic_number, atomic_numbers)))
    return colors


def equidistant_colors(some_list):
    u = np.unique(some_list)
    num_colors = len(u)

    hues = [i / (num_colors + 1) for i in range(num_colors)]

    # Convert HSL to RGB
    colors = [colorsys.hls_to_rgb(hue, 0.6, 0.6) for hue in hues]

    # Convert RGB to 8-bit integer values
    colors = [(int(r * 255), int(g * 255), int(b * 255), 1) for (r, g, b) in colors]

    return dict(zip(u, colors))


def color_chains(atomic_numbers, chain_ids):
    mask = atomic_numbers == 6
    colors = colors_from_elements(atomic_numbers)
    chain_color_dict = equidistant_colors(chain_ids)
    chain_colors = np.array(list(map(
        lambda x: chain_color_dict.get(x),
        chain_ids
    )))

    colors[mask] = chain_colors[mask]

    return colors / 255


iupac_colors_rgb = {
    "H": (255, 255, 255),   # Hydrogen
    "He": (217, 255, 255),  # Helium
    "Li": (204, 128, 255),  # Lithium
    "Be": (194, 255, 0),    # Beryllium
    "B": (255, 181, 181),   # Boron
    "C": (144, 144, 144),   # Carbon
    "N": (48, 80, 248),     # Nitrogen
    "O": (255, 13, 13),     # Oxygen
    "F": (144, 224, 80),    # Fluorine
    "Ne": (179, 227, 245),  # Neon
    "Na": (171, 92, 242),   # Sodium
    "Mg": (138, 255, 0),    # Magnesium
    "Al": (191, 166, 166),  # Aluminum
    "Si": (240, 200, 160),  # Silicon
    "P": (255, 128, 0),     # Phosphorus
    "S": (255, 255, 48),    # Sulfur
    "Cl": (31, 240, 31),    # Chlorine
    "K": (143, 64, 212),    # Potassium
    "Ar": (128, 209, 227),  # Argon
    "Ca": (61, 255, 0),     # Calcium
    "Sc": (230, 230, 230),  # Scandium
    "Ti": (191, 194, 199),  # Titanium
    "V": (166, 166, 171),   # Vanadium
    "Cr": (138, 153, 199),  # Chromium
    "Mn": (156, 122, 199),  # Manganese
    "Fe": (224, 102, 51),   # Iron
    "Ni": (199, 138, 138),  # Nickel
    "Co": (255, 217, 143),  # Cobalt
    "Cu": (200, 128, 51),   # Copper
    "Zn": (125, 128, 176),  # Zinc
    "Ga": (194, 143, 143),  # Gallium
    "Ge": (102, 143, 143),  # Germanium
    "As": (189, 128, 227),  # Arsenic
    "Se": (255, 161, 0),    # Selenium
    "Br": (166, 41, 41),    # Bromine
    "Kr": (92, 184, 209),   # Krypton
    "Rb": (112, 46, 176),   # Rubidium
    "Sr": (0, 255, 0),      # Strontium
    "Y": (148, 255, 255),   # Yttrium
    "Zr": (148, 224, 224),   # Zirconium
    "Nb": (115, 194, 201),   # Niobium
    "Mo": (84, 181, 181),   # Molybdenum
    "Tc": (59, 158, 158),   # Technetium
    "Ru": (36, 125, 125),   # Ruthenium
    "Rh": (10, 125, 140),   # Rhodium
    "Pd": (0, 105, 133),    # Palladium
    "Ag": (192, 192, 192),  # Silver
    "Cd": (255, 217, 143),  # Cadmium
    "In": (166, 117, 115),  # Indium
    "Sn": (102, 128, 128),  # Tin
    "Sb": (158, 99, 181),   # Antimony
    "Te": (212, 122, 0),    # Tellurium
    "I": (148, 0, 148),     # Iodine
    "Xe": (66, 158, 176),   # Xenon
    "Cs": (87, 23, 143),    # Cesium
    "Ba": (0, 201, 0),      # Barium
    "La": (112, 212, 255),  # Lanthanum
    "Ce": (255, 255, 199),  # Cerium
    "Pr": (217, 255, 199),  # Praseodymium
    "Nd": (199, 255, 199),  # Neodymium
    "Pm": (163, 255, 199),  # Promethium
    "Sm": (143, 255, 199),  # Samarium
    "Eu": (97, 255, 199),   # Europium
    "Gd": (69, 255, 199),   # Gadolinium
    "Tb": (48, 255, 199),   # Terbium
    "Dy": (31, 255, 199),   # Dysprosium
    "Ho": (0, 255, 156),    # Holmium
    "Er": (0, 230, 117),    # Erbium
    "Tm": (0, 212, 82),     # Thulium
    "Yb": (0, 191, 56),     # Ytterbium
    "Lu": (0, 171, 36),     # Lutetium
    "Hf": (77, 194, 255),   # Hafnium
    "Ta": (77, 166, 255),   # Tantalum
    "W": (33, 148, 214),    # Tungsten
    "Re": (38, 125, 171),    # Rhenium
    "Os": (38, 102, 150),    # Osmium
    "Ir": (23, 84, 135),    # Iridium
    "Pt": (208, 208, 224),  # Platinum
    "Au": (255, 209, 35),   # Gold
    "Hg": (184, 184, 208),  # Mercury
    "Tl": (166, 84, 77),    # Thallium
    "Pb": (87, 89, 97),     # Lead
    "Bi": (158, 79, 181),    # Bismuth
    "Th": (255, 161, 0),   # Thorium
    "Pa": (255, 161, 0),   # Protactinium
    "U": (255, 161, 0),    # Uranium
    "Np": (255, 161, 0),   # Neptunium
    "Pu": (255, 161, 0),   # Plutonium
    "Am": (255, 161, 0),   # Americium
    "Cm": (255, 161, 0),   # Curium
    "Bk": (255, 161, 0),   # Berkelium
    "Cf": (255, 161, 0),   # Californium
    "Es": (255, 161, 0),   # Einsteinium
    "Fm": (255, 161, 0),   # Fermium
    "Md": (255, 161, 0),   # Mendelevium
    "No": (255, 161, 0),   # Nobelium
    "Lr": (255, 161, 0),   # Lawrencium
    "Rf": (204, 0, 89),    # Rutherfordium
    "Db": (209, 0, 79),    # Dubnium
    "Sg": (217, 0, 69),    # Seaborgium
    "Bh": (224, 0, 56),    # Bohrium
    "Hs": (230, 0, 46),    # Hassium
    "Mt": (235, 0, 38),    # Meitnerium
    "Ds": (240, 0, 33),    # Darmstadtium
    "Rg": (241, 0, 30),    # Roentgenium
    "Cn": (242, 0, 26),    # Copernicium
    "Nh": (242, 0, 26),    # Nihonium
    "Fl": (242, 0, 26),    # Flerovium
    "Mc": (242, 0, 26),    # Moscovium
    "Lv": (242, 0, 26),    # Livermorium
    "Ts": (242, 0, 26),    # Tennessine
    "Og": (242, 0, 26)     # Oganesson
}
