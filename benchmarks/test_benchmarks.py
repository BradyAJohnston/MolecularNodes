"""
Benchmarks for loading structures and evaluating the node groups that style them.

Run locally with `uv run pytest benchmarks --codspeed`, which reports wall-clock times.
In CI they run under CodSpeed's CPU simulation, which tracks results over time and
reports changes on pull requests.
"""

import itertools
from pathlib import Path
import bpy
import pytest
import molecularnodes as mn
from molecularnodes.converters import universe_from_atoms
from molecularnodes.entities.molecule.reader import read_structure
from molecularnodes.nodes import geometry as mg

DATA_DIR = Path(__file__).parent.parent / "tests" / "data"

# ~4,000 atoms, available as .pdb, .cif and .bcif
STRUCTURE = "8H1B"
STYLES = ["spheres", "cartoon", "ribbon", "surface", "sticks", "ball_and_stick"]


# loading ---------------------------------------------------------------------


@pytest.mark.benchmark(group="load")
@pytest.mark.parametrize("suffix", ["pdb", "cif", "bcif"])
def test_load(run, suffix):
    run(mn.Molecule.load, DATA_DIR / f"{STRUCTURE}.{suffix}")


@pytest.mark.benchmark(group="load")
def test_load_large(run):
    # ~10,000 atoms
    run(mn.Molecule.load, DATA_DIR / "1cd3.bcif")


@pytest.mark.benchmark(group="parse")
def test_read_structure(benchmark):
    # parsing only, to separate it from building the Blender object and node tree
    benchmark(read_structure, DATA_DIR / f"{STRUCTURE}.bcif")


@pytest.mark.benchmark(group="parse")
def test_universe_from_atoms(benchmark):
    # converting the parsed biotite AtomArray into an MDAnalysis Universe
    array = read_structure(DATA_DIR / f"{STRUCTURE}.bcif").array
    benchmark(universe_from_atoms, array)


@pytest.mark.benchmark(group="load")
@pytest.mark.parametrize("style", STYLES)
def test_load_with_style(run, style):
    run(mn.Molecule.load, DATA_DIR / f"{STRUCTURE}.bcif", style=style)


@pytest.mark.benchmark(group="load")
def test_load_with_style_color_assembly(run):
    # a style that also adds color and assembly nodes, each from the asset library
    def func():
        mol = mn.Molecule.load(DATA_DIR / f"{STRUCTURE}.bcif")
        mol.add_style("cartoon", color="common", assembly=True)

    run(func)


@pytest.mark.benchmark(group="load-trajectory")
def test_load_trajectory(run):
    run(
        mn.Molecule.load,
        DATA_DIR / "md_ppr/box.gro",
        DATA_DIR / "md_ppr/first_5_frames.xtc",
    )


@pytest.mark.benchmark(group="trajectory-frame")
def test_trajectory_frame_change(benchmark):
    mn.Molecule.load(
        DATA_DIR / "md_ppr/box.gro", DATA_DIR / "md_ppr/first_5_frames.xtc"
    )
    # step to a different frame on every call, so each one has to update positions
    frames = itertools.cycle(range(5))
    benchmark(lambda: bpy.context.scene.frame_set(next(frames)))


# node evaluation -------------------------------------------------------------


def _evaluate(obj: bpy.types.Object):
    "Re-run the object's modifiers, which evaluates its Geometry Nodes tree."
    obj.update_tag()
    return obj.evaluated_get(bpy.context.evaluated_depsgraph_get())


@pytest.mark.benchmark(group="evaluate-style")
@pytest.mark.parametrize("style", STYLES)
def test_evaluate_style(benchmark, style):
    mol = mn.Molecule.load(DATA_DIR / f"{STRUCTURE}.bcif").add_style(style)
    _evaluate(mol.object)
    benchmark(_evaluate, mol.object)


COLORS = {
    "common": "common",
    "rainbow": lambda: mg.ColorRainbow(),
    "secondary_structure": lambda: mg.ColorSecondaryStructure(),
    "res_name": lambda: mg.ColorResName(),
    "goodsell": lambda: mg.ColorGoodsell(),
}


@pytest.mark.benchmark(group="evaluate-color")
@pytest.mark.parametrize("color", COLORS.values(), ids=COLORS.keys())
def test_evaluate_color(benchmark, color):
    # spheres are the cheapest style, so the color node makes up more of the time
    mol = mn.Molecule.load(DATA_DIR / f"{STRUCTURE}.bcif").add_style(
        "spheres", color=color
    )
    _evaluate(mol.object)
    benchmark(_evaluate, mol.object)
