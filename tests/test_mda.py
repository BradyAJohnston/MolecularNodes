import bpy
import MDAnalysis as mda
import pytest
from molecularnodes import MDAVis
from .constants import data_dir


@pytest.fixture
def mdavis():
    return MDAVis()


@pytest.fixture()
def universe():
    top = data_dir / "md_ppr/box.gro"
    traj = data_dir / "md_ppr/first_5_frames.xtc"
    u = mda.Universe(top, traj)
    return u


def test_add_universe(mdavis, universe):
    # add whole universe
    mdavis.add_universe(universe)
    assert "u0" in bpy.data.objects
    # add atomgroup
    ag = universe.select_atoms("name CA")
    mdavis.add_universe(ag, name="agUniverse")
    assert "agUniverse" in bpy.data.objects


def test_universe_attributes(mdavis, universe):
    vu = mdavis.add_universe(universe)
    # test blender object
    assert bpy.data.objects["u0"] == vu.object
    # test name changes
    assert vu.name == vu.object.name
    vu.name = "customName"
    assert bpy.data.objects["customName"] == vu.object
    # test key
    assert vu.key == "u0"


def test_subscriptable_universes(mdavis, universe):
    vu = mdavis.add_universe(universe, name="u")
    assert mdavis.universes["u"] == vu


def test_iterable_universes(mdavis, universe):
    vu1 = mdavis.add_universe(universe, name="u1")
    vu2 = mdavis.add_universe(universe, name="u2")
    ulist = []
    for u in mdavis.universes:
        ulist.append(u)
    assert vu1 in ulist
    assert vu2 in ulist


def test_operators():
    topo = str(data_dir / "md_ppr/box.gro")
    traj = str(data_dir / "md_ppr/first_5_frames.xtc")
    # add universe operator
    bpy.ops.mda.add_universe(topology=topo, trajectory=traj, name="u")
    assert "u" in bpy.data.objects
    # delete universe operator
    bpy.ops.mda.delete_universe()
    assert "u" not in bpy.data.objects
