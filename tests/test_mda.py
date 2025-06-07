import bpy
import MDAnalysis as mda
import pytest
from molecularnodes import MDAVis
from .constants import data_dir


@pytest.fixture
def mdavis():
    mdavis = MDAVis()
    yield mdavis
    mdavis.universes.clear_universes()


@pytest.fixture()
def universe():
    topo = data_dir / "md_ppr/box.gro"
    traj = data_dir / "md_ppr/first_5_frames.xtc"
    u = mda.Universe(topo, traj)
    return u


def test_add_universe(mdavis, universe):
    # add whole universe
    mdavis.universes.add_universe(universe)
    assert "u0" in bpy.data.objects
    # add atomgroup
    ag = universe.select_atoms("name CA")
    mdavis.universes.add_universe(ag, name="agUniverse")
    assert "agUniverse" in bpy.data.objects


def test_delete_universe(mdavis, universe):
    # add universe
    u = mdavis.universes.add_universe(universe)
    assert "u0" in bpy.data.objects
    # delete universe
    mdavis.universes.delete_universe(u)
    assert "u0" not in bpy.data.objects


def test_universe_attributes(mdavis, universe):
    vu = mdavis.universes.add_universe(universe)
    # test blender object
    assert bpy.data.objects["u0"] == vu.object
    # test name changes
    assert vu.name == vu.object.name
    vu.name = "customName"
    assert bpy.data.objects["customName"] == vu.object
    # test internal key
    assert vu._key == "u0"


def test_subscriptable_universes(mdavis, universe):
    vu = mdavis.universes.add_universe(universe, name="u")
    assert mdavis.universes["u"] == vu


def test_iterable_universes(mdavis, universe):
    vu1 = mdavis.universes.add_universe(universe, name="u1")
    vu2 = mdavis.universes.add_universe(universe, name="u2")
    ulist = []
    for u in mdavis.universes:
        ulist.append(u)
    assert vu1 in ulist
    assert vu2 in ulist
    assert len(mdavis.universes) == 2


def test_operators():
    topo = str(data_dir / "md_ppr/box.gro")
    traj = str(data_dir / "md_ppr/first_5_frames.xtc")
    # add universe operator
    bpy.ops.mda.add_universe(topology=topo, trajectory=traj, name="u")
    assert "u" in bpy.data.objects
    # frame selected universe
    bpy.ops.mda.frame_selected_universe(index=0)
    assert bpy.context.active_object == bpy.data.objects["u"]
    # delete universe operator
    bpy.ops.mda.delete_universe()
    assert "u" not in bpy.data.objects
