import bpy
import pytest
from databpy import ObjectTracker
import molecularnodes as mn


def test_trajectory_from_smiles():
    traj = mn.Trajectory.from_smiles("CCO", name="Ethanol", style=None)

    assert traj.name == "Ethanol"
    assert traj.object.name == "Ethanol"
    assert traj.universe.atoms.n_atoms == len(traj.object.data.vertices)
    assert bpy.context.scene.MNSession.get(traj.uuid) is traj


def test_import_smiles_basic():
    session = bpy.context.scene.MNSession

    with ObjectTracker() as tracker:
        result = bpy.ops.mn.import_smiles(
            smiles="C",
            name="Methane",
            setup_nodes=False,
        )
        traj = session.match(tracker.latest())

    assert result == {"FINISHED"}
    assert traj.name == "Methane"
    assert traj.object.name == "Methane"
    assert traj.universe.atoms.n_atoms == len(traj.object.data.vertices)
    assert bpy.context.view_layer.objects.active == traj.object


def test_import_smiles_empty():
    with pytest.raises(RuntimeError, match="SMILES string is empty"):
        bpy.ops.mn.import_smiles(smiles="   ", name="Empty")


def test_import_smiles_name_fallback():
    session = bpy.context.scene.MNSession

    with ObjectTracker() as tracker:
        result = bpy.ops.mn.import_smiles(
            smiles="C",
            name="   ",
            setup_nodes=False,
        )
        traj = session.match(tracker.latest())

    assert result == {"FINISHED"}
    assert traj.name == "SMILES"
    assert traj.object.name == "SMILES"


def test_import_smiles_parse_failure():
    with pytest.raises(RuntimeError, match="Failed to parse SMILES:"):
        bpy.ops.mn.import_smiles(smiles="not-a-smiles", name="Broken")
