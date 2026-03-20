import bpy
import MDAnalysis as mda
import molecularnodes.ui.ops as ops


def test_import_smiles_calls_session_add_trajectory(monkeypatch):
    calls = {}

    class DummySession:
        def add_trajectory(self, universe, name, style):
            calls["universe"] = universe
            calls["name"] = name
            calls["style"] = style

            class DummyTraj:
                object = None

            return DummyTraj()

    def fake_get_session():
        return DummySession()

    class DummyUniverse:
        pass

    def fake_from_smiles(smiles):
        calls["smiles"] = smiles
        return DummyUniverse()

    monkeypatch.setattr(ops, "get_session", fake_get_session)
    monkeypatch.setattr(
        mda.Universe, "from_smiles", staticmethod(fake_from_smiles), raising=False
    )

    result = bpy.ops.mn.import_smiles(
        smiles="C1CCCCC1",
        name="Cyclohexane",
        style="spheres",
        setup_nodes=False,
    )

    assert result == {"FINISHED"}
    assert calls["smiles"] == "C1CCCCC1"
    assert calls["name"] == "Cyclohexane"
    assert calls["style"] is None
    assert isinstance(calls["universe"], DummyUniverse)
