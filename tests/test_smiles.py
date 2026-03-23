import bpy
import MDAnalysis as mda
import pytest
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


def test_import_smiles_missing_support(monkeypatch):
    def fake_hasattr(_obj, _name):
        return False

    monkeypatch.setattr(ops, "hasattr", fake_hasattr)
    with pytest.raises(
        RuntimeError,
        match="MDAnalysis does not support SMILES in this version",
    ):
        bpy.ops.mn.import_smiles(smiles="C1CCCCC1", name="Cyclohexane")


def test_import_smiles_empty_input(monkeypatch):
    class DummyUniverse:
        pass

    def fake_from_smiles(_smiles):
        return DummyUniverse()

    monkeypatch.setattr(
        mda.Universe, "from_smiles", staticmethod(fake_from_smiles), raising=False
    )
    with pytest.raises(RuntimeError, match="SMILES string is empty"):
        bpy.ops.mn.import_smiles(smiles="   ", name="Empty")


def test_import_smiles_name_fallback(monkeypatch):
    calls = {}

    class DummySession:
        def add_trajectory(self, _universe, name, style):
            calls["name"] = name
            calls["style"] = style

            class DummyTraj:
                object = None

            return DummyTraj()

    def fake_get_session():
        return DummySession()

    class DummyUniverse:
        pass

    def fake_from_smiles(_smiles):
        return DummyUniverse()

    monkeypatch.setattr(ops, "get_session", fake_get_session)
    monkeypatch.setattr(
        mda.Universe, "from_smiles", staticmethod(fake_from_smiles), raising=False
    )

    result = bpy.ops.mn.import_smiles(smiles="C1CCCCC1", name="   ")
    assert result == {"FINISHED"}
    assert calls["name"] == "SMILES"


def test_import_smiles_parse_failure(monkeypatch):
    def mock_fail(_smiles):
        raise ValueError("invalid smiles")

    monkeypatch.setattr(mda.Universe, "from_smiles", mock_fail, raising=False)
    with pytest.raises(RuntimeError, match="Failed to parse SMILES: invalid smiles"):
        bpy.ops.mn.import_smiles(smiles="bad", name="Broken")
