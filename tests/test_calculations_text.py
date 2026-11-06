import bpy
import numpy as np
import pytest
import molecularnodes as mn
from .constants import data_dir
from .utils import mn_ops, scene

CODE = """
class ZPosition(BaseCalculation):
    def calculate(self):
        return self.universe.atoms.positions[:, 2]


class IsHeavy(BaseCalculation):
    def calculate(self):
        return self.universe.atoms.masses > 2
"""


@pytest.fixture
def mol():
    return mn.Molecule.load(
        data_dir / "md_ppr/box.gro", data_dir / "md_ppr/first_5_frames.xtc"
    )


def new_text(name: str, code: str) -> bpy.types.Text:
    text = bpy.data.texts.new(name)
    text.from_string(code)
    return text


def test_from_text_adds_each_subclass(mol):
    text = new_text("calcs.py", CODE)
    assert mol.calculations.from_text(text) == ["ZPosition", "IsHeavy"]

    # stored straight away, before any frame change
    attributes = mol.object.data.attributes
    assert attributes["ZPosition"].data_type == "FLOAT"
    assert attributes["IsHeavy"].data_type == "BOOLEAN"
    assert mol.object.mn.calculations["ZPosition"].text == text

    for frame in (1, 3):
        scene().frame_set(frame)
        atoms = mol.atoms
        assert np.allclose(mol.named_attribute("ZPosition"), atoms.positions[:, 2])
        assert np.array_equal(mol.named_attribute("IsHeavy"), atoms.masses > 2)


def test_from_text_rejects_wrong_length(mol):
    text = new_text(
        "bad.py",
        CODE
        + """
class OneAtom(BaseCalculation):
    def calculate(self):
        return self.universe.atoms.positions[2, :]
""",
    )
    with pytest.raises(ValueError, match="'OneAtom' returned 3 values"):
        mol.calculations.from_text(text)
    # nothing from the text is added when one of its calculations fails
    assert "ZPosition" not in mol.calculations
    assert "OneAtom" not in mol.calculations


def test_from_text_without_subclasses(mol):
    with pytest.raises(ValueError, match="No BaseCalculation subclasses"):
        mol.calculations.from_text(new_text("empty.py", "x = 1"))


def test_reload_replaces_calculation(mol):
    text = new_text("calcs.py", CODE)
    mn_ops.calculations_from_text(uuid=mol.uuid, text=text.name)
    text.from_string(CODE.replace("positions[:, 2]", "positions[:, 0]"))
    mn_ops.calculations_from_text(uuid=mol.uuid, text=text.name)
    assert np.allclose(mol.named_attribute("ZPosition"), mol.atoms.positions[:, 0])


def test_reload_keeps_updating_setting(mol):
    text = new_text("calcs.py", CODE)
    mol.calculations.from_text(text)
    mol.calculations.set_updating("ZPosition", False)
    mn_ops.calculations_from_text(uuid=mol.uuid, text=text.name)
    assert mol.calculations.updating("ZPosition") is False
    assert mol.calculations.updating("IsHeavy") is True


def test_operator_reports_missing_text(mol):
    with pytest.raises(RuntimeError, match="Choose a text block"):
        mn_ops.calculations_from_text(uuid=mol.uuid, text="missing.py")


def save_and_reopen(mol, tmp_path, updating=True) -> mn.Molecule:
    mol.calculations.from_text(new_text("calcs.py", CODE), updating=updating)
    uuid = mol.uuid
    filepath = tmp_path / "calcs.blend"
    bpy.ops.wm.save_as_mainfile(filepath=str(filepath))

    session = mn.session.get_session()
    session.clear()
    bpy.ops.wm.open_mainfile(filepath=str(filepath))

    restored = session.get(uuid)
    assert isinstance(restored, mn.Molecule)
    scene().frame_set(2)
    return restored


@pytest.fixture
def auto_run():
    """Turn on Blender's Auto Run Python Scripts for one test."""
    preferences = bpy.context.preferences
    assert preferences is not None
    before = preferences.filepaths.use_scripts_auto_execute
    preferences.filepaths.use_scripts_auto_execute = True
    yield preferences
    preferences.filepaths.use_scripts_auto_execute = before


@pytest.mark.parametrize("updating", [True, False])
def test_reopening_waits_for_reload(mol, tmp_path, updating):
    restored = save_and_reopen(mol, tmp_path, updating)
    # without Auto Run Python Scripts the text isn't run when the file opens
    assert "ZPosition" not in restored.calculations
    assert restored.frame_attributes["ZPosition"] == "CALCULATION"

    mn_ops.calculations_from_text(uuid=restored.uuid, text="calcs.py")
    assert set(restored.calculations) >= {"ZPosition", "IsHeavy"}
    assert restored.calculations.updating("ZPosition") is updating
    # Reload stores the values for the current frame straight away
    assert np.allclose(
        restored.named_attribute("ZPosition"), restored.atoms.positions[:, 2]
    )


@pytest.mark.parametrize("updating", [True, False])
def test_reopening_with_auto_run(mol, tmp_path, auto_run, updating):
    restored = save_and_reopen(mol, tmp_path, updating)
    assert set(restored.calculations) >= {"ZPosition", "IsHeavy"}
    assert restored.calculations.updating("ZPosition") is updating
    assert np.allclose(
        restored.named_attribute("ZPosition"), restored.atoms.positions[:, 2]
    )


def test_reopening_in_excluded_path(mol, tmp_path, auto_run):
    excluded = auto_run.autoexec_paths.new()
    excluded.path = str(tmp_path)
    try:
        restored = save_and_reopen(mol, tmp_path)
    finally:
        auto_run.autoexec_paths.remove(excluded)
    assert "ZPosition" not in restored.calculations


def test_remove_unloaded_calculation(mol, tmp_path):
    restored = save_and_reopen(mol, tmp_path)
    mn_ops.frame_attribute_remove(uuid=restored.uuid, name="ZPosition")
    assert "ZPosition" not in restored.frame_attributes
