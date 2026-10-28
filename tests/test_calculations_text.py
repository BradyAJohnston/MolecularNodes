import bpy
import numpy as np
import pytest
import molecularnodes as mn
from .constants import data_dir

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
        bpy.context.scene.frame_set(frame)
        atoms = mol.universe.atoms
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
    bpy.ops.mn.calculations_from_text(uuid=mol.uuid, text=text.name)
    text.from_string(CODE.replace("positions[:, 2]", "positions[:, 0]"))
    bpy.ops.mn.calculations_from_text(uuid=mol.uuid, text=text.name)
    assert np.allclose(
        mol.named_attribute("ZPosition"), mol.universe.atoms.positions[:, 0]
    )


def test_operator_reports_missing_text(mol):
    with pytest.raises(RuntimeError, match="Choose a text block"):
        bpy.ops.mn.calculations_from_text(uuid=mol.uuid, text="missing.py")


@pytest.mark.parametrize("updating", [True, False])
def test_calculations_from_text_survive_reopening(mol, tmp_path, updating):
    mol.calculations.from_text(new_text("calcs.py", CODE), updating=updating)
    uuid = mol.uuid
    filepath = tmp_path / "calcs.blend"
    bpy.ops.wm.save_as_mainfile(filepath=str(filepath))

    session = mn.session.get_session()
    session.clear()
    bpy.ops.wm.open_mainfile(filepath=str(filepath))

    restored = session.get(uuid)
    assert restored is not None
    bpy.context.scene.frame_set(2)
    assert set(restored.calculations) >= {"ZPosition", "IsHeavy"}
    assert restored.calculations.updating("ZPosition") is updating
    assert np.allclose(
        restored.named_attribute("ZPosition"),
        restored.universe.atoms.positions[:, 2],
    )
