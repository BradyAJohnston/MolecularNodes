import bpy
import MolecularNodes as mn
import pytest
from .utils import get_verts, realize_intances

codes = ['6n2y', '4ozs', '8H1B', '1BNA']

@pytest.mark.parametrize("preset", [1, 2, 3, 4])
@pytest.mark.parametrize("code", codes)
def test_preset_style(snapshot, preset, code):
    obj = mn.load.molecule_rcsb(code, starting_style=f"presets")
    node = obj.modifiers['MolecularNodes'].node_group.nodes['MN_style_presets']
    node.inputs['Preset'].default_value = preset
    node.inputs['Atom: Eevee / Cycles'].default_value = True
    realize_intances(obj)
    verts = get_verts(obj, apply_modifiers=True, float_decimals=3)
    snapshot.assert_match(verts, f"preset_{preset}_{code}.txt")