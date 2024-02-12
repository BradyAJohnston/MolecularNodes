import bpy
import numpy as np
import pytest
import itertools
import tempfile
import molecularnodes as mn
from .constants import (
    data_dir,
    codes,
    attributes
)
from .utils import get_verts, sample_attribute_to_string, sample_attribute

mn.unregister()
mn.register()

styles = ['preset_1', 'cartoon', 'ribbon',
          'spheres', 'surface', 'ball_and_stick']


def useful_function(snapshot, style, code, assembly, cache_dir=None):
    obj = mn.io.fetch(
        code,
        style=style,
        build_assembly=assembly,
        cache_dir=cache_dir
    ).object
    node = mn.blender.nodes.get_style_node(obj)

    if 'EEVEE' in node.inputs.keys():
        node.inputs['EEVEE'].default_value = True

    mn.blender.nodes.realize_instances(obj)
    dont_realise = style == 'cartoon' and code == '1BNA'
    for att in attributes:
        snapshot.assert_match(
            sample_attribute_to_string(obj, att, evaluate=dont_realise),
            f"{att}.txt"
        )


with tempfile.TemporaryDirectory() as temp_dir:
    @pytest.mark.parametrize("assembly, code, style", itertools.product([False], codes, styles))
    def test_style_1(snapshot, assembly, code, style):
        useful_function(snapshot, style, code, assembly, cache_dir=temp_dir)

    # have to test a subset of styles with the biological assembly.
    # testing some of the heavier styles run out of memory and fail on github actions
    @pytest.mark.parametrize("assembly, code, style", itertools.product([True], codes, ['cartoon', 'surface', 'ribbon']))
    def test_style_2(snapshot, assembly, code, style):
        useful_function(snapshot, style, code, assembly, cache_dir=temp_dir)

    @pytest.mark.parametrize("code, format", itertools.product(codes, ['mmtf', 'cif', 'pdb']))
    def test_download_format(code, format):
        mol = mn.io.fetch(code, format=format, style=None).object
        scene = bpy.context.scene
        scene.MN_pdb_code = code
        scene.MN_import_node_setup = False
        scene.MN_import_format_download = format
        names = [o.name for o in bpy.data.objects]
        bpy.ops.mn.import_wwpdb()

        for o in bpy.data.objects:
            if o.name not in names:
                mol2 = o

        def verts(object):
            return mn.blender.obj.get_attribute(object, 'position')
        assert np.isclose(verts(mol), verts(mol2)).all()


def test_local_pdb(snapshot):
    molecules = [
        mn.io.load(data_dir / f'1l58.{ext}', style='spheres')
        for ext in ('cif', 'pdb')
    ]
    molecules.append(mn.io.fetch('1l58', format='mmtf'))
    for att in ['position']:
        for mol in molecules:
            snapshot.assert_match(
                sample_attribute_to_string(
                    mol, att, evaluate=False, precision=3),
                f'1L58_{att}'
            )


def test_rcsb_nmr(snapshot):
    mol = mn.io.fetch('2M6Q', style='cartoon')
    assert len(mol.frames.objects) == 10
    assert mol.object.modifiers['MolecularNodes'].node_group.nodes['MN_animate_value'].inputs['To Max'].default_value == 9

    snapshot.assert_match(
        sample_attribute(mol, 'position', evaluate=True, as_string=True),
        'position.txt'
    )

    pos_1 = mol.get_attribute('position', evaluate=True)
    bpy.context.scene.frame_set(100)
    pos_2 = mol.get_attribute('position', evaluate=True)
    assert (pos_1 != pos_2).all()


def test_load_small_mol(snapshot):
    mol = mn.io.load(data_dir / "ASN.cif")
    for att in ['position', 'bond_type']:
        snapshot.assert_match(sample_attribute(
            mol, att, as_string=True), f'{att}.txt')


def test_rcsb_cache(snapshot):
    from pathlib import Path
    import tempfile
    import os
    # we want to make sure cached files are freshly downloaded, but
    # we don't want to delete our entire real cache
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        test_cache = Path(temp_dir)

        # Run the test
        obj_1 = mn.io.fetch('6BQN', style='cartoon', cache_dir=test_cache)
        file = os.path.join(test_cache, '6BQN.mmtf')
        assert os.path.exists(file)

        obj_2 = mn.io.fetch('6BQN', style='cartoon', cache_dir=test_cache)
        assert (sample_attribute(
            obj_1, 'position') == sample_attribute(obj_2, 'position')).all()
