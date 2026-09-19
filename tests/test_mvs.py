import json
import warnings
import numpy as np
import pytest
from molviewspec import create_builder
import molecularnodes as mn
from molecularnodes.entities import mvs
from tests.utils import GeometrySet
from .constants import data_dir


def _write_mvsj(tmp_path, builder, name="scene.mvsj"):
    path = tmp_path / name
    path.write_text(builder.get_state().dumps())
    return path


@pytest.fixture
def basic_builder():
    """A scene with two representations of a local structure file."""
    builder = create_builder()
    structure = builder.download(url="1BNA.bcif").parse(format="bcif").model_structure()
    structure.component(selector="nucleic").representation(type="cartoon").color(
        color="#00ff00"
    )
    structure.component(selector="water").representation(type="spacefill")
    return builder


def test_mvs_basic_import(tmp_path, basic_builder, snapshot):
    # the .mvsj must sit next to the structure it references relatively
    path = _write_mvsj(data_dir, basic_builder)
    try:
        molecules = mvs.load(path)
    finally:
        path.unlink()

    assert len(molecules) == 1
    mol = molecules[0]
    assert mol.name == "1BNA"

    # two style branches were created
    styles = [
        node
        for node in mol.node_group.nodes
        if getattr(node, "node_tree", None) and "Style" in node.node_tree.name
    ]
    assert sorted(node.node_tree.name for node in styles) == [
        "Style Cartoon",
        "Style Spheres",
    ]

    # the nucleic component's selection and baked color were stored
    attributes = mol.list_attributes(drop_hidden=False)
    selection_names = [name for name in attributes if name.startswith("mvs_selection")]
    color_names = [name for name in attributes if name.startswith("mvs_color")]
    assert len(selection_names) == 2
    assert len(color_names) == 1

    nucleic = mol.named_attribute("is_nucleic").astype(bool)
    assert np.array_equal(mol.named_attribute(selection_names[0]).astype(bool), nucleic)
    baked = mol.named_attribute(color_names[0])
    assert np.allclose(baked[nucleic, :3], (0.0, 1.0, 0.0), atol=1e-3)

    assert snapshot == GeometrySet(mol.object)


def test_mvs_expressions(tmp_path):
    builder = create_builder()
    structure = builder.download(url="1cd3.cif").parse(format="mmcif").model_structure()
    structure.component(
        selector={"auth_asym_id": "B", "beg_auth_seq_id": 1, "end_auth_seq_id": 20}
    ).representation(type="ball_and_stick")

    path = _write_mvsj(data_dir, builder)
    try:
        molecules = mvs.load(path)
    finally:
        path.unlink()

    mol = molecules[0]
    selection_names = [
        name
        for name in mol.list_attributes(drop_hidden=False)
        if name.startswith("mvs_selection")
    ]
    assert len(selection_names) == 1
    mask = mol.named_attribute(selection_names[0]).astype(bool)

    atoms = mol.universe.atoms
    expected = (atoms.chainIDs == "B") & (atoms.resids >= 1) & (atoms.resids <= 20)
    assert np.array_equal(mask, expected)
    assert mask.any()


def test_mvs_label_asym_selection(tmp_path):
    # label_asym_id selection resolved from the secondary parse of the cif
    builder = create_builder()
    structure = builder.download(url="1BNA.bcif").parse(format="bcif").model_structure()
    structure.component(selector={"label_asym_id": "A"}).representation(type="cartoon")

    path = _write_mvsj(data_dir, builder)
    try:
        molecules = mvs.load(path)
    finally:
        path.unlink()

    mol = molecules[0]
    selection_names = [
        name
        for name in mol.list_attributes(drop_hidden=False)
        if name.startswith("mvs_selection")
    ]
    mask = mol.named_attribute(selection_names[0]).astype(bool)
    assert mask.any()
    assert not mask.all()


def test_mvs_opacity_and_camera(tmp_path):
    import bpy

    builder = create_builder()
    builder.canvas(background_color="#101020")
    structure = builder.download(url="1BNA.bcif").parse(format="bcif").model_structure()
    (
        structure.component(selector="all")
        .representation(type="surface")
        .opacity(opacity=0.4)
    )
    builder.camera(target=(0, 0, 0), position=(0, 0, 100), up=(0, 1, 0))

    # the builder does not expose `near` yet, but the schema defines it
    state = json.loads(builder.get_state().dumps())
    for node in state["root"]["children"]:
        if node["kind"] == "camera":
            node["params"]["near"] = 30
    path = data_dir / "scene.mvsj"
    path.write_text(json.dumps(state))
    try:
        molecules = mvs.load(path)
    finally:
        path.unlink()

    mol = molecules[0]
    color_names = [
        name
        for name in mol.list_attributes(drop_hidden=False)
        if name.startswith("mvs_color")
    ]
    baked = mol.named_attribute(color_names[0])
    assert np.allclose(baked[:, 3], 0.4, atol=1e-3)

    # camera placed at the scaled position, looking down -Z toward the origin,
    # with the near clipping plane from the camera node's `near`
    camera = bpy.context.scene.camera
    assert np.allclose(np.array(camera.location), (0, 0, 10), atol=1e-4)
    assert camera.data.clip_start == pytest.approx(3.0)


def test_mvs_focus_clips_to_component(tmp_path):
    import bpy

    builder = create_builder()
    structure = builder.download(url="1BNA.bcif").parse(format="bcif").model_structure()
    component = structure.component(selector={"auth_seq_id": 8})
    component.representation(type="ball_and_stick")
    component.focus()
    structure.component(selector="all").representation(type="cartoon")

    path = _write_mvsj(data_dir, builder)
    try:
        (mol,) = mvs.load(path)
    finally:
        path.unlink()

    # the camera looks down -Z (the MVS focus default) at the residue, with the
    # near plane pushed up to the front of its bounding sphere so anything in
    # front of it is clipped away
    camera = bpy.context.scene.camera
    basis_forward = np.array(
        camera.matrix_world.to_3x3() @ __import__("mathutils").Vector((0, 0, -1))
    )
    assert np.allclose(basis_forward, (0, 0, -1), atol=1e-5)

    selections = [
        name
        for name in mol.list_attributes(drop_hidden=False)
        if name.startswith("mvs_selection")
    ]
    mask = mol.named_attribute(selections[0]).astype(bool)
    points = mol.named_attribute("position")[mask]
    center = (points.min(axis=0) + points.max(axis=0)) / 2
    depth = float((center - np.array(camera.location)) @ basis_forward)
    assert 0 < camera.data.clip_start < depth


def test_mvs_unsupported_warns(tmp_path):
    builder = create_builder()
    structure = builder.download(url="1BNA.bcif").parse(format="bcif").model_structure()
    component = structure.component(selector="all")
    component.representation(type="cartoon")
    component.label(text="a label")

    path = _write_mvsj(data_dir, builder)
    try:
        with pytest.warns(mvs.MVSImportWarning, match="label"):
            molecules = mvs.load(path)
    finally:
        path.unlink()
    assert len(molecules) == 1


def test_mvs_missing_source_warns(tmp_path):
    builder = create_builder()
    builder.download(url="does_not_exist.bcif").parse(format="bcif").model_structure()
    path = _write_mvsj(tmp_path, builder)
    with pytest.warns(mvs.MVSImportWarning, match="not found"):
        molecules = mvs.load(path)
    assert molecules == []


def test_mvs_entity_reload_roundtrip(tmp_path):
    """Imported molecules are ordinary entities: styled, evaluated, framed."""
    builder = create_builder()
    structure = builder.download(url="1BNA.bcif").parse(format="bcif").model_structure()
    structure.component(selector="nucleic").representation(type="cartoon")

    path = _write_mvsj(data_dir, builder)
    try:
        (mol,) = mvs.load(path)
    finally:
        path.unlink()

    geo = GeometrySet(mol.object)
    assert geo.mesh is not None
    assert mn.entities.base.EntityType.MOLECULE == mol._entity_type


# ---- official examples --------------------------------------------------- #
# Vendored from https://molstar.org/mol-view-spec/ (the landing-page examples
# at landing/public/examples/*/state.mvsj and test-data/colab_examples in the
# molstar/mol-view-spec repository). They download their structures from
# wwPDB, so a first run needs network access; downloads are cached.

EXAMPLES = {
    "colab_components": 1,
    "colab_geometrical": 0,
    "colab_labels": 1,
    "colab_minimal": 1,
    "colab_superimpose": 2,
    "colab_volumetric": 1,
    "landing_annotations": 1,
    "landing_basic": 1,
    "landing_components": 1,
    "landing_label": 1,
    "landing_primitives": 0,
    "landing_superposition": 2,
    "landing_symmetry": 1,
    "landing_volumes": 1,
}


@pytest.mark.parametrize("name,n_structures", EXAMPLES.items())
def test_mvs_official_examples(name, n_structures):
    """Every official example imports without error, with the expected number
    of structures; unsupported parts only ever warn."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", mvs.MVSImportWarning)
        molecules = mvs.load(data_dir / "mvs" / f"{name}.mvsj")
    assert len(molecules) == n_structures
    for mol in molecules:
        assert mol.object is not None
        assert len(mol.universe.atoms) > 0


@pytest.mark.parametrize(
    "name,pattern",
    [
        ("landing_label", "'label' nodes"),
        ("landing_superposition", "'transform' nodes"),
        ("landing_volumes", "'volume' nodes"),
        ("landing_annotations", "'component_from_uri' nodes"),
        ("landing_symmetry", "structure type 'symmetry'"),
    ],
)
def test_mvs_official_examples_warn(name, pattern):
    with pytest.warns(mvs.MVSImportWarning, match=pattern):
        mvs.load(data_dir / "mvs" / f"{name}.mvsj")


def test_mvs_components_example_selections():
    """The landing 'components' example end to end: five styled branches, with
    static and label-identifier selections resolved against the structure."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", mvs.MVSImportWarning)
        (mol,) = mvs.load(data_dir / "mvs" / "landing_components.mvsj")

    styles = [
        node.node_tree.name
        for node in mol.node_group.nodes
        if getattr(node, "node_tree", None) and "Style" in node.node_tree.name
    ]
    assert sorted(styles) == [
        "Style Ball and Stick",
        "Style Ball and Stick",
        "Style Ball and Stick",
        "Style Cartoon",
        "Style Cartoon",
    ]

    attributes = mol.list_attributes(drop_hidden=False)
    selections = [
        mol.named_attribute(name).astype(bool)
        for name in attributes
        if name.startswith("mvs_selection")
    ]
    colors = [name for name in attributes if name.startswith("mvs_color")]
    assert len(colors) == 5

    # the 'protein' and 'nucleic' components match the import-time attributes
    is_peptide = mol.named_attribute("is_peptide").astype(bool)
    is_nucleic = mol.named_attribute("is_nucleic").astype(bool)
    assert any(np.array_equal(mask, is_peptide) for mask in selections)
    assert any(np.array_equal(mask, is_nucleic) for mask in selections)

    # the single-residue {label_asym_id, label_seq_id} selections resolve to
    # exactly one residue each
    single_residue = [
        mask
        for mask in selections
        if 0 < mask.sum() < 50 and len(np.unique(mol.universe.atoms.resids[mask])) == 1
    ]
    assert len(single_residue) >= 2
