import gzip
import bpy
import numpy as np
from nodebpy.nodes.geometry import RealizeInstances
import molecularnodes as mn
from molecularnodes.entities.reload import reload_entity
from molecularnodes.nodes.geometry import EnsembleInstance
from .constants import data_dir

# `placements.json` and `packed.gro.gz` were written by `bentopy-pack` and
# `bentopy-render` (bentopy 1.0.1) from `input.bent` in the same directory
bentopy_dir = data_dir / "bentopy"


def read_gro_positions(path) -> np.ndarray:
    with gzip.open(path, "rt") as f:
        lines = f.read().splitlines()
    n_atoms = int(lines[1])
    return np.array(
        [
            [float(a[20:28]), float(a[28:36]), float(a[36:44])]
            for a in lines[2:][:n_atoms]
        ]
    )


def test_load_bentopy():
    ens = mn.entities.ensemble.Bentopy.load(
        bentopy_dir / "placements.json", node_setup=False
    )
    assert ens._entity_type == mn.entities.base.EntityType.ENSEMBLE_BENTOPY
    assert ens.props.entity_type == ens._entity_type.value
    assert ens.props.chain_ids == ["lyso", "dna"]
    assert len(ens.object.data.vertices) == 12
    obj_names = sorted(obj.name for obj in ens.instance_collection.objects)
    assert obj_names == ["0000_lyso", "0001_dna"]

    with ens.tree.reset() as (atoms, join):
        instance = EnsembleInstance(instances=ens.instance_collection, as_points=False)
        atoms >> instance >> RealizeInstances() >> join

    # bentopy works in nm, which is one Blender unit at MN's world scale
    expected = read_gro_positions(bentopy_dir / "packed.gro.gz")
    pos_eval = ens.named_attribute("position", evaluate=True)
    assert pos_eval.shape == expected.shape
    # .gro positions are written to 3 decimal places
    np.testing.assert_allclose(pos_eval, expected, atol=1e-3)


def test_load_bentopy_tutorial():
    # Tutorial 1 of "Bentopy: from simple packing to building cellular models"
    # (cgmartini.nl): 650 Martini lysozymes in a 40 nm box. `lysozyme.pdb` is from
    # the tutorial's `tutorial_files.tar.gz`, `simple_packing.bent` is copied from
    # the tutorial and `placements.json` was written from it by `bentopy-pack`
    # (bentopy 1.0.1). The full render is ~200k atoms, so `packed_subset.gro.gz`
    # keeps only the `bentopy-render` lines for the instances listed below.
    tutorial_dir = bentopy_dir / "tutorial"
    ens = mn.entities.ensemble.Bentopy.load(
        tutorial_dir / "placements.json", node_setup=False
    )
    assert ens.props.chain_ids == ["LYZ"]
    assert len(ens.object.data.vertices) == 650

    with ens.tree.reset() as (atoms, join):
        instance = EnsembleInstance(instances=ens.instance_collection, as_points=False)
        atoms >> instance >> RealizeInstances() >> join

    n_atoms = 304
    pos_eval = ens.named_attribute("position", evaluate=True)
    assert pos_eval.shape == (650 * n_atoms, 3)
    picked = [0, 1, 324, 649]
    pos_picked = pos_eval.reshape(650, n_atoms, 3)[picked].reshape(-1, 3)
    expected = read_gro_positions(tutorial_dir / "packed_subset.gro.gz")
    np.testing.assert_allclose(pos_picked, expected, atol=1e-3)


def test_op_import_and_reload_bentopy():
    bpy.ops.mn.import_ensemble(
        ensemble_type="bentopy",
        filepath=str(bentopy_dir / "placements.json"),
        node_setup=True,
    )
    obj = bpy.data.objects["placements"]
    assert obj.mn.entity_type == "ensemble-bentopy"

    entity = reload_entity(obj)
    assert isinstance(entity, mn.entities.ensemble.Bentopy)
    assert entity.object == obj
    assert len(entity.transformations) == 12


def test_load_bentopy_twice_keeps_instances_separate():
    path = bentopy_dir / "placements.json"
    first = mn.entities.ensemble.Bentopy.load(path, name="packed")
    second = mn.entities.ensemble.Bentopy.load(path, name="packed")

    # the second object is renamed by Blender, and its instances must follow it
    assert second.object.name != first.object.name
    assert first.instance_collection != second.instance_collection
    for ens in (first, second):
        assert len(ens.instance_collection.objects) == 2
        trees = {obj.modifiers[0].node_group for obj in ens.instance_collection.objects}
        assert len(trees) == 1
        assert len(next(iter(trees)).interface.items_tree) == 2

    first_trees = {o.modifiers[0].node_group for o in first.instance_collection.objects}
    second_trees = {
        o.modifiers[0].node_group for o in second.instance_collection.objects
    }
    assert first_trees.isdisjoint(second_trees)
