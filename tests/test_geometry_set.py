import databpy as db
import numpy as np
import molecularnodes as mn
from .utils import GeometrySet


def test_get_set(snapshot):
    code = "1BNA"
    mol = mn.Molecule.fetch(code).add_style("ribbon")
    geom = GeometrySet(mol.object)
    assert snapshot == geom


def test_spheres_style(snapshot):
    code = "1BNA"
    mol = mn.Molecule.fetch(code).add_style("spheres")
    geom = GeometrySet(mol.object)
    assert snapshot == geom


def test_specific_attribute():
    code = "1BNA"
    mol = mn.Molecule.fetch(code).add_style("ribbon")
    geom = GeometrySet(mol.object)

    mesh = geom.geom.mesh
    b_factor_attr = db.Attribute(mesh.attributes["b_factor"])
    b_factor = b_factor_attr.as_array()

    assert b_factor.shape[0] == len(mesh.vertices)
    assert b_factor.dtype == np.float32
    assert b_factor.min() > 0
    assert b_factor.max() < 200

    instances = geom.instances
    if ".reference_index" in instances.attributes:
        ref_attr = db.Attribute(instances.attributes[".reference_index"])
        ref_index = ref_attr.as_array()
        n_unique = len(np.unique(ref_index))
        assert n_unique == 1


def test_attribute_comparison():
    code = "1BNA"

    mol_ribbon = mn.Molecule.fetch(code).add_style("ribbon")
    geom_ribbon = GeometrySet(mol_ribbon.object)

    mol_spheres = mn.Molecule.fetch(code).add_style("spheres")
    geom_spheres = GeometrySet(mol_spheres.object)

    pc = geom_spheres.geom.pointcloud
    pc_count = len(db.Attribute(pc.attributes["position"]))

    instances = geom_ribbon.instances
    inst_count = len(db.Attribute(instances.attributes["position"]))

    assert pc_count > inst_count

    pc_atomic = db.Attribute(pc.attributes["atomic_number"]).as_array()
    inst_atomic = db.Attribute(instances.attributes["atomic_number"]).as_array()

    assert len(np.unique(pc_atomic)) > 1
    assert len(np.unique(inst_atomic)) == 1
    assert inst_atomic[0] == 6
