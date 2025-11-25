import databpy as db
import molecularnodes as mn
from .utils import get_geometry_set


def test_get_set():
    code = "1BNA"
    mol = mn.Molecule.fetch(code).add_style("ribbon")
    geom = get_geometry_set(mol.object)
    print(geom)
    print(db.Attribute(geom.mesh.attributes["position"]))
    print(geom.pointcloud)
    print(geom.curves)
    print(geom.volume)
    print(geom.grease_pencil)
    pc = geom.instances_pointcloud()
    print(db.Attribute(pc.attributes["instance_transform"]).as_array())
    print(db.Attribute(pc.attributes[".reference_index"]).as_array())
    assert len(geom) == 1
