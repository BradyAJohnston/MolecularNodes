import numpy as np
from nodebpy.nodes.geometry import CurveToMesh, SetHandleType, SetSplineType
import molecularnodes as mn
from molecularnodes.nodes import nodes
from molecularnodes.nodes.assets import AtomsToCACurves


def test_atoms_to_ca_splines():
    """
    This test is for if the chains are still connected as continuous, even though they are
    separate by a larger distance
    """
    mol = mn.Molecule.fetch("1HQM")
    with mol.tree as tree:
        atoms, join = tree.reset()
        atca = AtomsToCACurves()
        (
            atoms
            >> atca
            >> SetSplineType(spline_type="BEZIER")
            >> SetHandleType.free()
            >> CurveToMesh()
            >> join
        )

    pos_pre = mol.named_attribute("position", evaluate=True)
    atca.i.threshold.default_value = 100
    pos_post = mol.named_attribute("position", evaluate=True)
    assert not np.allclose(pos_pre.shape, pos_post.shape)
