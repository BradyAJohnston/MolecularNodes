import MolecularNodes as mn
import numpy as np

def test_entity_id():
    mol = mn.load.molecule_rcsb('1cd3')
    ents = mn.obj.get_attribute(mol, 'entity_id')
    
    assert np.all(np.isin(ents, np.array(range(4))))
    