import MolecularNodes as mn
import numpy as np

def test_entity_id():
    mol = mn.load.molecule_rcsb('1cd3')
    ents = mn.obj.get_attribute(mol, 'entity_id')
    
    assert np.all(np.isin(ents, np.array(range(4))))

def test_entity_id_range():
    mol = mn.load.molecule_rcsb('6n2y')
    ents = mn.obj.get_attribute(mol, 'entity_id')
    chain_idx = mn.obj.get_attribute(mol, 'chain_id')
    chain_id = np.array(mol['chain_id_unique'])[chain_idx]
    
    assert max(ents) == 8
    assert np.all(np.array(['A', 'B', 'C']) == np.unique(chain_id[ents == 0]))