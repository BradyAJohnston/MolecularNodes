from abc import ABCMeta

from ... import blender as bl

class Ensemble(metaclass=ABCMeta):
    
    def create_model(self, name='CellPack', node_setup: bool=True, world_scale: float=0.01, fraction: float=1.0):
        "Create data object and instancing collection for the Ensemble."
        data_model = self._create_data_object(name=f'{name}')
        data_model['chain_ids'] = self.chain_ids
        instance_collection = self._create_object_instances(name=name, node_setup=node_setup)
        
        self._setup_node_tree(data_model, instance_collection, fraction=fraction)
        
        return data_model

    
    def _create_data_object(self, name='DataObject'):
        "Create the data object for the Ensemble"
        data_object = bl.obj.create_data_object(self.transformations, name=name, collection=bl.coll.mn())
        
        # set chain ids and other unique object info
        return data_object
    
    def _setup_node_tree(self, model, collection, name='CellPack', fraction=1.0, as_points = False):
        # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
        mod = bl.nodes.get_mod(model)
        
        group = bl.nodes.new_group(name = f"MN_ensemble_{name}", fallback=False)
        mod.node_group = group
        
        node_pack = bl.nodes.add_custom(group, 'MN_pack_instances', location=[-100,0])
        node_pack.inputs['Collection'].default_value = collection
        node_pack.inputs['Fraction'].default_value = fraction
        node_pack.inputs['As Points'].default_value = as_points
        
        link = group.links.new
        link(bl.nodes.get_input(group).outputs[0], node_pack.inputs[0])
        link(node_pack.outputs[0], bl.nodes.get_output(group).inputs[0])