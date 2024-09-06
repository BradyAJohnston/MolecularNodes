import os
import json
from pathlib import Path

import numpy as np
import bpy

from .ensemble import Ensemble
from .cif import CIF
from ..molecule import molecule
from ... import blender as bl
from ... import color


class CellPack(Ensemble):
    def __init__(self, file_path, remove_space=False):
        super().__init__(file_path)
        self.file_type = self._file_type()
        self.data = self._read(self.file_path, remove_space)
        self.array = self.data.array
        # look up color_palette of entity_id
        wpath = os.path.dirname(os.path.abspath(self.file_path))
        self.color_palette = os.path.join(wpath, "color_palette.json")
        self.color_entity = {}
        if os.path.exists(self.color_palette):
            self.color_palette = json.load(open(self.color_palette, 'r'))
        for entity in np.unique(self.array.entity_id):
            ename = self.data.entities[entity]
            if ename in self.color_palette:
                self.color_entity[entity] = np.array([
                    self.color_palette[ename]['x']/255.0,
                    self.color_palette[ename]['y']/255.0,
                    self.color_palette[ename]['z']/255.0,
                    1.0])
            else:
                self.color_entity[entity] = color.random_rgb(int(entity))
        self.transformations = self.data.assemblies(as_array=True)
        self.chain_ids = self.array.asym_id
        self.entity_ids = np.unique(self.array.entity_id)
        self.entity_chains = {}
        for i, entity in enumerate(self.entity_ids):
            symids = self.array.asym_id[self.array.entity_id == entity]
            self.entity_chains[entity] = np.unique(symids)

    def create_transparent_material(self, name="MN Transparent"):
        # Create a new material
        material_name = name
        material = bpy.data.materials.new(name=material_name)

        # Enable 'Use Nodes'
        material.use_nodes = True

        # Clear all default nodes
        nodes = material.node_tree.nodes
        nodes.clear()

        # Add a Material Output node
        output_node = nodes.new(type='ShaderNodeOutputMaterial')
        output_node.location = (300, 0)

        # Add a Transparent BSDF node
        transparent_node = nodes.new(type='ShaderNodeBsdfTransparent')
        transparent_node.location = (0, 0)

        # Connect the Transparent BSDF node to the Material Output node
        material.node_tree.links.new(transparent_node.outputs['BSDF'], output_node.inputs['Surface'])

        # Optionally set the color of the transparent BSDF
        transparent_node.inputs['Color'].default_value = (1, 1, 1, 1)  # RGBA
        return material

    def create_object(
        self,
        name="CellPack",
        node_setup: bool = True,
        world_scale: float = 0.01,
        fraction: float = 1.0,
    ):
        self.data_object = self._create_data_object(name=f"{name}")
        self._create_object_instances(name=name, node_setup=node_setup)

        self._setup_node_tree(fraction=fraction)

        return self.data_object

    def _file_type(self):
        return Path(self.file_path).suffix.strip(".")

    def _read(self, file_path, remove_space=False):
        "Read a Cellpack File"
        data = CIF(file_path, remove_space)
        return data

    def _create_object_instances(
        self, name: str = "CellPack", node_setup: bool = True
    ) -> bpy.types.Collection:
        collection = bl.coll.cellpack(name)
        for i, chain in enumerate(np.unique(self.array.asym_id)):
            # print(f"Creating chain {chain}...")
            chain_atoms = self.array[self.array.asym_id == chain]
            model, coll_none = molecule._create_object(
                array=chain_atoms,
                name=f"{str(i).rjust(4, '0')}_{chain}",
                collection=collection,
            )
            # random color per chain
            # could also do by entity, + chain-lighten + atom-lighten
            entity = chain_atoms.entity_id[0]
            color_entity = self.color_entity[entity]
            # color.random_rgb(int(entity))
            # lighten for each chain
            nc = len(self.entity_chains[entity])
            ci = np.where(self.entity_chains[entity] == chain)[0][0] * 2
            color_chain = color.Lab.lighten_color(color_entity,
                                        (float(ci) / nc))
            colors = np.tile(color_chain, (len(chain_atoms), 1))
            bl.mesh.store_named_attribute(
                model,
                name="Color",
                data=colors,
                data_type="FLOAT_COLOR",
                overwrite=True
            )

            if node_setup:
                bl.nodes.create_starting_node_tree(
                    model, name=f"MN_pack_instance_{name}", color=None
                )

        self.data_collection = collection

        return collection

    def _create_data_object(self, name="DataObject"):
        data_object = bl.mesh.create_data_object(
            self.transformations, name=name, collection=bl.coll.mn()
        )

        data_object["chain_ids"] = self.chain_ids

        return data_object

    def _setup_node_tree(self, name="CellPack", fraction=1.0, as_points=False):
        mod = bl.nodes.get_mod(self.data_object)

        group = bl.nodes.new_group(name=f"MN_ensemble_{name}", fallback=False)
        mod.node_group = group

        node_pack = bl.nodes.add_custom(
            group, 'Ensemble Instance', location=[-100, 0])
        node_pack.inputs['Instances'].default_value = self.data_collection
        # node_pack.inputs['Fraction'].default_value = fraction
        # node_pack.inputs['As Points'].default_value = as_points

        # Create the GeometryNodeIsViewport node
        node_is_viewport = group.nodes.new('GeometryNodeIsViewport')
        node_is_viewport.location = (-490.0, -240.0)

        # Create the GeometryNodeSwitch node
        node_switch = group.nodes.new('GeometryNodeSwitch')
        node_switch.location = (-303.0, -102.0)
        # Set the input type of the switch node to FLOAT
        node_switch.input_type = 'FLOAT'
        # Set the true and false values of the switch node
        node_switch.inputs[1].default_value = 1.0
        node_switch.inputs[2].default_value = 0.1

        group.links.new(node_is_viewport.outputs[0], node_switch.inputs[0])

        group.links.new(node_switch.outputs[0], node_pack.inputs['Fraction'])

        group.links.new(node_is_viewport.outputs[0], node_pack.inputs['As Points'])

        # createa a plane primitive node
        node_plane = group.nodes.new('GeometryNodeMeshGrid')
        node_plane.location = (-1173, 252)

        # create a geomtry transform node
        node_transform = group.nodes.new('GeometryNodeTransform')
        node_transform.location = (-947, 245)
        # change mesh translation
        node_transform.inputs[1].default_value = (3.0, 0.0, 0.0)
        # change mesh rotation
        node_transform.inputs[2].default_value = (0.0, 3.14/2.0, 0.0)
        # change mesh scale
        node_transform.inputs[3].default_value = (50.0, 50.0, 1.0)
        # link the plane to the transform node
        group.links.new(node_plane.outputs[0], node_transform.inputs[0])

        # create transparent material and setMaterial node
        material = self.create_transparent_material()
        node_set_material = group.nodes.new('GeometryNodeSetMaterial')
        node_set_material.location = (-100, 289)
        group.links.new(node_transform.outputs[0], node_set_material.inputs[0])
        node_set_material.inputs[2].default_value = material
        # create the join geoemtry node
        node_join = group.nodes.new('GeometryNodeJoinGeometry')
        node_join.location = (151, 122)
        group.links.new(node_set_material.outputs[0], node_join.inputs[0])
        group.links.new(node_pack.outputs[0], node_join.inputs[0])

        # create a geomtry proximity node and link the plane to it
        node_proximity = group.nodes.new('GeometryNodeProximity')
        node_proximity.location = (-586, 269)
        group.links.new(node_transform.outputs[0], node_proximity.inputs[0])

        # get the position attribute node
        node_position = group.nodes.new('GeometryNodeInputPosition')
        node_position.location = (-796, 86)

        # link it to the posistion sample in proximity
        group.links.new(node_position.outputs[0], node_proximity.inputs[2])

        # create a compare node that take the distance from the proximity node
        # and compare it to be greter than 2.0
        node_compare = group.nodes.new('FunctionNodeCompare')
        node_compare.location = (-354, 316)
        node_compare.data_type = 'FLOAT'
        node_compare.operation = 'GREATER_THAN'
        node_compare.inputs[1].default_value = 2.0
        # do the link
        group.links.new(node_proximity.outputs[1], node_compare.inputs[0])

        # link the outpot of the compare node to the selection node_pack
        group.links.new(node_compare.outputs[0], node_pack.inputs['Selection'])

        link = group.links.new
        link(bl.nodes.get_input(group).outputs[0], node_pack.inputs[0])
        link(node_join.outputs[0], bl.nodes.get_output(group).inputs[0])
