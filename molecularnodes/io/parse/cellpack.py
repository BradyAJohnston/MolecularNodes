from pathlib import Path

import numpy as np
import bpy

from .ensemble import Ensemble
from .bcif import BCIF
from .cif import CIF
from ..parse import molecule
from ... import blender as bl
from ... import color


class CellPack(Ensemble):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.file_type = self._file_type()
        self.data = self._read(self.file_path)
        self.array = self.data.array
        self.transformations = self.data.assemblies(as_array=True)
        self.chain_ids = self.data.chain_ids

    def create_model(
            self,
            name='CellPack',
            node_setup: bool = True,
            world_scale: float = 0.01,
            fraction: float = 1.0
    ):
        self.data_object = self._create_data_object(name=f'{name}')
        self._create_object_instances(
            name=name,
            node_setup=node_setup
        )

        self._setup_node_tree(fraction=fraction)

        return self.data_object

    def _file_type(self):
        return Path(self.file_path).suffix.strip(".")

    def _read(self, file_path):
        "Read a Cellpack File"
        suffix = Path(file_path).suffix

        if suffix in (".bin", ".bcif"):
            data = BCIF(file_path)
        elif suffix == ".cif":
            data = CIF(file_path)
        else:
            raise ValueError(f"Invalid file format: '{suffix}")

        return data

    def _create_object_instances(
            self,
            name: str = 'CellPack',
            node_setup: bool = True
    ) -> bpy.types.Collection:
        collection = bl.coll.cellpack(name)

        if self.file_type == "cif":
            array = self.array[0]
        else:
            array = self.array
        for i, chain in enumerate(np.unique(array.chain_id)):
            chain_atoms = array[array.chain_id == chain]
            model, coll_none = molecule._create_model(
                array=chain_atoms,
                name=f"{str(i).rjust(4, '0')}_{chain}",
                collection=collection
            )

            colors = np.tile(color.random_rgb(i), (len(chain_atoms), 1))
            bl.obj.set_attribute(
                model,
                name="Color",
                data=colors,
                type="FLOAT_COLOR",
                overwrite=True
            )

            if node_setup:
                bl.nodes.create_starting_node_tree(
                    model,
                    name=f"MN_pack_instance_{name}",
                    set_color=False
                )

        self.data_collection = collection

        return collection

    def _create_data_object(self, name='DataObject'):
        data_object = bl.obj.create_data_object(
            self.transformations,
            name=name,
            collection=bl.coll.mn()
        )

        data_object['chain_ids'] = self.chain_ids

        return data_object

    def _setup_node_tree(
            self,
            name='CellPack',
            fraction=1.0,
            as_points=False
    ):
        mod = bl.nodes.get_mod(self.data_object)

        group = bl.nodes.new_group(name=f"MN_ensemble_{name}", fallback=False)
        mod.node_group = group

        node_pack = bl.nodes.add_custom(
            group, 'MN_pack_instances', location=[-100, 0])
        node_pack.inputs['Collection'].default_value = self.data_collection
        node_pack.inputs['Fraction'].default_value = fraction
        node_pack.inputs['As Points'].default_value = as_points

        link = group.links.new
        link(bl.nodes.get_input(group).outputs[0], node_pack.inputs[0])
        link(node_pack.outputs[0], bl.nodes.get_output(group).inputs[0])
