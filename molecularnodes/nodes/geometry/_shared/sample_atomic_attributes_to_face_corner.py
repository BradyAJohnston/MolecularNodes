# Node group ".Sample Atomic Attributes to Face Corner" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    GeometrySocket,
    IntegerSocket,
    SocketAccessor,
)
from nodebpy.types import InputGeometry, InputInteger
from ..color import Color


class SampleAtomicAttributesToFaceCorner(CustomGeometryGroup):
    """
    .Sample Atomic Attributes to Face Corner

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    sample_atoms : InputGeometry
        Sample Atoms
    index : InputInteger
        Index

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.sample_atoms : GeometrySocket
        Sample Atoms
    i.index : IntegerSocket
        Index

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = ".Sample Atomic Attributes to Face Corner"
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        sample_atoms: GeometrySocket
        """Sample Atoms"""
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        sample_atoms: InputGeometry = None,
        index: InputInteger = 0,
    ):
        super().__init__(
            **{"Geometry": geometry, "Sample Atoms": sample_atoms, "Index": index}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        sample_atoms = tree.inputs.geometry("Sample Atoms")
        index = tree.inputs.integer("Index", 0, default_input="INDEX")
        geometry_1 = tree.outputs.geometry("Geometry")

        repeat_zone = g.RepeatZone(9)
        geometry_2 = repeat_zone.items.geometry("Geometry", geometry)
        index_switch = g.IndexSwitch.string(
            repeat_zone.iteration,
            (
                "atomic_number",
                "atom_id",
                "atom_name",
                "res_id",
                "res_name",
                "chain_id",
                "entity_id",
                "sec_struct",
                "ures_id",
            ),
        )
        sample_index = g.SampleIndex(
            geometry=sample_atoms,
            value=g.NamedAttribute.integer(index_switch).o.attribute,
            index=index,
            data_type="INT",
        )
        sample_index_1 = g.SampleIndex(
            geometry=sample_atoms, value=Color(), index=index, data_type="FLOAT_COLOR"
        )
        store_named_attribute = g.StoreNamedAttribute.corner.integer(
            geometry_2.current, name=index_switch, value=sample_index
        )
        store_named_attribute >> geometry_2.next
        repeat_zone_1 = g.RepeatZone(3)
        geometry_3 = repeat_zone_1.items.geometry("Geometry", geometry_2.result)
        index_switch_1 = g.IndexSwitch.string(
            repeat_zone_1.iteration, ("vdw_radii", "mass", "b_factor")
        )
        sample_index_2 = g.SampleIndex(
            geometry=sample_atoms,
            value=g.NamedAttribute.float(index_switch_1).o.attribute,
            index=index,
        )
        store_named_attribute_1 = g.StoreNamedAttribute.corner.float(
            geometry_3.current, name=index_switch_1, value=sample_index_2
        )
        store_named_attribute_1 >> geometry_3.next
        repeat_zone_2 = g.RepeatZone(6)
        geometry_4 = repeat_zone_2.items.geometry("Geometry", geometry_3.result)
        index_switch_2 = g.IndexSwitch.string(
            repeat_zone_2.iteration,
            (
                "is_alpha_carbon",
                "is_side_chain",
                "is_backbone",
                "is_solvent",
                "is_nucleic",
                "is_peptide",
            ),
        )
        sample_index_3 = g.SampleIndex(
            geometry=sample_atoms,
            value=g.NamedAttribute.boolean(index_switch_2).o.attribute,
            index=index,
            data_type="BOOLEAN",
        )
        store_named_attribute_2 = g.StoreNamedAttribute.corner.boolean(
            geometry_4.current, name=index_switch_2, value=sample_index_3
        )
        store_named_attribute_2 >> geometry_4.next
        (
            geometry_4.result
            >> g.StoreNamedAttribute.corner.color(name="Color", value=sample_index_1)
            >> geometry_1
        )
