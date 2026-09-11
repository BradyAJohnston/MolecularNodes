# Node-group asset "Sample Nearest Atoms" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputGeometry
from .atomic_number import AtomicNumber
from .chain_id import ChainID
from .color import Color
from .residue_id import ResidueID
from .residue_name import ResidueName


class SampleNearestAtoms(AssetGeometryGroup):
    """
    Sample Nearest Atoms

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges

    Outputs
    -------
    o.color : ColorSocket
        Color
    o.b_factor : FloatSocket
        b_factor
    o.atomic_number : IntegerSocket
        atomic_number
    o.chain_id : IntegerSocket
        chain_id
    o.res_id : IntegerSocket
        res_id
    o.res_name : IntegerSocket
        res_name
    """

    _name = "Sample Nearest Atoms"
    _asset_name = "Sample Nearest Atoms"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.sample_nearest_atoms"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """Color"""
        b_factor: FloatSocket
        atomic_number: IntegerSocket
        chain_id: IntegerSocket
        res_id: IntegerSocket
        res_name: IntegerSocket

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
    ):
        super().__init__(**{"Atoms": atoms})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        color = tree.outputs.color("Color", (0.0, 0.0, 0.0, 0.0))
        b_factor = tree.outputs.float("b_factor")
        atomic_number = tree.outputs.integer("atomic_number")
        chain_id = tree.outputs.integer("chain_id")
        res_id = tree.outputs.integer("res_id")
        res_name = tree.outputs.integer("res_name")

        sample_nearest = g.SampleNearest.point(atoms)
        sample_index = g.SampleIndex(
            geometry=atoms,
            value=g.NamedAttribute.float("b_factor").o.attribute,
            index=sample_nearest,
        )
        sample_index_1 = g.SampleIndex(
            geometry=atoms, value=AtomicNumber(), index=sample_nearest, data_type="INT"
        )
        sample_index_2 = g.SampleIndex(
            geometry=atoms, value=ChainID(), index=sample_nearest, data_type="INT"
        )
        sample_index_3 = g.SampleIndex(
            geometry=atoms, value=ResidueID(), index=sample_nearest, data_type="INT"
        )
        sample_index_4 = g.SampleIndex(
            geometry=atoms, value=ResidueName(), index=sample_nearest, data_type="INT"
        )
        sample_index_5 = g.SampleIndex(
            geometry=atoms, value=Color(), index=sample_nearest, data_type="FLOAT_COLOR"
        )

        sample_index_5 >> color
        sample_index >> b_factor
        sample_index_1 >> atomic_number
        sample_index_2 >> chain_id
        sample_index_3 >> res_id
        sample_index_4 >> res_name


ASSET = SampleNearestAtoms

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
