# Node-group asset 'Separate Polymers' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputGeometry
from .is_nucleic import IsNucleic
from .is_peptide import IsPeptide


class SeparatePolymers(AssetGeometryGroup):
    """
    Separate Polymers

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
    o.peptide : GeometrySocket
        Peptide
    o.nucleic : GeometrySocket
        Nucleic
    o.other : GeometrySocket
        Other
    """

    _name = "Separate Polymers"
    _asset_name = "Separate Polymers"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.separate_polymers"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""

    class _Outputs(SocketAccessor):
        peptide: GeometrySocket
        """Peptide"""
        nucleic: GeometrySocket
        """Nucleic"""
        other: GeometrySocket
        """Other"""

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

    def _build_group(self, tree):
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        peptide = tree.outputs.geometry("Peptide")
        nucleic = tree.outputs.geometry("Nucleic")
        other = tree.outputs.geometry("Other")

        separate_geometry = g.SeparateGeometry.point(atoms, IsPeptide().o.selection)
        separate_geometry_1 = g.SeparateGeometry.point(
            separate_geometry.o.inverted, IsNucleic().o.selection
        )

        separate_geometry >> peptide
        separate_geometry_1 >> nucleic
        separate_geometry_1.o.inverted >> other


ASSET = SeparatePolymers

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
