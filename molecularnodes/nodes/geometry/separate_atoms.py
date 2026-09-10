# Node-group asset 'Separate Atoms' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputGeometry
from .evaluate_on_atoms import EvaluateOnAtoms


class SeparateAtoms(AssetGeometryGroup):
    """
    Separate Atoms

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection field for which atoms to separate

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection field for which atoms to separate

    Outputs
    -------
    o.atoms : GeometrySocket
        The selected atoms
    o.inverted : GeometrySocket
        The parts of the geometry not in the selection
    o.index : IntegerSocket
        Index of the point before being separated.
    """

    _name = "Separate Atoms"
    _asset_name = "Separate Atoms"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.separate_atoms"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection field for which atoms to separate"""

    class _Outputs(SocketAccessor):
        atoms: GeometrySocket
        """The selected atoms"""
        inverted: GeometrySocket
        """The parts of the geometry not in the selection"""
        index: IntegerSocket
        """Index of the point before being separated."""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
    ):
        super().__init__(**{"Atoms": atoms, "Selection": selection})

    def _build_group(self, tree):
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection field for which atoms to separate",
            hide_value=True,
        )
        atoms_1 = tree.outputs.geometry("Atoms", description="The selected atoms")
        inverted = tree.outputs.geometry(
            "Inverted", description="The parts of the geometry not in the selection"
        )
        index = tree.outputs.integer(
            "Index", description="Index of the point before being separated."
        )

        closure_zone = g.ClosureZone()
        atoms_2 = closure_zone.inputs.geometry("Atoms")
        geometry = closure_zone.outputs.geometry("Geometry")
        capture = g.CaptureAttribute.point(geometry=atoms_2)
        capture.items.integer("Value")
        capture.o.geometry >> geometry
        _group = EvaluateOnAtoms()
        capture_1 = g.CaptureAttribute.point(geometry=atoms)
        value = capture_1.items.integer("Value", g.Index())
        separate_geometry = capture_1.o.geometry >> g.SeparateGeometry.point(
            selection=selection
        )

        separate_geometry >> atoms_1
        separate_geometry.o.inverted >> inverted
        value.output >> index


ASSET = SeparateAtoms

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
