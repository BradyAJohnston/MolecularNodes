# Node-group asset 'Set Color' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    ColorSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputColor, InputGeometry


class SetColor(AssetGeometryGroup):
    """
    Set Color

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this node to
    color : InputColor
        Color to apply to the selected atoms

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.color : ColorSocket
        Color to apply to the selected atoms

    Outputs
    -------
    o.atoms : GeometrySocket
        Atomic geometry with an updated `Color` attribute
    """

    _name = "Set Color"
    _asset_name = "Set Color"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.set_color"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        color: ColorSocket
        """Color to apply to the selected atoms"""

    class _Outputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry with an updated `Color` attribute"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        color: InputColor = None,
    ):
        super().__init__(**{"Atoms": atoms, "Selection": selection, "Color": color})

    def _build_group(self, tree):
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        color = tree.inputs.color(
            "Color",
            (0.161517, 0.623961, 0.195602, 1.0),
            description="Color to apply to the selected atoms",
        )
        atoms_1 = tree.outputs.geometry(
            "Atoms", description="Atomic geometry with an updated `Color` attribute"
        )

        (
            atoms
            >> g.StoreNamedAttribute.point.color(
                selection=selection, name="Color", value=color
            )
            >> atoms_1
        )


ASSET = SetColor

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
