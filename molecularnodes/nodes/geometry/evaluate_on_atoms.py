# Node-group asset "Evaluate on Atoms" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    ClosureSocket,
    GeometrySocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputClosure, InputGeometry, InputMenu
from ._shared.mn_ensure_ures_id import MN_ensure_ures_id
from .check_geometry import CheckGeometry
from .get_geometry_atoms import GetGeometryAtoms


class EvaluateOnAtoms(AssetGeometryGroup):
    """
    Evaluate on Atoms

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to get the bundle of
    selection : InputBoolean
        The parts of the geometry that go into the first output
    closure : InputClosure
        Closure
    result : InputMenu | Literal["Geometry", "Bundle"]
        Where to store the result of the closure. Bundle overwrite the existing `MN/Atoms` bundle. Geometry passes along the bundle and joins the resulting geometry into the output.

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to get the bundle of
    i.selection : BooleanSocket
        The parts of the geometry that go into the first output
    i.closure : ClosureSocket
        Closure
    i.result : MenuSocket
        Where to store the result of the closure. Bundle overwrite the existing `MN/Atoms` bundle. Geometry passes along the bundle and joins the resulting geometry into the output.

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Evaluate on Atoms"
    _asset_name = "Evaluate on Atoms"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to get the bundle of"""
        selection: BooleanSocket
        """The parts of the geometry that go into the first output"""
        closure: ClosureSocket
        """Closure"""
        result: MenuSocket
        """Where to store the result of the closure. Bundle overwrite the existing `MN/Atoms` bundle. Geometry passes along the bundle and joins the resulting geometry into the output."""

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
        selection: InputBoolean = True,
        closure: InputClosure = None,
        result: InputMenu | Literal["Geometry", "Bundle"] = "Geometry",
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Selection": selection,
                "Closure": closure,
                "Result": result,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to get the bundle of"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="The parts of the geometry that go into the first output",
            hide_value=True,
        )
        closure = tree.inputs.closure("Closure")
        result = tree.inputs.menu(
            "Result",
            description="Where to store the result of the closure. Bundle overwrite the existing `MN/Atoms` bundle. Geometry passes along the bundle and joins the resulting geometry into the output.",
            expanded=True,
            optional_label=True,
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        group = GetGeometryAtoms(geometry=geometry)
        with g.Frame("Evaluate the closure on the 'Atoms'"):
            group_1 = CheckGeometry(
                geometry=group.o.atoms,
                message="'Atoms' contains no geometry, check your node connections.",
            )
            evaluate_closure = g.EvaluateClosure(closure)
            evaluate_closure.inputs.geometry(
                "Atoms",
                (
                    MN_ensure_ures_id(input=group_1)
                    >> g.SeparateGeometry.point(selection=selection)
                ).o.selection,
            )
            geometry_2 = evaluate_closure.outputs.geometry("Geometry")
        menu_switch = g.MenuSwitch.integer(result, {"Geometry": 0, "Bundle": 1})
        store_bundle_item = g.StoreBundleItem.geometry(
            group.o.bundle,
            "MN/Atoms",
            g.IndexSwitch.geometry(menu_switch.o.output, (group.o.atoms, geometry_2)),
        )
        (
            g.JoinGeometry(
                geometry=(
                    group.o.geometry,
                    g.IndexSwitch.geometry(menu_switch.o.output, (geometry_2, None)),
                )
            )
            >> g.SetGeometryBundle(bundle=store_bundle_item)
            >> geometry_1
        )
        _string = g.String(
            string="Evaluate the 'Closure' on the 'Atoms' geometry, join the resulting geometry with the input geometry.\n\nIf the input 'Geometry' doesn't contain a 'MN/Atoms' bundle then we treat the input geometry as the 'Atoms'."
        )

        result.default_value = "Geometry"


ASSET = EvaluateOnAtoms

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
