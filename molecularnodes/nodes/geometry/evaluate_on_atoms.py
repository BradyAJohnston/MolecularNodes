# Node-group asset "Evaluate on Atoms" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
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
    result : InputMenu | Literal["Atoms", "Representation"]
        What the closure returns. 'Atoms': the closure runs on the input geometry and its result is the output. 'Representation' (styles): the closure runs on the `MN/Atoms` bundle if the input has one, its result is joined to the input's existing representation, and the atoms are stored in the `MN/Atoms` bundle so later styles can still get them

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to get the bundle of
    i.selection : BooleanSocket
        The parts of the geometry that go into the first output
    i.closure : ClosureSocket
        Closure
    i.result : MenuSocket
        What the closure returns. 'Atoms': the closure runs on the input geometry and its result is the output. 'Representation' (styles): the closure runs on the `MN/Atoms` bundle if the input has one, its result is joined to the input's existing representation, and the atoms are stored in the `MN/Atoms` bundle so later styles can still get them

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Evaluate on Atoms"
    _asset_name = "Evaluate on Atoms"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to get the bundle of"""
        selection: BooleanSocket
        """The parts of the geometry that go into the first output"""
        closure: ClosureSocket
        """Closure"""
        result: MenuSocket
        """What the closure returns. 'Atoms': the closure runs on the input geometry and its result is the output. 'Representation' (styles): the closure runs on the `MN/Atoms` bundle if the input has one, its result is joined to the input's existing representation, and the atoms are stored in the `MN/Atoms` bundle so later styles can still get them"""

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
        result: InputMenu | Literal["Atoms", "Representation"] = "Atoms",
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
            "Atoms",
            description="What the closure returns. 'Atoms': the closure runs on the input geometry and its result is the output. 'Representation' (styles): the closure runs on the `MN/Atoms` bundle if the input has one, its result is joined to the input's existing representation, and the atoms are stored in the `MN/Atoms` bundle so later styles can still get them",
            expanded=True,
            optional_label=True,
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        get_geometry_atoms = GetGeometryAtoms(geometry=geometry)
        menu_switch = g.MenuSwitch.integer(result, {"Atoms": 0, "Representation": 1})
        with g.Frame("Evaluate the closure on the 'Atoms'"):
            check_geometry = CheckGeometry(
                geometry=g.IndexSwitch.geometry(
                    menu_switch.o.output, (geometry, get_geometry_atoms.o.atoms)
                ),
                message="'Atoms' contains no geometry, check your node connections.",
            )
            evaluate_closure = g.EvaluateClosure(closure)
            evaluate_closure.inputs.geometry(
                (
                    MN_ensure_ures_id(input=check_geometry)
                    >> g.SeparateGeometry.point(selection=selection)
                ).o.selection,
                "Atoms",
            )
            geometry_2 = evaluate_closure.outputs.geometry("Geometry")
        store_bundle_item = g.StoreBundleItem.geometry(
            get_geometry_atoms.o.bundle, "MN/Atoms", get_geometry_atoms.o.atoms
        )
        set_geometry_bundle = g.JoinGeometry(
            (get_geometry_atoms.o.geometry, geometry_2.output)
        ).o.geometry.set_bundle(store_bundle_item)
        (
            g.IndexSwitch.geometry(
                menu_switch.o.output, (geometry_2.output, set_geometry_bundle)
            )
            >> geometry_1
        )
        _string = g.String(
            "'Atoms' result: the closure runs on the input geometry and its result is the output.\n\n'Representation' result (styles): the 'Atoms' are the input's 'MN/Atoms' bundle, or the input geometry when it has none. The result is joined to the input's representation and the 'Atoms' are stored in the 'MN/Atoms' bundle, so styles can be stacked."
        )


ASSET = EvaluateOnAtoms

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
