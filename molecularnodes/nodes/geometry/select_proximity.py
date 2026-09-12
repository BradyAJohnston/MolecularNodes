# Node-group asset "Select Proximity" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry
from ._shared.mn_units import MNUnits
from .select_res_whole import SelectResWhole


class SelectProximity(AssetGeometryGroup):
    """
    Select Proximity

    Parameters
    ----------
    target_atoms : InputGeometry
        The atoms to measure the distance from.
    subset : InputBoolean
        Subset of input atoms to use for proximity calculation
    expand : InputBoolean
        Include an entire residue if even a single atom is within the threshold
    distance_a : InputFloat
        Cutoff distance for the selection in Angstroms

    Inputs
    ------
    i.target_atoms : GeometrySocket
        The atoms to measure the distance from.
    i.subset : BooleanSocket
        Subset of input atoms to use for proximity calculation
    i.expand : BooleanSocket
        Include an entire residue if even a single atom is within the threshold
    i.distance_a : FloatSocket
        Cutoff distance for the selection in Angstroms

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    o.inverted : BooleanSocket
        The inverse of the calculated selection
    """

    _name = "Select Proximity"
    _asset_name = "Select Proximity"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_proximity"}

    class _Inputs(SocketAccessor):
        target_atoms: GeometrySocket
        """The atoms to measure the distance from."""
        subset: BooleanSocket
        """Subset of input atoms to use for proximity calculation"""
        expand: BooleanSocket
        """Include an entire residue if even a single atom is within the threshold"""
        distance_a: FloatSocket
        """Cutoff distance for the selection in Angstroms"""

    class _Outputs(SocketAccessor):
        selection: BooleanSocket
        """The calculated selection"""
        inverted: BooleanSocket
        """The inverse of the calculated selection"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        target_atoms: InputGeometry = None,
        subset: InputBoolean = True,
        expand: InputBoolean = False,
        distance_a: InputFloat = 5.0,
    ):
        super().__init__(
            **{
                "Target Atoms": target_atoms,
                "Subset": subset,
                "Expand": expand,
                "Distance (A)": distance_a,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        target_atoms = tree.inputs.geometry(
            "Target Atoms", description="The atoms to measure the distance from."
        )
        subset = tree.inputs.boolean(
            "Subset",
            True,
            description="Subset of input atoms to use for proximity calculation",
            hide_value=True,
        )
        expand = tree.inputs.boolean(
            "Expand",
            False,
            description="Include an entire residue if even a single atom is within the threshold",
        )
        distance_a = tree.inputs.float(
            "Distance (A)",
            5.0,
            description="Cutoff distance for the selection in Angstroms",
            min_value=0.0,
            max_value=10_000.0,
        )
        selection = tree.outputs.boolean(
            "Selection", description="The calculated selection"
        )
        inverted = tree.outputs.boolean(
            "Inverted", description="The inverse of the calculated selection"
        )

        geometry_proximity = (
            target_atoms
            >> g.SeparateGeometry.point(selection=subset)
            >> g.GeometryProximity(target_element="POINTS")
        )
        boolean_math = g.BooleanMath.subtract(
            geometry_proximity.o.distance < MNUnits(value=distance_a).o.angstrom,
            g.Switch.boolean(
                g.AccumulateField.point.integer(~subset).o.total, true=subset
            ),
        )
        group = SelectResWhole(selection=boolean_math, expand=expand)
        ~group.o.selection >> inverted

        group >> selection


ASSET = SelectProximity

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
