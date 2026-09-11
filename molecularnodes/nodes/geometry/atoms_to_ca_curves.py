# Node-group asset "Atoms to CA Curves" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from ._shared.mn_bs_smooth import MN_bs_smooth
from ._shared.mn_init_tmp_attributes import MN_init_tmp_attributes
from ._shared.mn_topo_assign_backbone import MN_topo_assign_backbone
from .angstrom_to_world import AngstromToWorld
from .atoms_to_curves import AtomsToCurves
from .backbone_vectors import BackboneVectors
from .is_alpha_carbon import IsAlphaCarbon


class AtomsToCACurves(AssetGeometryGroup):
    """
    Atoms to CA Curves

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this node to
    bs_smoothing : InputFloat
        BS Smoothing
    threshold : InputFloat
        Distance (Angstroms) over which subsequent CA points are treated as a new chain

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.bs_smoothing : FloatSocket
        BS Smoothing
    i.threshold : FloatSocket
        Distance (Angstroms) over which subsequent CA points are treated as a new chain

    Outputs
    -------
    o.curves : GeometrySocket
        Curves
    """

    _name = "Atoms to CA Curves"
    _asset_name = "Atoms to CA Curves"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.atoms_to_ca_curves"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        bs_smoothing: FloatSocket
        """BS Smoothing"""
        threshold: FloatSocket
        """Distance (Angstroms) over which subsequent CA points are treated as a new chain"""

    class _Outputs(SocketAccessor):
        curves: GeometrySocket
        """Curves"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        bs_smoothing: InputFloat = 1.0,
        threshold: InputFloat = 4.5,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "BS Smoothing": bs_smoothing,
                "Threshold": threshold,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        bs_smoothing = tree.inputs.float(
            "BS Smoothing", 1.0, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        threshold = tree.inputs.float(
            "Threshold",
            4.5,
            description="Distance (Angstroms) over which subsequent CA points are treated as a new chain",
            min_value=0.0,
            max_value=10_000.0,
        )
        curves = tree.outputs.geometry("Curves")

        with g.Frame("Turn backbone points to curves"):
            group = AtomsToCurves(
                atoms=MN_topo_assign_backbone(atoms=atoms).o.atoms,
                selection=IsAlphaCarbon(and_=selection).o.selection,
                cutoff=threshold,
            )
        position = g.Position()
        set_curve_normal = MN_init_tmp_attributes(
            geometry=group
            >> g.StoreNamedAttribute.point.integer(name="tmp_idx", value=g.Index())
        ) >> g.SetCurveNormal(
            normal=BackboneVectors(method="Read").o.normal, mode="Free"
        )
        vector_math = position.o.position.point.at(
            g.PointsOfCurve().o.point_index
        ).distance(
            position.o.position.point.at(g.PointsOfCurve(sort_index=-1).o.point_index)
        )
        (
            MN_bs_smooth(geometry=set_curve_normal, factor=bs_smoothing, iterations=1)
            >> g.SeparateGeometry.spline(selection=g.SplineLength().o.point_count > 1)
            >> g.SetSplineCyclic(cyclic=vector_math < AngstromToWorld(angstrom=4.0))
            >> curves
        )


ASSET = AtomsToCACurves

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
