# Node-group asset 'Visualize Chi Angle' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
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
from ._shared.mn_pivot_peptide import MN_pivot_peptide
from ._shared.mn_units import MNUnits
from .dihedral_chi_angle import DihedralChiAngle
from .visualize_angle import VisualizeAngle


class VisualizeChiAngle(AssetGeometryGroup):
    """
    Visualize Chi Angle

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    selection : InputBoolean
        Selection
    factor : InputFloat
        Factor
    value : InputFloat
        A value which will be scaled appropriately for the world
    radius : InputFloat
        Radius

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.selection : BooleanSocket
        Selection
    i.factor : FloatSocket
        Factor
    i.value : FloatSocket
        A value which will be scaled appropriately for the world
    i.radius : FloatSocket
        Radius

    Outputs
    -------
    o.mesh : GeometrySocket
        Mesh
    o.curve : GeometrySocket
        Curve
    """

    _name = "Visualize Chi Angle"
    _asset_name = "Visualize Chi Angle"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        selection: BooleanSocket
        """Selection"""
        factor: FloatSocket
        """Factor"""
        value: FloatSocket
        """A value which will be scaled appropriately for the world"""
        radius: FloatSocket
        """Radius"""

    class _Outputs(SocketAccessor):
        mesh: GeometrySocket
        """Mesh"""
        curve: GeometrySocket
        """Curve"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        selection: InputBoolean = True,
        factor: InputFloat = 0.5,
        value: InputFloat = 0.61,
        radius: InputFloat = 0.2,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Selection": selection,
                "Factor": factor,
                "Value": value,
                "Radius": radius,
            }
        )

    def _build_group(self, tree):
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        factor = tree.inputs.float(
            "Factor", 0.5, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        value = tree.inputs.float(
            "Value",
            0.61,
            description="A value which will be scaled appropriately for the world",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        radius = tree.inputs.float("Radius", 0.2, min_value=0.0, subtype="DISTANCE")
        mesh = tree.outputs.geometry("Mesh")
        curve = tree.outputs.geometry("Curve")

        group = DihedralChiAngle()
        capture = g.CaptureAttribute.point(geometry=geometry)
        angle = capture.items.float("Angle", group.o.angle)
        bc = capture.items.vector("BC", group.o.axis)
        output = capture.items.vector("Output", group.o.up)
        vector_math = (
            bc.output
            * g.Mix(factor_float=factor, b_float=1.0, clamp_factor=True).o.result_float
            + g.Position()
        )
        group_1 = VisualizeAngle(
            points=capture.o.geometry,
            selection=(g.EdgesOfVertex().o.total > 1)
            & (selection & MN_pivot_peptide()),
            position=vector_math,
            angle=angle.output * -1.0,
            length=MNUnits(value=value).o.angstrom,
            up=output.output,
            axis=bc.output,
            radius=radius,
        )

        group_1 >> mesh
        group_1.o.curve >> curve


ASSET = VisualizeChiAngle

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
