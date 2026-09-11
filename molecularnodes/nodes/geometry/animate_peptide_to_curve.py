# Node-group asset "Animate Peptide to Curve" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputGeometry
from ._shared.mn_utils_aa_atom_pos import MN_utils_aa_atom_pos


class MN_utils_curve_resample(CustomGeometryGroup):
    _name = ".MN_utils_curve_resample"
    _tree_properties = {"node_tool_idname": "geometry.mn_utils_curve_resample"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        offset = tree.inputs.float(
            "Offset", 2.3, min_value=-10_000.0, max_value=10_000.0
        )
        length = tree.inputs.float("Length", 0.36, min_value=0.01, subtype="DISTANCE")
        field_float = tree.inputs.float("Field Float", 0.0, hide_value=True)
        field_int = tree.inputs.integer("Field Int", 0, hide_value=True)
        field_vec = tree.inputs.vector("Field Vec", (0.0, 0.0, 0.0), hide_value=True)
        geometry_1 = tree.outputs.geometry("Geometry")
        position = tree.outputs.vector("Position")
        tangent = tree.outputs.vector("Tangent")
        normal = tree.outputs.vector("Normal")
        field_float_1 = tree.outputs.float("Field Float")
        field_int_1 = tree.outputs.integer("Field Int")
        field_vec_1 = tree.outputs.vector("Field Vec")

        capture = g.CaptureAttribute.curve(
            geometry=g.ResampleCurve(
                curve=geometry, length=length, mode="Length", count=18
            )
        )
        value = capture.items.integer("Value", g.Index())
        accumulate_field = length.point.trailing(value.output)
        switch = g.Compare.float.equal(offset, 0.0, 0.001).o.result.switch.float(
            (accumulate_field + offset).wrap(0.0, g.SplineLength().o.length),
            accumulate_field,
        )
        sample_curve = g.SampleCurve(
            curves=capture.o.geometry,
            value=field_float,
            length=switch,
            curve_index=value.output,
            mode="LENGTH",
        )
        sample_curve_1 = g.SampleCurve(
            curves=capture.o.geometry,
            value=field_int,
            length=switch,
            curve_index=value.output,
            mode="LENGTH",
            data_type="INT",
        )
        sample_curve_2 = g.SampleCurve(
            curves=capture.o.geometry,
            value=field_vec,
            length=switch,
            curve_index=value.output,
            mode="LENGTH",
            data_type="FLOAT_VECTOR",
        )
        (
            g.ResampleCurve(
                curve=capture.o.geometry,
                count=g.SplineLength().o.point_count - 1.0,
                length=0.1,
            )
            >> g.SetPosition(position=sample_curve.o.position)
            >> geometry_1
        )

        sample_curve.o.position >> position
        sample_curve.o.tangent >> tangent
        sample_curve.o.normal >> normal
        sample_curve >> field_float_1
        sample_curve_1 >> field_int_1
        sample_curve_2 >> field_vec_1


class AnimatePeptideToCurve(AssetGeometryGroup):
    """
    Animate Peptide to Curve

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    curve : InputGeometry
        Curve
    offset : InputFloat
        Offset
    start : InputFloat
        Start
    end : InputFloat
        End
    rotate : InputFloat
        Rotate
    twist : InputFloat
        Twist

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.curve : GeometrySocket
        Curve
    i.offset : FloatSocket
        Offset
    i.start : FloatSocket
        Start
    i.end : FloatSocket
        End
    i.rotate : FloatSocket
        Rotate
    i.twist : FloatSocket
        Twist

    Outputs
    -------
    o.atoms : GeometrySocket
        Atoms
    """

    _name = "Animate Peptide to Curve"
    _asset_name = "Animate Peptide to Curve"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.animate_peptide_to_curve"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        curve: GeometrySocket
        """Curve"""
        offset: FloatSocket
        """Offset"""
        start: FloatSocket
        """Start"""
        end: FloatSocket
        """End"""
        rotate: FloatSocket
        """Rotate"""
        twist: FloatSocket
        """Twist"""

    class _Outputs(SocketAccessor):
        atoms: GeometrySocket
        """Atoms"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        curve: InputGeometry = None,
        offset: InputFloat = 0.0,
        start: InputFloat = 0.0,
        end: InputFloat = 1.0,
        rotate: InputFloat = 0.5,
        twist: InputFloat = 1.0,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Curve": curve,
                "Offset": offset,
                "Start": start,
                "End": end,
                "Rotate": rotate,
                "Twist": twist,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        curve = tree.inputs.geometry("Curve")
        offset = tree.inputs.float(
            "Offset", 0.0, min_value=-10_000.0, max_value=10_000.0
        )
        start = tree.inputs.float(
            "Start", 0.0, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        end = tree.inputs.float(
            "End", 1.0, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        rotate = tree.inputs.float(
            "Rotate", 0.5, min_value=-10_000.0, max_value=10_000.0
        )
        twist = tree.inputs.float("Twist", 1.0, min_value=-10_000.0, max_value=10_000.0)
        atoms_1 = tree.outputs.geometry("Atoms")

        group = MN_utils_curve_resample(
            Geometry=g.TrimCurve(curve=curve, start=start, end=end),
            Offset=offset,
            Length=g.Math.divide(3.13, 100.0),
        )
        domain_size = g.DomainSize(geometry=group, component="CURVE")
        with g.Frame("Initial setup and alignment of amino acids"):
            group_1 = MN_utils_aa_atom_pos(atom_name=1)
            align_rotation_to_vector = g.AlignRotationToVector(
                vector=MN_utils_aa_atom_pos(atom_name=4).o.position, pivot_axis="X"
            )
            position = g.Position()
            align_rotation_to_vector_1 = g.AlignRotationToVector(
                vector=group_1.o.position
                - MN_utils_aa_atom_pos(atom_name=3).o.position,
                axis="X",
            )
            vector_rotate = g.VectorRotate(
                vector=position,
                rotation=align_rotation_to_vector_1,
                rotation_type="EULER_XYZ",
                invert=True,
            )
            vector_rotate_1 = g.VectorRotate(
                vector=position,
                rotation=align_rotation_to_vector,
                rotation_type="EULER_XYZ",
                invert=True,
            )
            set_position = (
                atoms
                >> g.SeparateGeometry.point(
                    selection=group_1.o.group_index < domain_size.o.point_count
                )
                >> g.SetPosition(offset=group_1.o.position * -1.0)
                >> g.SetPosition(position=vector_rotate)
                >> g.SetPosition(position=vector_rotate_1)
                >> g.SetPosition(
                    position=g.VectorRotate(
                        vector=position,
                        angle=-0.5574582,
                        rotation_type="Y_AXIS",
                        invert=True,
                    )
                )
            )
        sample_index = g.SampleIndex(
            geometry=group,
            value=group.o.position,
            index=group_1.o.group_index,
            data_type="FLOAT_VECTOR",
        )
        sample_index_1 = g.SampleIndex(
            geometry=group,
            value=group.o.tangent,
            index=group_1.o.group_index,
            data_type="FLOAT_VECTOR",
        )
        sample_index_2 = g.SampleIndex(
            geometry=g.SetCurveTilt(
                curve=group, tilt=(math.pi / 2 * twist).point.leading() + rotate
            ),
            value=g.CurveTilt(),
            index=group_1.o.group_index,
        )
        with g.Frame("Placing and Aligning AA Along the Curve"):
            position_1 = g.Position()
            vector_rotate_2 = g.VectorRotate.euler(
                position_1,
                group_1.o.position,
                g.AlignRotationToVector(vector=sample_index_1.o.value * -1.0, axis="X"),
            )
            vector_rotate_3 = g.VectorRotate(
                vector=position_1,
                center=group_1.o.position,
                axis=sample_index_1,
                angle=sample_index_2,
            )
            (
                set_position
                >> g.SetPosition(offset=sample_index.o.value - group_1.o.position)
                >> g.SetPosition(position=vector_rotate_2)
                >> g.SetPosition(position=vector_rotate_3)
                >> atoms_1
            )


ASSET = AnimatePeptideToCurve

ASSET_METADATA = {
    "catalog_id": "85730213-4c2e-469f-b333-52ac53adf274",
}
