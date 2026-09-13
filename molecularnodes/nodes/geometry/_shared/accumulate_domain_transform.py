# Node group ".Accumulate Domain Transform" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    IntegerSocket,
    MatrixSocket,
    MenuSocket,
    SocketAccessor,
)
from nodebpy.types import InputInteger, InputMatrix, InputMenu
from .domain_switch_matrix import DomainSwitchMatrix


class AccumulateDomainTransform(CustomGeometryGroup):
    """
    .Accumulate Domain Transform

    Parameters
    ----------
    domain : InputMenu | Literal["Point", "Edge", "Face", "Face Corner", "Spline", "Instance"]
        Domain
    transform : InputMatrix
        Transform
    group_id : InputInteger
        Group ID

    Inputs
    ------
    i.domain : MenuSocket
        Domain
    i.transform : MatrixSocket
        Transform
    i.group_id : IntegerSocket
        Group ID

    Outputs
    -------
    o.leading : MatrixSocket
        Leading
    o.trailling : MatrixSocket
        Trailling
    """

    _name = ".Accumulate Domain Transform"
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry._accumulate_domain_transform"}

    class _Inputs(SocketAccessor):
        domain: MenuSocket
        """Domain"""
        transform: MatrixSocket
        """Transform"""
        group_id: IntegerSocket
        """Group ID"""

    class _Outputs(SocketAccessor):
        leading: MatrixSocket
        """Leading"""
        trailling: MatrixSocket
        """Trailling"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        domain: InputMenu
        | Literal[
            "Point", "Edge", "Face", "Face Corner", "Spline", "Instance"
        ] = "Point",
        transform: InputMatrix = None,
        group_id: InputInteger = 0,
    ):
        super().__init__(
            **{"Domain": domain, "Transform": transform, "Group ID": group_id}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        domain = tree.inputs.menu("Domain", optional_label=True, hide_value=True)
        transform = tree.inputs.matrix("Transform", hide_value=True)
        group_id = tree.inputs.integer("Group ID", 0, hide_value=True)
        leading = tree.outputs.matrix("Leading")
        trailling = tree.outputs.matrix("Trailling")

        accumulate_field = g.AccumulateField.point.transform(transform, group_id)
        accumulate_field_1 = g.AccumulateField.edge.transform(transform, group_id)
        accumulate_field_2 = g.AccumulateField.face.transform(transform, group_id)
        accumulate_field_3 = g.AccumulateField.corner.transform(transform, group_id)
        accumulate_field_4 = g.AccumulateField.spline.transform(transform, group_id)
        accumulate_field_5 = g.AccumulateField.instance.transform(transform, group_id)
        (
            DomainSwitchMatrix(
                menu=domain,
                point=accumulate_field.o.leading,
                edge=accumulate_field_1.o.leading,
                face=accumulate_field_2.o.leading,
                face_corner=accumulate_field_3.o.leading,
                spline=accumulate_field_4.o.leading,
                instance=accumulate_field_5.o.leading,
            )
            >> leading
        )
        (
            DomainSwitchMatrix(
                menu=domain,
                point=accumulate_field.o.trailing,
                edge=accumulate_field_1.o.trailing,
                face=accumulate_field_2.o.trailing,
                face_corner=accumulate_field_3.o.trailing,
                spline=accumulate_field_4.o.trailing,
                instance=accumulate_field_5.o.trailing,
            )
            >> trailling
        )

        domain.default_value = "Point"
