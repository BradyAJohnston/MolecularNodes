# Node-group asset "Periodic Box" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    MatrixSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputFloat


class PeriodicBox(AssetGeometryGroup):
    """
    Periodic Box

    Parameters
    ----------
    update : InputBoolean
        Set the box lengths and angles from the trajectory on every frame change, as the molecule updates its positions. Turn off to set them by hand. As a modifier input, reads the "True" attribute by default.
    a : InputFloat
        Length of box vector A in world units, where 1 m is 1 nm
    b : InputFloat
        Length of box vector B in world units, where 1 m is 1 nm
    c : InputFloat
        Length of box vector C in world units, where 1 m is 1 nm
    alpha : InputFloat
        Angle between box vectors B and C, in degrees
    beta : InputFloat
        Angle between box vectors A and C, in degrees
    gamma : InputFloat
        Angle between box vectors A and B, in degrees

    Inputs
    ------
    i.update : BooleanSocket
        Set the box lengths and angles from the trajectory on every frame change, as the molecule updates its positions. Turn off to set them by hand
    i.a : FloatSocket
        Length of box vector A in world units, where 1 m is 1 nm
    i.b : FloatSocket
        Length of box vector B in world units, where 1 m is 1 nm
    i.c : FloatSocket
        Length of box vector C in world units, where 1 m is 1 nm
    i.alpha : FloatSocket
        Angle between box vectors B and C, in degrees
    i.beta : FloatSocket
        Angle between box vectors A and C, in degrees
    i.gamma : FloatSocket
        Angle between box vectors A and B, in degrees

    Outputs
    -------
    o.a : VectorSocket
        Box vector A in world units, along the X axis
    o.b : VectorSocket
        Box vector B in world units, in the XY plane
    o.c : VectorSocket
        Box vector C in world units
    o.matrix : MatrixSocket
        The box vectors A, B and C as columns, which transforms fractional coordinates of the box into world coordinates
    """

    _name = "Periodic Box"
    _asset_name = "Periodic Box"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "VECTOR"

    class _Inputs(SocketAccessor):
        update: BooleanSocket
        """Set the box lengths and angles from the trajectory on every frame change, as the molecule updates its positions. Turn off to set them by hand"""
        a: FloatSocket
        """Length of box vector A in world units, where 1 m is 1 nm"""
        b: FloatSocket
        """Length of box vector B in world units, where 1 m is 1 nm"""
        c: FloatSocket
        """Length of box vector C in world units, where 1 m is 1 nm"""
        alpha: FloatSocket
        """Angle between box vectors B and C, in degrees"""
        beta: FloatSocket
        """Angle between box vectors A and C, in degrees"""
        gamma: FloatSocket
        """Angle between box vectors A and B, in degrees"""

    class _Outputs(SocketAccessor):
        a: VectorSocket
        """Box vector A in world units, along the X axis"""
        b: VectorSocket
        """Box vector B in world units, in the XY plane"""
        c: VectorSocket
        """Box vector C in world units"""
        matrix: MatrixSocket
        """The box vectors A, B and C as columns, which transforms fractional coordinates of the box into world coordinates"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        update: InputBoolean = True,
        a: InputFloat = 0.0,
        b: InputFloat = 0.0,
        c: InputFloat = 0.0,
        alpha: InputFloat = 0.0,
        beta: InputFloat = 0.0,
        gamma: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Update": update,
                "a": a,
                "b": b,
                "c": c,
                "alpha": alpha,
                "beta": beta,
                "gamma": gamma,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        _update = tree.inputs.boolean(
            "Update",
            True,
            description="Set the box lengths and angles from the trajectory on every frame change, as the molecule updates its positions. Turn off to set them by hand",
            default_attribute="True",
        )
        with tree.inputs.panel("Lengths"):
            a = tree.inputs.float(
                "a",
                0.0,
                description="Length of box vector A in world units, where 1 m is 1 nm",
                subtype="DISTANCE",
            )
            b = tree.inputs.float(
                "b",
                0.0,
                description="Length of box vector B in world units, where 1 m is 1 nm",
                subtype="DISTANCE",
            )
            c_ = tree.inputs.float(
                "c",
                0.0,
                description="Length of box vector C in world units, where 1 m is 1 nm",
                subtype="DISTANCE",
            )
        with tree.inputs.panel("Angles"):
            alpha = tree.inputs.float(
                "alpha",
                0.0,
                description="Angle between box vectors B and C, in degrees",
            )
            beta = tree.inputs.float(
                "beta", 0.0, description="Angle between box vectors A and C, in degrees"
            )
            gamma = tree.inputs.float(
                "gamma",
                0.0,
                description="Angle between box vectors A and B, in degrees",
            )
        a_1 = tree.outputs.vector(
            "A", description="Box vector A in world units, along the X axis"
        )
        b_1 = tree.outputs.vector(
            "B", description="Box vector B in world units, in the XY plane"
        )
        c__1 = tree.outputs.vector("C", description="Box vector C in world units")
        matrix = tree.outputs.matrix(
            "Matrix",
            description="The box vectors A, B and C as columns, which transforms fractional coordinates of the box into world coordinates",
        )

        with g.Frame("A Vector"):
            combine_xyz = g.CombineXYZ(a)
        with g.Frame("Angles To Radians"):
            math_1 = alpha.to_radians()
            math_2 = beta.to_radians()
            math_3 = gamma.to_radians()
        with g.Frame("Trig Ops"):
            math_4 = math_1.cos()
            math_5 = math_2.cos()
            math_6 = math_3.cos()
            math_7 = math_3.sin()
        with g.Frame("B Vector"):
            with g.Frame("B Y Component"):
                math_8 = b * math_7
            with g.Frame("B X Component"):
                math_9 = b * math_6
            combine_xyz_1 = g.CombineXYZ(x=math_9, y=math_8)
        with g.Frame("C Vector"):
            with g.Frame("C X Component"):
                math_10 = c_ * math_5
            with g.Frame("C Y Component"):
                math_11 = c_ * ((math_4 - math_5 * math_6) / math_7)
            with g.Frame("C Z Component"):
                math_12 = (c_**2.0 - math_10**2.0 - math_11**2.0).sqrt()
            combine_xyz_2 = g.CombineXYZ(x=math_10, y=math_11, z=math_12)
        combine_matrix = g.CombineMatrix(
            column_1_row_1=a,
            column_2_row_1=math_9,
            column_2_row_2=math_8,
            column_3_row_1=math_10,
            column_3_row_2=math_11,
            column_3_row_3=math_12,
        )

        combine_xyz >> a_1
        combine_xyz_1 >> b_1
        combine_xyz_2 >> c__1
        combine_matrix >> matrix


ASSET = PeriodicBox

ASSET_METADATA = {
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
