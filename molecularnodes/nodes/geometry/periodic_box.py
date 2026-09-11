# Node-group asset "Periodic Box" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputFloat
from .angstrom_to_world import AngstromToWorld


class PeriodicBox(AssetGeometryGroup):
    """
    Periodic Box

    Parameters
    ----------
    update : InputBoolean
        Update the box lengths and angles with the simulation. This does not change anything directly inside of the node tree, but the mda.Universe that updates the positions will also update this node.
    a : InputFloat
        a
    b : InputFloat
        b
    c : InputFloat
        c
    alpha : InputFloat
        alpha
    beta : InputFloat
        beta
    gamma : InputFloat
        gamma

    Inputs
    ------
    i.update : BooleanSocket
        Update the box lengths and angles with the simulation. This does not change anything directly inside of the node tree, but the mda.Universe that updates the positions will also update this node.
    i.a : FloatSocket
        a
    i.b : FloatSocket
        b
    i.c : FloatSocket
        c
    i.alpha : FloatSocket
        alpha
    i.beta : FloatSocket
        beta
    i.gamma : FloatSocket
        gamma

    Outputs
    -------
    o.a : VectorSocket
        A
    o.b : VectorSocket
        B
    o.c : VectorSocket
        C
    """

    _name = "Periodic Box"
    _asset_name = "Periodic Box"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "VECTOR"

    class _Inputs(SocketAccessor):
        update: BooleanSocket
        """Update the box lengths and angles with the simulation. This does not change anything directly inside of the node tree, but the mda.Universe that updates the positions will also update this node."""
        a: FloatSocket
        b: FloatSocket
        c: FloatSocket
        alpha: FloatSocket
        beta: FloatSocket
        gamma: FloatSocket

    class _Outputs(SocketAccessor):
        a: VectorSocket
        """A"""
        b: VectorSocket
        """B"""
        c: VectorSocket
        """C"""

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
            description="Update the box lengths and angles with the simulation. This does not change anything directly inside of the node tree, but the mda.Universe that updates the positions will also update this node.",
            default_attribute="True",
        )
        with tree.inputs.panel("Lengths"):
            a = tree.inputs.float("a", 0.0)
            b = tree.inputs.float("b", 0.0)
            c_ = tree.inputs.float("c", 0.0)
        with tree.inputs.panel("Angles"):
            alpha = tree.inputs.float("alpha", 0.0)
            beta = tree.inputs.float("beta", 0.0)
            gamma = tree.inputs.float("gamma", 0.0)
        a_1 = tree.outputs.vector("A")
        b_1 = tree.outputs.vector("B")
        c__1 = tree.outputs.vector("C")

        with g.Frame("A Vector"):
            combine_xyz = g.CombineXYZ(x=AngstromToWorld(angstrom=a))
        with g.Frame("Angles To Radians"):
            math_1 = alpha.to_radians()
            math_2 = beta.to_radians()
            math_3 = gamma.to_radians()
        with g.Frame("Trig Ops"):
            math_4 = math_3.cos()
            math_5 = math_2.cos()
            math_6 = math_3.sin()
            math_7 = math_1.cos()
        with g.Frame("B Vector"):
            group = AngstromToWorld(angstrom=b)
            with g.Frame("B Y Component"):
                math_8 = group.o.world * math_6
            with g.Frame("B X Component"):
                math_9 = group.o.world * math_4
            combine_xyz_1 = g.CombineXYZ(x=math_9, y=math_8)
        with g.Frame("C Vector"):
            group_1 = AngstromToWorld(angstrom=c_)
            with g.Frame("C X Component"):
                math_10 = group_1.o.world * math_5
            with g.Frame("C Y Component"):
                math_11 = group_1.o.world * ((math_7 - math_5 * math_4) / math_6)
            with g.Frame("C Z Component"):
                math_12 = (group_1.o.world**2.0 - math_10**2.0 - math_11**2.0).sqrt()
            combine_xyz_2 = g.CombineXYZ(x=math_10, y=math_11, z=math_12)

        combine_xyz >> a_1
        combine_xyz_1 >> b_1
        combine_xyz_2 >> c__1


ASSET = PeriodicBox

ASSET_METADATA = {
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
