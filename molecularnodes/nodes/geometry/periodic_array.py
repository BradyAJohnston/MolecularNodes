# Node-group asset "Periodic Array" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry, InputInteger
from .lattice_grid import LatticeGrid
from .periodic_box import PeriodicBox


class PeriodicArray(AssetGeometryGroup):
    """
    Periodic Array

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    update : InputBoolean
        Update the box lengths and angles with the simulation
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
    x : InputInteger
        Number of images along the a vector
    y : InputInteger
        Number of images along the b vector
    z : InputInteger
        Number of images along the c vector

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.update : BooleanSocket
        Update the box lengths and angles with the simulation
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
    i.x : IntegerSocket
        Number of images along the a vector
    i.y : IntegerSocket
        Number of images along the b vector
    i.z : IntegerSocket
        Number of images along the c vector

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    """

    _name = "Periodic Array"
    _asset_name = "Periodic Array"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        update: BooleanSocket
        """Update the box lengths and angles with the simulation"""
        a: FloatSocket
        b: FloatSocket
        c: FloatSocket
        alpha: FloatSocket
        beta: FloatSocket
        gamma: FloatSocket
        x: IntegerSocket
        """Number of images along the a vector"""
        y: IntegerSocket
        """Number of images along the b vector"""
        z: IntegerSocket
        """Number of images along the c vector"""

    class _Outputs(SocketAccessor):
        instances: GeometrySocket
        """Instances"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        update: InputBoolean = True,
        a: InputFloat = 0.0,
        b: InputFloat = 0.0,
        c: InputFloat = 0.0,
        alpha: InputFloat = 0.0,
        beta: InputFloat = 0.0,
        gamma: InputFloat = 0.0,
        x: InputInteger = 3,
        y: InputInteger = 3,
        z: InputInteger = 3,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Update": update,
                "a": a,
                "b": b,
                "c": c,
                "alpha": alpha,
                "beta": beta,
                "gamma": gamma,
                "X": x,
                "Y": y,
                "Z": z,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        with tree.inputs.panel("Periodic Box", default_closed=True):
            update = tree.inputs.boolean(
                "Update",
                True,
                description="Update the box lengths and angles with the simulation",
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
        with tree.inputs.panel("Count"):
            x = tree.inputs.integer(
                "X",
                3,
                description="Number of images along the a vector",
                min_value=1,
                max_value=10000,
            )
            y = tree.inputs.integer(
                "Y",
                3,
                description="Number of images along the b vector",
                min_value=1,
                max_value=10000,
            )
            z = tree.inputs.integer(
                "Z",
                3,
                description="Number of images along the c vector",
                min_value=1,
                max_value=10000,
            )
        instances = tree.outputs.geometry("Instances")

        group = PeriodicBox(
            update=update, a=a, b=b, c=c_, alpha=alpha, beta=beta, gamma=gamma
        )
        (
            LatticeGrid(a=group.o.a, b=group.o.b, c=group.o.c, x=x, y=y, z=z)
            >> g.InstanceOnPoints(instance=geometry)
            >> instances
        )


ASSET = PeriodicArray

ASSET_METADATA = {
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
