# Node-group asset "Periodic Array" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
        Geometry to copy into each periodic image
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
    x : InputInteger
        Number of images along box vector A
    y : InputInteger
        Number of images along box vector B
    z : InputInteger
        Number of images along box vector C

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to copy into each periodic image
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
    i.x : IntegerSocket
        Number of images along box vector A
    i.y : IntegerSocket
        Number of images along box vector B
    i.z : IntegerSocket
        Number of images along box vector C

    Outputs
    -------
    o.instances : GeometrySocket
        A copy of the geometry in each periodic image, starting at the original and going along box vectors A, B and C
    """

    _name = "Periodic Array"
    _asset_name = "Periodic Array"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to copy into each periodic image"""
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
        x: IntegerSocket
        """Number of images along box vector A"""
        y: IntegerSocket
        """Number of images along box vector B"""
        z: IntegerSocket
        """Number of images along box vector C"""

    class _Outputs(SocketAccessor):
        instances: GeometrySocket
        """A copy of the geometry in each periodic image, starting at the original and going along box vectors A, B and C"""

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
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to copy into each periodic image"
        )
        with tree.inputs.panel("Periodic Box", default_closed=True):
            update = tree.inputs.boolean(
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
                    "beta",
                    0.0,
                    description="Angle between box vectors A and C, in degrees",
                )
                gamma = tree.inputs.float(
                    "gamma",
                    0.0,
                    description="Angle between box vectors A and B, in degrees",
                )
        with tree.inputs.panel("Count"):
            x = tree.inputs.integer(
                "X",
                3,
                description="Number of images along box vector A",
                min_value=1,
                max_value=10000,
            )
            y = tree.inputs.integer(
                "Y",
                3,
                description="Number of images along box vector B",
                min_value=1,
                max_value=10000,
            )
            z = tree.inputs.integer(
                "Z",
                3,
                description="Number of images along box vector C",
                min_value=1,
                max_value=10000,
            )
        instances = tree.outputs.geometry(
            "Instances",
            description="A copy of the geometry in each periodic image, starting at the original and going along box vectors A, B and C",
        )

        periodic_box = PeriodicBox(
            update=update, a=a, b=b, c=c_, alpha=alpha, beta=beta, gamma=gamma
        )
        (
            LatticeGrid(
                a=periodic_box.o.a,
                b=periodic_box.o.b,
                c=periodic_box.o.c,
                x=x,
                y=y,
                z=z,
            )
            >> g.InstanceOnPoints(instance=geometry)
            >> instances
        )


ASSET = PeriodicArray

ASSET_METADATA = {
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
