# Node-group asset "Symmetry Cubic" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    GeometrySocket,
    MenuSocket,
    PackageLibrary,
    RotationSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    InputFloat,
    InputGeometry,
    InputMenu,
    InputRotation,
    InputVector,
)
from ._shared.symmetry_instance import SymmetryInstance


class SymmetryCubic(AssetGeometryGroup):
    """
    Build a tetrahedral (T, 12 copies), octahedral (O, 24) or icosahedral (I, 60) assembly from one copy

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to replicate over the group
    group : InputMenu
        Which cubic point group to build
    orientation : InputRotation
        Rotation of the group's axes. At zero the two-folds lie along X, Y and Z and a three-fold along (1, 1, 1); for icosahedral a five-fold lies along (0, 1, 1.618)
    centre : InputVector
        Point all the symmetry axes pass through
    animate : InputFloat
        0 places a copy back on the original, 1 builds the full symmetry. Evaluated on each copy, so a field such as Stagger Value moves the copies one after another

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to replicate over the group
    i.group : MenuSocket
        Which cubic point group to build
    i.orientation : RotationSocket
        Rotation of the group's axes. At zero the two-folds lie along X, Y and Z and a three-fold along (1, 1, 1); for icosahedral a five-fold lies along (0, 1, 1.618)
    i.centre : VectorSocket
        Point all the symmetry axes pass through
    i.animate : FloatSocket
        0 places a copy back on the original, 1 builds the full symmetry. Evaluated on each copy, so a field such as Stagger Value moves the copies one after another

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    """

    _name = "Symmetry Cubic"
    _asset_name = "Symmetry Cubic"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Build a tetrahedral (T, 12 copies), octahedral (O, 24) or icosahedral (I, 60) assembly from one copy"
    }

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to replicate over the group"""
        group: MenuSocket
        """Which cubic point group to build"""
        orientation: RotationSocket
        """Rotation of the group's axes. At zero the two-folds lie along X, Y and Z and a three-fold along (1, 1, 1); for icosahedral a five-fold lies along (0, 1, 1.618)"""
        centre: VectorSocket
        """Point all the symmetry axes pass through"""
        animate: FloatSocket
        """0 places a copy back on the original, 1 builds the full symmetry. Evaluated on each copy, so a field such as Stagger Value moves the copies one after another"""

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
        group: InputMenu = "",
        orientation: InputRotation = None,
        centre: InputVector = None,
        animate: InputFloat = 1.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Group": group,
                "Orientation": orientation,
                "Centre": centre,
                "Animate": animate,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to replicate over the group"
        )
        group = tree.inputs.menu(
            "Group", "", description="Which cubic point group to build"
        )
        orientation = tree.inputs.rotation(
            "Orientation",
            (0.0, 0.0, 0.0),
            description="Rotation of the group's axes. At zero the two-folds lie along X, Y and Z and a three-fold along (1, 1, 1); for icosahedral a five-fold lies along (0, 1, 1.618)",
        )
        centre = tree.inputs.vector(
            "Centre",
            (0.0, 0.0, 0.0),
            description="Point all the symmetry axes pass through",
            subtype="XYZ",
        )
        animate = tree.inputs.float(
            "Animate",
            1.0,
            description="0 places a copy back on the original, 1 builds the full symmetry. Evaluated on each copy, so a field such as Stagger Value moves the copies one after another",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        instances = tree.outputs.geometry("Instances")

        index = g.Index()
        menu_switch = g.MenuSwitch.integer(
            group,
            {
                "Tetrahedral": (0, "T: 12 copies, as in Dps and dodecin"),
                "Octahedral": (1, "O: 24 copies, as in ferritin"),
                "Icosahedral": (2, "I: 60 copies, as in small virus capsids"),
            },
        )
        with g.Frame("Cubic operators"):
            index_switch = g.IndexSwitch.vector(
                index.o.index % 4,
                ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            )
            index_switch_1 = g.IndexSwitch.vector(
                menu_switch.o.output,
                ((0.0, 0.0, 1.0), (0.0, 0.0, 1.0), (0.0, 1.0, 1.618034)),
            )
            math_1 = (
                index.o.index
                // 12
                * g.IndexSwitch.float(
                    menu_switch.o.output, (0.0, math.pi / 2, math.tau / 5)
                )
            )
            axis_angle_to_rotation = g.AxisAngleToRotation(
                axis=g.RotateVector(rotation=orientation, vector=(1.0, 1.0, 1.0)),
                angle=index.o.index // 4 % 3 * (math.tau / 3),
            )
            _string = g.String(
                string="Every cubic group is built from the tetrahedral group T, which is the four rotations of D2 (identity and the two-folds about X, Y and Z) followed by 0, 120 or 240 degrees about the three-fold (1, 1, 1). Copy k uses D2 element k mod 4 and three-fold power (k div 4) mod 3. Octahedral adds a quarter turn about Z and icosahedral a fifth turn about the five-fold (0, 1, phi), raised to the power k div 12. Each axis is first turned by Orientation."
            )
            axis_angle_to_rotation_1 = g.AxisAngleToRotation(
                axis=index_switch.o.output.rotate(orientation),
                angle=(index.o.index % 4 > 0).switch.float(true=math.pi),
            )
            rotate_rotation = axis_angle_to_rotation_1.o.rotation.rotate(
                axis_angle_to_rotation
            ).rotate(
                g.AxisAngleToRotation(
                    axis=index_switch_1.o.output.rotate(orientation), angle=math_1
                )
            )
        (
            SymmetryInstance(
                geometry=geometry,
                points=g.Points(
                    count=g.IndexSwitch.integer(menu_switch.o.output, (12, 24, 60))
                ),
                rotation=rotate_rotation,
                centre=centre,
                animate=animate,
            )
            >> instances
        )


ASSET = SymmetryCubic

ASSET_METADATA = {
    "description": "Build a tetrahedral (T, 12 copies), octahedral (O, 24) or icosahedral (I, 60) assembly from one copy",
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
