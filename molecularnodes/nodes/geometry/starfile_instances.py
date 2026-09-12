# Node-group asset "Starfile Instances" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MaterialSocket,
    MenuSocket,
    ObjectSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputMaterial,
    InputMenu,
    InputObject,
)
from .fallback_geometry import FallbackGeometry
from .fallback_matrix import FallbackMatrix
from .fallback_rotation import FallbackRotation
from .primitive_gimbal import PrimitiveGimbal
from .rotation_cistem import RotationCisTEM
from .rotation_relion import RotationRELION


class StarfileInstances(AssetGeometryGroup):
    """
    Starfile Instances

    Parameters
    ----------
    points : InputGeometry
        Points
    selection : InputBoolean
        Becomes the output value if it is chosen by the menu input
    image : InputInteger
        The ID of the image that should be shown
    menu : InputMenu | Literal["Object", "Geometry"]
        Menu
    geometry : InputGeometry
        Becomes the output value if it is chosen by the menu input
    object : InputObject
        The object that should be placed at each instance
    pixel_scale : InputFloat
        The data in .ndjson files are in pixel coordinates so require scaling to fit overall data coordinates.
    instance_scale : InputFloat
        Scale of the instances
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.points : GeometrySocket
        Points
    i.selection : BooleanSocket
        Becomes the output value if it is chosen by the menu input
    i.image : IntegerSocket
        The ID of the image that should be shown
    i.menu : MenuSocket
        Menu
    i.geometry : GeometrySocket
        Becomes the output value if it is chosen by the menu input
    i.object : ObjectSocket
        The object that should be placed at each instance
    i.pixel_scale : FloatSocket
        The data in .ndjson files are in pixel coordinates so require scaling to fit overall data coordinates.
    i.instance_scale : FloatSocket
        Scale of the instances
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    """

    _name = "Starfile Instances"
    _asset_name = "Starfile Instances"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.starfile_instances"}

    class _Inputs(SocketAccessor):
        points: GeometrySocket
        """Points"""
        selection: BooleanSocket
        """Becomes the output value if it is chosen by the menu input"""
        image: IntegerSocket
        """The ID of the image that should be shown"""
        menu: MenuSocket
        """Menu"""
        geometry: GeometrySocket
        """Becomes the output value if it is chosen by the menu input"""
        object: ObjectSocket
        """The object that should be placed at each instance"""
        pixel_scale: FloatSocket
        """The data in .ndjson files are in pixel coordinates so require scaling to fit overall data coordinates."""
        instance_scale: FloatSocket
        """Scale of the instances"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""

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
        points: InputGeometry = None,
        selection: InputBoolean = True,
        image: InputInteger = 0,
        menu: InputMenu | Literal["Object", "Geometry"] = "Object",
        geometry: InputGeometry = None,
        object: InputObject = None,
        pixel_scale: InputFloat = 1.0,
        instance_scale: InputFloat = 1.0,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Points": points,
                "Selection": selection,
                "Image": image,
                "Menu": menu,
                "Geometry": geometry,
                "Object": object,
                "Pixel Scale": pixel_scale,
                "Instance Scale": instance_scale,
                "Material": material,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        points = tree.inputs.geometry("Points")
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Becomes the output value if it is chosen by the menu input",
            hide_value=True,
            structure_type="FIELD",
        )
        image = tree.inputs.integer(
            "Image",
            0,
            description="The ID of the image that should be shown",
            min_value=0,
            structure_type="SINGLE",
            force_non_field=True,
        )
        menu = tree.inputs.menu("Menu", expanded=True, optional_label=True)
        geometry = tree.inputs.geometry(
            "Geometry",
            description="Becomes the output value if it is chosen by the menu input",
        )
        object = tree.inputs.object(
            "Object",
            description="The object that should be placed at each instance",
            optional_label=True,
        )
        pixel_scale = tree.inputs.float(
            "Pixel Scale",
            1.0,
            description="The data in .ndjson files are in pixel coordinates so require scaling to fit overall data coordinates.",
            min_value=-10_000.0,
            max_value=10_000.0,
            structure_type="FIELD",
        )
        instance_scale = tree.inputs.float(
            "Instance Scale", 1.0, description="Scale of the instances"
        )
        material = tree.inputs.material(
            "Material",
            description="Material to apply to the resulting geometry",
            optional_label=True,
        )
        instances = tree.outputs.geometry("Instances")

        with g.Frame("Selection"):
            boolean_math = (
                g.Compare.integer.equal(
                    image, g.NamedAttribute.integer("image_id").o.attribute
                ).o.result
                & selection
            )
        with g.Frame("Instance"):
            menu_switch = g.MenuSwitch.geometry(
                menu,
                {
                    "Object": g.ObjectInfo(object=object).o.geometry,
                    "Geometry": geometry,
                },
            )
            transform_geometry = PrimitiveGimbal(
                vertices=3, material=material
            ) >> g.TransformGeometry(scale=g.Value(2.0))
            group = FallbackGeometry(geometry=menu_switch, fallback=transform_geometry)
        with g.Frame("Scale Pixel Coordinates"):
            set_position = points >> g.SetPosition(
                position=g.Position().o.position * pixel_scale
            )
        with g.Frame("Instances"):
            with g.Frame("Rotation"):
                group_1 = RotationCisTEM()
                group_2 = FallbackRotation(
                    name="rotation",
                    fallback=group_1.o.is_valid.switch.rotation(
                        RotationRELION().o.rotation, group_1.o.rotation
                    ),
                )
                group_3 = FallbackMatrix(name="transform", fallback=group_2)
            (
                set_position
                >> g.InstanceOnPoints(
                    selection=boolean_math,
                    instance=group,
                    rotation=group_3,
                    scale=instance_scale,
                )
                >> instances
            )

        menu.default_value = "Object"


ASSET = StarfileInstances

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
