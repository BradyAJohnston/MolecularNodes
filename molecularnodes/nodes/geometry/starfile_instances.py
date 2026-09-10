# Node-group asset 'Starfile Instances' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
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
    InputGeometry,
    InputInteger,
    InputMaterial,
    InputMenu,
    InputObject,
)
from .fallback_geometry import FallbackGeometry
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
    point_selection : InputMenu | Literal["Image", "Selection"]
        Select points to be instanced base don their `image_id` or a boolean selection field.
    selection : InputBoolean
        Becomes the output value if it is chosen by the menu input
    image : InputInteger
        The ID of the image that should be shown
    menu : InputMenu | Literal["Instance", "Simple"]
        Menu
    instance : InputObject
        The object that should be placed at each instance
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.points : GeometrySocket
        Points
    i.point_selection : MenuSocket
        Select points to be instanced base don their `image_id` or a boolean selection field.
    i.selection : BooleanSocket
        Becomes the output value if it is chosen by the menu input
    i.image : IntegerSocket
        The ID of the image that should be shown
    i.menu : MenuSocket
        Menu
    i.instance : ObjectSocket
        The object that should be placed at each instance
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    """

    _name = "Starfile Instances"
    _asset_name = "Starfile Instances"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.starfile_instances"}

    class _Inputs(SocketAccessor):
        points: GeometrySocket
        """Points"""
        point_selection: MenuSocket
        """Select points to be instanced base don their `image_id` or a boolean selection field."""
        selection: BooleanSocket
        """Becomes the output value if it is chosen by the menu input"""
        image: IntegerSocket
        """The ID of the image that should be shown"""
        menu: MenuSocket
        """Menu"""
        instance: ObjectSocket
        """The object that should be placed at each instance"""
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
        point_selection: InputMenu | Literal["Image", "Selection"] = "Image",
        selection: InputBoolean = True,
        image: InputInteger = 0,
        menu: InputMenu | Literal["Instance", "Simple"] = "Instance",
        instance: InputObject = None,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Points": points,
                "Point Selection": point_selection,
                "Selection": selection,
                "Image": image,
                "Menu": menu,
                "Instance": instance,
                "Material": material,
            }
        )

    def _build_group(self, tree):
        points = tree.inputs.geometry("Points")
        point_selection = tree.inputs.menu(
            "Point Selection",
            description="Select points to be instanced base don their `image_id` or a boolean selection field.",
            expanded=True,
            optional_label=True,
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Becomes the output value if it is chosen by the menu input",
            hide_value=True,
        )
        image = tree.inputs.integer(
            "Image",
            0,
            description="The ID of the image that should be shown",
            min_value=0,
        )
        menu = tree.inputs.menu("Menu", expanded=True, optional_label=True)
        instance = tree.inputs.object(
            "Instance",
            description="The object that should be placed at each instance",
            optional_label=True,
        )
        material = tree.inputs.material(
            "Material",
            description="Material to apply to the resulting geometry",
            optional_label=True,
        )
        instances = tree.outputs.geometry("Instances")

        with g.Frame("Instance"):
            transform_geometry = PrimitiveGimbal(
                vertices=3, material=material
            ) >> g.TransformGeometry(scale=g.Value(10.0))
            group = FallbackGeometry(
                geometry=g.ObjectInfo(object=instance).o.geometry,
                fallback=transform_geometry,
            )
            menu_switch = g.MenuSwitch.geometry(
                menu, {"Instance": group, "Simple": transform_geometry}
            )
        with g.Frame("Rotation"):
            group_1 = RotationCisTEM()
            group_2 = FallbackRotation(
                name="rotation",
                fallback=group_1.o.is_valid.switch.rotation(
                    RotationRELION().o.rotation, group_1.o.rotation
                ),
            )
        with g.Frame("Selection"):
            menu_switch_1 = g.MenuSwitch.boolean(
                point_selection,
                {
                    "Image": g.Compare.integer.equal(
                        image, g.NamedAttribute.integer("image_id").o.attribute
                    ),
                    "Selection": selection,
                },
            )
        (
            points
            >> g.InstanceOnPoints(
                selection=menu_switch_1.o.output, instance=menu_switch, rotation=group_2
            )
            >> instances
        )

        point_selection.default_value = "Image"
        menu.default_value = "Instance"


ASSET = StarfileInstances

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}

DATABLOCK_DEPENDENCIES = {
    "materials": ("MN Default.old",),
}
