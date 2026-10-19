# Node-group asset "Orient to Best View" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    MatrixSocket,
    ObjectSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry, InputObject
from .best_view import BestView


class OrientToBestView(AssetGeometryGroup):
    """
    Rotate the geometry about the selection so the camera sees it from its Best View

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to rotate. Its meshes are what can block the view
    atoms : InputGeometry
        Atoms holding the points that should be seen
    selection : InputBoolean
        The atoms that should be seen, such as a ligand
    active_camera : InputBoolean
        Orient towards the scene's active camera rather than Camera
    camera : InputObject
        Object to orient towards
    factor : InputFloat
        Blend from the original orientation to the best view
    smooth : InputBoolean
        Ease towards each new best view over the simulation instead of jumping. Only advances during playback
    half_life : InputFloat
        Seconds for the remaining gap to the best view to halve
    stickiness : InputFloat
        Favour best views near the one shown: the share of the selected atoms a turn of 180° must gain to be taken. Works with or without Smooth, and only advances during playback

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to rotate. Its meshes are what can block the view
    i.atoms : GeometrySocket
        Atoms holding the points that should be seen
    i.selection : BooleanSocket
        The atoms that should be seen, such as a ligand
    i.active_camera : BooleanSocket
        Orient towards the scene's active camera rather than Camera
    i.camera : ObjectSocket
        Object to orient towards
    i.factor : FloatSocket
        Blend from the original orientation to the best view
    i.smooth : BooleanSocket
        Ease towards each new best view over the simulation instead of jumping. Only advances during playback
    i.half_life : FloatSocket
        Seconds for the remaining gap to the best view to halve
    i.stickiness : FloatSocket
        Favour best views near the one shown: the share of the selected atoms a turn of 180° must gain to be taken. Works with or without Smooth, and only advances during playback

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    o.transform : MatrixSocket
        Transform applied to the geometry
    """

    _name = "Orient to Best View"
    _asset_name = "Orient to Best View"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Rotate the geometry about the selection so the camera sees it from its Best View"
    }

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to rotate. Its meshes are what can block the view"""
        atoms: GeometrySocket
        """Atoms holding the points that should be seen"""
        selection: BooleanSocket
        """The atoms that should be seen, such as a ligand"""
        active_camera: BooleanSocket
        """Orient towards the scene's active camera rather than Camera"""
        camera: ObjectSocket
        """Object to orient towards"""
        factor: FloatSocket
        """Blend from the original orientation to the best view"""
        smooth: BooleanSocket
        """Ease towards each new best view over the simulation instead of jumping. Only advances during playback"""
        half_life: FloatSocket
        """Seconds for the remaining gap to the best view to halve"""
        stickiness: FloatSocket
        """Favour best views near the one shown: the share of the selected atoms a turn of 180° must gain to be taken. Works with or without Smooth, and only advances during playback"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        transform: MatrixSocket
        """Transform applied to the geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        active_camera: InputBoolean = True,
        camera: InputObject = None,
        factor: InputFloat = 1.0,
        smooth: InputBoolean = False,
        half_life: InputFloat = 0.25,
        stickiness: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Atoms": atoms,
                "Selection": selection,
                "Active Camera": active_camera,
                "Camera": camera,
                "Factor": factor,
                "Smooth": smooth,
                "Half Life": half_life,
                "Stickiness": stickiness,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry",
            description="Geometry to rotate. Its meshes are what can block the view",
        )
        atoms = tree.inputs.geometry(
            "Atoms", description="Atoms holding the points that should be seen"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="The atoms that should be seen, such as a ligand",
            hide_value=True,
        )
        active_camera = tree.inputs.boolean(
            "Active Camera",
            True,
            description="Orient towards the scene's active camera rather than Camera",
        )
        camera = tree.inputs.object(
            "Camera", description="Object to orient towards", optional_label=True
        )
        factor = tree.inputs.float(
            "Factor",
            1.0,
            description="Blend from the original orientation to the best view",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        with tree.inputs.panel("Smoothing", default_closed=True):
            smooth = tree.inputs.boolean(
                "Smooth",
                False,
                description="Ease towards each new best view over the simulation instead of jumping. Only advances during playback",
            )
            half_life = tree.inputs.float(
                "Half Life",
                0.25,
                description="Seconds for the remaining gap to the best view to halve",
                min_value=0.0,
            )
            stickiness = tree.inputs.float(
                "Stickiness",
                0.0,
                description="Favour best views near the one shown: the share of the selected atoms a turn of 180° must gain to be taken. Works with or without Smooth, and only advances during playback",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
        geometry_1 = tree.outputs.geometry("Geometry")
        transform = tree.outputs.matrix(
            "Transform", description="Transform applied to the geometry"
        )

        best_view = BestView(target=atoms, selection=selection, occluders=geometry)
        with g.Frame("Rotate so the best view faces the camera"):
            _string = g.String(
                "The camera's rotation times the inverse of the best view's rotation carries the best view onto the camera, rotating about the centre of the selection so it stays put."
            )
            object_info = g.ObjectInfo(
                object=active_camera.switch.object(camera, g.ActiveCamera()),
                transform_space="RELATIVE",
            )
            mix = factor.mix.rotation(
                g.Rotation((0.0, 0.0, 0.0)),
                best_view.o.rotation.invert().rotate(object_info.o.rotation),
            )
            combine_transform = g.CombineTransform(
                translation=best_view.o.centre - best_view.o.centre.rotate(mix),
                rotation=mix,
            )
        _string_1 = g.String(
            "Each simulation step moves a = lerp(a, b, 1 - 2^(-dt / Half Life)) towards the target, so the gap halves every Half Life seconds whatever the frame rate. Translation and rotation are mixed separately, the rotation by slerp. Without Smooth each step lands on the target. The first frame starts on the target."
        )
        _string_2 = g.String(
            "The camera looks down its -Z. Undoing the rotation applied so far gives the direction it currently sees the unrotated geometry from, which Stickiness favours. With Stickiness, its +Y undone the same way is the Up for the next view, so the view rolls as little as it can. World Up would flip the roll each time the held view crosses near Up."
        )
        simulation_zone = g.SimulationZone()
        translation = simulation_zone.items.vector(
            combine_transform.o.transform.translation, "Translation"
        )
        rotation = simulation_zone.items.rotation(
            combine_transform.o.transform.rotation, "Rotation"
        )
        invert_rotation = rotation.current.invert()
        rotate_vector = (
            g.Vector((0.0, 0.0, -1.0))
            .o.vector.rotate(object_info.o.rotation)
            .rotate(invert_rotation)
        )
        switch = smooth.switch.float(
            1.0,
            1.0 - 2.0 ** (simulation_zone.delta_time * -1.0 / half_life.max(0.000001)),
        )
        rotate_vector_1 = (
            g.Vector((0.0, 1.0, 0.0))
            .o.vector.rotate(object_info.o.rotation)
            .rotate(invert_rotation)
        )
        best_view_1 = BestView(
            target=atoms,
            selection=selection,
            occluders=geometry,
            up=(stickiness > 0.0).switch.vector((0.0, 0.0, 1.0), rotate_vector_1),
            stickiness=stickiness,
            previous_direction=rotate_vector,
        )
        mix_1 = factor.mix.rotation(
            g.Rotation((0.0, 0.0, 0.0)),
            best_view_1.o.rotation.invert().rotate(object_info.o.rotation),
        )
        mix_2 = switch.mix.vector(
            translation.current,
            best_view_1.o.centre - best_view_1.o.centre.rotate(mix_1),
        )
        mix_2 >> translation.next
        switch.mix.rotation(rotation.current, mix_1) >> rotation.next
        switch_1 = (smooth | (stickiness > 0.0)).switch.matrix(
            combine_transform,
            g.CombineTransform(
                translation=translation.result, rotation=rotation.result
            ),
        )
        geometry >> g.TransformGeometry(transform=switch_1, mode="Matrix") >> geometry_1

        switch_1 >> transform


ASSET = OrientToBestView

ASSET_METADATA = {
    "description": "Rotate the geometry about the selection so the camera sees it from its Best View",
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
