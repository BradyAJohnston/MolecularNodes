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
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    MatrixSocket,
    ObjectSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry, InputObject
from .best_view import BestView
from .decay import Decay
from .fallback_camera import FallbackCamera


class StickyTransform(CustomGeometryGroup):
    _name = "Sticky Transform"
    _color_tag = "CONVERTER"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        smooth = tree.inputs.boolean(
            "Smooth",
            False,
            description="Ease towards each new best view over the simulation instead of jumping. Only advances during playback",
        )
        target = tree.inputs.geometry(
            "Target", description="Atoms holding the points that should be seen"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="The atoms that should be seen, such as a ligand",
            hide_value=True,
        )
        occluders = tree.inputs.geometry(
            "Occluders",
            description="Geometry to rotate. Its meshes are what can block the view",
        )
        transform = tree.inputs.matrix("Transform")
        stickiness = tree.inputs.float(
            "Stickiness",
            0.0,
            description="Favour best views near the one shown: the share of the selected atoms a turn of 180° must gain to be taken. Works with or without Smooth, and only advances during playback",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        delta_time = tree.inputs.float("Delta Time", 0.0)
        factor = tree.inputs.float(
            "Factor",
            1.0,
            description="Blend from the original orientation to the best view",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        half_life = tree.inputs.float(
            "Half Life",
            0.25,
            description="Seconds for the remaining gap to the best view to halve",
            min_value=0.0,
            subtype="TIME_ABSOLUTE",
        )
        camera_transform = tree.inputs.matrix("Camera Transform")
        transform_1 = tree.outputs.matrix("Transform")

        invert_rotation = transform.rotation.invert()
        switch = smooth.switch.float(
            1.0, Decay(delta_time=delta_time, half_life=half_life)
        )
        switch_1 = (stickiness > 0.0).switch.vector(
            (0.0, 0.0, 1.0),
            g.Vector((0.0, 1.0, 0.0))
            .o.vector.rotate(camera_transform)
            .rotate(invert_rotation),
        )
        best_view = BestView(
            target=target,
            selection=selection,
            occluders=occluders,
            up=switch_1,
            stickiness=stickiness,
            previous_direction=g.Vector((0.0, 0.0, -1.0))
            .o.vector.rotate(camera_transform)
            .rotate(invert_rotation),
        )
        mix = factor.mix.rotation(
            g.Rotation((0.0, 0.0, 0.0)),
            best_view.o.rotation.invert().rotate(camera_transform),
        )
        mix_1 = switch.mix.vector(
            transform.translation, best_view.o.centre - best_view.o.centre.rotate(mix)
        )
        combine_transform = g.CombineTransform(
            translation=mix_1, rotation=switch.mix.rotation(transform.rotation, mix)
        )

        combine_transform >> transform_1


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
    camera : InputObject
        Object to orient towards
    factor : InputFloat
        Blend from the original orientation to the best view
    smoothing : InputBoolean
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
    i.camera : ObjectSocket
        Object to orient towards
    i.factor : FloatSocket
        Blend from the original orientation to the best view
    i.smoothing : BooleanSocket
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
        camera: ObjectSocket
        """Object to orient towards"""
        factor: FloatSocket
        """Blend from the original orientation to the best view"""
        smoothing: BooleanSocket
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
        camera: InputObject = None,
        factor: InputFloat = 1.0,
        smoothing: InputBoolean = False,
        half_life: InputFloat = 0.25,
        stickiness: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Atoms": atoms,
                "Selection": selection,
                "Camera": camera,
                "Factor": factor,
                "Smoothing": smoothing,
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
            smoothing = tree.inputs.boolean(
                "Smoothing",
                False,
                description="Ease towards each new best view over the simulation instead of jumping. Only advances during playback",
                is_panel_toggle=True,
            )
            half_life = tree.inputs.float(
                "Half Life",
                0.25,
                description="Seconds for the remaining gap to the best view to halve",
                min_value=0.0,
                subtype="TIME_ABSOLUTE",
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

        closure_zone = g.ClosureZone()
        transform_1 = closure_zone.inputs.matrix("Transform")
        delta_time = closure_zone.inputs.float("Delta Time")
        transform_2 = closure_zone.outputs.matrix("Transform")
        sticky_transform = StickyTransform(
            Smooth=smoothing,
            Target=atoms,
            Selection=selection,
            Occluders=geometry,
            Transform=transform_1.output,
            Stickiness=stickiness,
            **{"Delta Time": delta_time.output},
            Factor=factor,
            **{
                "Half Life": half_life,
                "Camera Transform": FallbackCamera(camera=camera).o.transform,
            },
        )
        sticky_transform >> transform_2.input
        evaluate_closure = g.EvaluateClosure(closure_zone.closure)
        evaluate_closure.inputs.matrix(name="Transform")
        evaluate_closure.inputs.float(1_000_000.0, "Delta Time")
        transform_3 = evaluate_closure.outputs.matrix("Transform")
        simulation_zone = g.SimulationZone()
        transform_4 = simulation_zone.items.matrix(transform_3.output, "Transform")
        evaluate_closure_1 = g.EvaluateClosure(closure_zone.closure)
        evaluate_closure_1.inputs.matrix(transform_4.current, "Transform")
        evaluate_closure_1.inputs.float(simulation_zone.delta_time, "Delta Time")
        transform_5 = evaluate_closure_1.outputs.matrix("Transform")
        transform_5.output >> transform_4.next
        _string = g.String(
            "Each simulation step moves a = lerp(a, b, 1 - 2^(-dt / Half Life)) towards the target, so the gap halves every Half Life seconds whatever the frame rate. Translation and rotation are mixed separately, the rotation by slerp. Without Smooth each step lands on the target. The first frame starts on the target: outside the simulation the step is given a Delta Time long enough for the gap to close."
        )
        _string_1 = g.String(
            "The camera looks down its -Z. Undoing the rotation applied so far gives the direction it currently sees the unrotated geometry from, which Stickiness favours. With Stickiness, its +Y undone the same way is the Up for the next view, so the view rolls as little as it can. World Up would flip the roll each time the held view crosses near Up."
        )
        switch = (~smoothing | (stickiness > 0.0)).switch.matrix(
            transform_4.result, transform_3.output
        )
        geometry >> g.TransformGeometry(transform=switch, mode="Matrix") >> geometry_1

        switch >> transform


ASSET = OrientToBestView

ASSET_METADATA = {
    "description": "Rotate the geometry about the selection so the camera sees it from its Best View",
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
