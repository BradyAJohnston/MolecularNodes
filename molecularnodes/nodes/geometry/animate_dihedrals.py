# Node-group asset "Animate Dihedrals" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    CollectionSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputCollection, InputFloat, InputGeometry
from ._shared.animate_collection_pick import AnimateCollectionPick
from ._shared.animate_fraction import AnimateFraction
from .dihedral_chi_angle import DihedralChiAngle
from .dihedral_phi import DihedralPhi
from .dihedral_psi import DihedralPsi
from .sample_position import SamplePosition
from .set_chi_angle import SetChiAngle
from .set_phi_psi_angle import SetPhiPsiAngle


class SampleMixAngle(CustomGeometryGroup):
    _name = "Sample Mix Angle"
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Sample a float value from two different geometries and mix from A to B",
        "node_tool_idname": "geometry.sample_mix_angle",
    }

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        a = tree.inputs.geometry("A", description="Geometry A to sample and mix from")
        b = tree.inputs.geometry("B", description="Geometry B to sample and mix to")
        angle = tree.inputs.float(
            "Angle",
            0.0,
            description="Float field to sample and mix",
            hide_value=True,
            subtype="ANGLE",
        )
        factor = tree.inputs.float(
            "Factor",
            0.5,
            description="Amount to mix from A to B",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        index = tree.inputs.integer(
            "Index",
            0,
            description="`Index` on the geometries to sample from",
            default_input="INDEX",
        )
        value = tree.outputs.float(
            "Value", description="The final mixed float", subtype="ANGLE"
        )

        sample_index = g.SampleIndex(geometry=a, value=angle, index=index)
        sample_index_1 = g.SampleIndex(geometry=b, value=angle, index=index)
        mix = g.Mix(
            factor_float=factor,
            a_float=sample_index,
            b_float=sample_index_1,
            clamp_factor=True,
        )
        vector = g.Vector(vector=(0.0, 0.0, 1.0))
        mix_1 = g.Mix(
            factor_float=factor,
            a_rotation=g.AxisAngleToRotation(axis=vector, angle=sample_index),
            b_rotation=g.AxisAngleToRotation(axis=vector, angle=sample_index_1),
            data_type="ROTATION",
            clamp_factor=True,
        )
        _rotation_to_axis_angle = g.RotationToAxisAngle(
            rotation=mix_1.o.result_rotation
        )

        mix >> value


class AnimateDihedrals(AssetGeometryGroup):
    """
    Animate Dihedrals

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this node to
    frames : InputCollection
        Collection which holds the frames of the trajectory
    smoother_step : InputBoolean
        Ease in and out of the individual frames if interpolating
    interpolate : InputBoolean
        Whether to interpolate between frames of a trajectory or snap
    frame : InputFloat
        Which frame to select from the collection. The fraction component of the float is how much to interpolate between the current and next frame

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.frames : CollectionSocket
        Collection which holds the frames of the trajectory
    i.smoother_step : BooleanSocket
        Ease in and out of the individual frames if interpolating
    i.interpolate : BooleanSocket
        Whether to interpolate between frames of a trajectory or snap
    i.frame : FloatSocket
        Which frame to select from the collection. The fraction component of the float is how much to interpolate between the current and next frame

    Outputs
    -------
    o.atoms : GeometrySocket
        Atomic geometry with new positions based on the trajectory
    """

    _name = "Animate Dihedrals"
    _asset_name = "Animate Dihedrals"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.animate_dihedrals"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        frames: CollectionSocket
        """Collection which holds the frames of the trajectory"""
        smoother_step: BooleanSocket
        """Ease in and out of the individual frames if interpolating"""
        interpolate: BooleanSocket
        """Whether to interpolate between frames of a trajectory or snap"""
        frame: FloatSocket
        """Which frame to select from the collection. The fraction component of the float is how much to interpolate between the current and next frame"""

    class _Outputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry with new positions based on the trajectory"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        frames: InputCollection = None,
        smoother_step: InputBoolean = False,
        interpolate: InputBoolean = True,
        frame: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Frames": frames,
                "Smoother Step": smoother_step,
                "Interpolate": interpolate,
                "Frame": frame,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        frames = tree.inputs.collection(
            "Frames",
            description="Collection which holds the frames of the trajectory",
            optional_label=True,
        )
        smoother_step = tree.inputs.boolean(
            "Smoother Step",
            False,
            description="Ease in and out of the individual frames if interpolating",
        )
        interpolate = tree.inputs.boolean(
            "Interpolate",
            True,
            description="Whether to interpolate between frames of a trajectory or snap",
        )
        frame = tree.inputs.float(
            "Frame",
            0.0,
            description="Which frame to select from the collection. The fraction component of the float is how much to interpolate between the current and next frame",
            min_value=0.0,
            max_value=10_000.0,
        )
        atoms_1 = tree.outputs.geometry(
            "Atoms",
            description="Atomic geometry with new positions based on the trajectory",
        )

        group = AnimateCollectionPick(collection=frames, item=frame)
        group_1 = AnimateFraction(
            interpolate=interpolate, smoother_step=smoother_step, float=frame
        )
        object_info = g.ObjectInfo(object=g.SelfObject())
        set_position = g.SetPosition(
            geometry=atoms, position=SamplePosition(geometry=group.o.current)
        )
        set_position_1 = g.SetPosition(
            geometry=atoms, position=SamplePosition(geometry=group.o.next)
        )
        group_2 = SampleMixAngle(
            A=set_position,
            B=set_position_1,
            Angle=DihedralChiAngle().o.angle,
            Factor=group_1,
        )
        group_3 = SampleMixAngle(
            A=set_position,
            B=set_position_1,
            Angle=DihedralPhi(menu="Read").o.phi,
            Factor=group_1,
        )
        group_4 = SampleMixAngle(
            A=set_position,
            B=set_position_1,
            Angle=DihedralPsi(method="Read").o.psi,
            Factor=group_1,
        )
        group_5 = SetChiAngle(
            geometry=atoms,
            selection=selection,
            x1=group_2,
            x2=group_2,
            x3=group_2,
            x4=group_2,
            x5=group_2,
        )
        group_6 = SetPhiPsiAngle(geometry=group_5, phi=group_3, psi=group_4)
        group_6.node.mute = True
        (
            group_6
            >> g.TransformGeometry(
                translation=object_info.o.location,
                rotation=object_info.o.rotation,
                scale=abs(object_info.o.scale),
            )
            >> atoms_1
        )


ASSET = AnimateDihedrals

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
