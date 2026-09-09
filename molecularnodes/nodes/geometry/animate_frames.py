# Node-group asset 'Animate Frames' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    CollectionSocket,
    FloatSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputCollection, InputFloat, InputGeometry
from ._shared.animate_collection_pick import AnimateCollectionPick
from ._shared.animate_fraction import AnimateFraction
from .sample_mix_float import SampleMixFloat
from .sample_mix_vector import SampleMixVector
from .sample_position import SamplePosition
from .set_ures_id import SetUResID


class AnimateFrames(AssetGeometryGroup):
    """
    Animate Frames

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
    o.all_frames : GeometrySocket
        Each frame from the collection of coordinates is output in the a single mesh. They now contain an additional attribute `frame_id` to specify which structure they are from
    """

    _name = "Animate Frames"
    _asset_name = "Animate Frames"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.animate_frames"}

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
        all_frames: GeometrySocket
        """Each frame from the collection of coordinates is output in the a single mesh. They now contain an additional attribute `frame_id` to specify which structure they are from"""

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

    def _build_group(self, tree):
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
        all_frames = tree.outputs.geometry(
            "All Frames",
            description="Each frame from the collection of coordinates is output in the a single mesh. They now contain an additional attribute `frame_id` to specify which structure they are from",
        )

        with g.Frame("Interpolation position based on frames from collection"):
            group = AnimateFraction(
                interpolate=interpolate, smoother_step=smoother_step, float=frame
            )
            group_1 = AnimateCollectionPick(collection=frames, item=frame)
            group_2 = SampleMixFloat(
                a=group_1.o.current,
                b=group_1.o.next,
                value=g.NamedAttribute.float("b_factor").o.attribute,
                factor=group,
            )
            store_named_attribute = (
                atoms
                >> g.SetPosition(
                    selection=selection,
                    position=SampleMixVector(
                        a=group_1.o.current, b=group_1.o.next, factor=group
                    ),
                )
                >> g.StoreNamedAttribute.point.float(name="b_factor", value=group_2)
            )
        with g.Frame("Have to copy original structure as frames only contain position"):
            collection_info = g.CollectionInfo(
                collection=frames, separate_children=True
            )
            duplicate_elements = g.GeometryToInstance(
                atoms
            ) >> g.DuplicateElements.instance(
                amount=g.DomainSize(
                    geometry=collection_info, component="INSTANCES"
                ).o.instance_count
            )
            sample_index = g.SampleIndex(
                geometry=g.RealizeInstances(
                    geometry=collection_info, realize_to_point_domain=True
                ),
                value=g.Position(),
                index=g.Index(),
                data_type="FLOAT_VECTOR",
            )
            set_position = (
                duplicate_elements
                >> g.RealizeInstances(realize_to_point_domain=True)
                >> g.StoreNamedAttribute.point.integer(
                    name="frame_id", value=duplicate_elements.o.duplicate_index
                )
                >> g.SetPosition(position=sample_index)
            )
            group_3 = SetUResID(geometry=set_position)
        _group_4 = SamplePosition()
        with g.Frame("Set Transform to align with structure"):
            object_info = g.ObjectInfo(object=g.SelfObject())
            (
                group_3
                >> g.TransformGeometry(transform=object_info.o.transform, mode="Matrix")
                >> all_frames
            )
            (
                store_named_attribute
                >> g.TransformGeometry(transform=object_info.o.transform, mode="Matrix")
                >> atoms_1
            )


ASSET = AnimateFrames

ASSET_METADATA = {
    "catalog_id": "85730213-4c2e-469f-b333-52ac53adf274",
}
