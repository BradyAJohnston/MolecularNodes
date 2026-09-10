# Node-group asset 'Assembly Instance' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    ObjectSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputObject,
)
from .assembly_id import AssemblyID
from .chain_id import ChainID
from .split_to_centred_instances import SplitToCentredInstances
from .symmetry_id import SymmetryID


class AssemblyInstance(AssetGeometryGroup):
    """
    Assembly Instance

    Parameters
    ----------
    geometry : InputGeometry
        Geometry that is instanced on the points
    selection : InputBoolean
        Selection
    data_object : InputObject
        Data Object
    assembly_id : InputInteger
        Assembly ID
    realize_all : InputBoolean
        Realize all levels of nested instances for a top-level instances. Overrides the value of the Depth input
    position : InputFloat
        Amount of mixing between the A and B inputs
    rotation : InputFloat
        Amount of mixing between the A and B inputs

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry that is instanced on the points
    i.selection : BooleanSocket
        Selection
    i.data_object : ObjectSocket
        Data Object
    i.assembly_id : IntegerSocket
        Assembly ID
    i.realize_all : BooleanSocket
        Realize all levels of nested instances for a top-level instances. Overrides the value of the Depth input
    i.position : FloatSocket
        Amount of mixing between the A and B inputs
    i.rotation : FloatSocket
        Amount of mixing between the A and B inputs

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    o.chain_id : IntegerSocket
        chain_id
    """

    _name = "Assembly Instance"
    _asset_name = "Assembly Instance"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry that is instanced on the points"""
        selection: BooleanSocket
        """Selection"""
        data_object: ObjectSocket
        """Data Object"""
        assembly_id: IntegerSocket
        """Assembly ID"""
        realize_all: BooleanSocket
        """Realize all levels of nested instances for a top-level instances. Overrides the value of the Depth input"""
        position: FloatSocket
        """Amount of mixing between the A and B inputs"""
        rotation: FloatSocket
        """Amount of mixing between the A and B inputs"""

    class _Outputs(SocketAccessor):
        instances: GeometrySocket
        """Instances"""
        chain_id: IntegerSocket

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        selection: InputBoolean = True,
        data_object: InputObject = None,
        assembly_id: InputInteger = 1,
        realize_all: InputBoolean = False,
        position: InputFloat = 1.0,
        rotation: InputFloat = 1.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Selection": selection,
                "Data Object": data_object,
                "Assembly ID": assembly_id,
                "Realize All": realize_all,
                "Position": position,
                "Rotation": rotation,
            }
        )

    def _build_group(self, tree):
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry that is instanced on the points"
        )
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        data_object = tree.inputs.object("Data Object", optional_label=True)
        assembly_id = tree.inputs.integer("Assembly ID", 1, min_value=1)
        realize_all = tree.inputs.boolean(
            "Realize All",
            False,
            description="Realize all levels of nested instances for a top-level instances. Overrides the value of the Depth input",
        )
        with tree.inputs.panel("Scale Transform", default_closed=True):
            position = tree.inputs.float(
                "Position",
                1.0,
                description="Amount of mixing between the A and B inputs",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
            rotation = tree.inputs.float(
                "Rotation",
                1.0,
                description="Amount of mixing between the A and B inputs",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
        instances = tree.outputs.geometry("Instances")
        chain_id = tree.outputs.integer("chain_id")

        with g.Frame():
            group = ChainID()
            group_1 = SplitToCentredInstances(geometry=geometry, group_id=group)
        with g.Frame("Select required assemblies"):
            separate_geometry = g.SeparateGeometry.point(
                g.ObjectInfo(object=data_object).o.geometry,
                g.Compare.integer.equal(assembly_id, AssemblyID()).o.result & selection,
            )
        with g.Frame("Mix Transform values"):
            attribute = g.NamedAttribute.input_4x4_matrix("transform").o.attribute
            mix = g.Mix(
                factor_float=rotation,
                b_rotation=attribute.rotation,
                data_type="ROTATION",
                clamp_factor=True,
            )
            mix_1 = g.Mix(
                factor_float=position,
                b_vector=attribute.translation,
                data_type="VECTOR",
                clamp_factor=True,
            )
        capture = g.CaptureAttribute.point(geometry=separate_geometry.o.selection)
        position_1 = capture.items.vector("Position", mix_1.o.result_vector)
        rotation_1 = capture.items.rotation("Rotation", mix.o.result_rotation)
        sym_id = capture.items.integer("sym_id", SymmetryID())
        instance_on_points = (
            capture.o.geometry
            >> g.SetPosition(position=position_1.output)
            >> g.InstanceOnPoints(
                instance=group_1,
                instance_index=group,
                rotation=rotation_1.output,
                pick_instance=True,
            )
        )
        capture_1 = g.CaptureAttribute.instance(geometry=instance_on_points)
        new_chain_id = capture_1.items.integer("new_chain_id", g.Index())
        (
            capture_1.o.geometry
            >> g.RealizeInstances(realize_all=realize_all)
            >> g.StoreNamedAttribute.point.integer(name="sym_id", value=sym_id.output)
            >> instances
        )

        new_chain_id.output >> chain_id


ASSET = AssemblyInstance

ASSET_METADATA = {
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
