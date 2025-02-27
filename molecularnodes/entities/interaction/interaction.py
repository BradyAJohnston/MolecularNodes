import bpy
import databpy
from databpy import Domains, AttributeTypes
import numpy as np
from bpy.types import Operator
from bpy.props import StringProperty
from mathutils import Vector
import json

from ..base import EntityType, MolecularEntity


class Interaction(MolecularEntity):
    @property
    def vertices(self):
        return self._vertices

    @vertices.setter
    def vertices(self, value):
        self._vertices = np.array(value)

    @property
    def edges(self):
        return self._edges

    @edges.setter
    def edges(self, value):
        self._edges = np.array(value)

    def __init__(self):
        super().__init__()
        self._entity_type = EntityType.INTERACTION
        self.frame_dict = {}
        self.trajectory_name = ""
        self.vertices = np.array([])
        self.edges = np.array([])
        self.EDGE_COLORS = {
            0: [1.0, 0.0, 0.0, 1.0],  # Red
            1: [0.0, 1.0, 0.0, 1.0],  # Green
            2: [0.0, 0.0, 1.0, 1.0],  # Blue
            3: [1.0, 1.0, 0.0, 1.0],  # Yellow
            4: [1.0, 0.5, 0.0, 1.0],  # Orange
            5: [0.5, 0.0, 1.0, 1.0],  # Purple
            6: [0.0, 1.0, 1.0, 1.0],  # Cyan
            7: [1.0, 0.0, 0.5, 1.0],  # Magenta
            "default": [0.5, 0.5, 0.5, 1.0],  # Gray
        }

        # Define interaction types as a class constant
        self.INTERACTION_TYPES = ["CationPi",
                                  "Hydrogen_Bond",
                                  "PiStacking",
                                  "Salt_Bridge",
                                  "Metal_Bond",
                                  "Halogen_Bond",
                                  "Hydrophobic",
                                  "VdWContact"]

    def setup_geometry_nodes(self) -> None:
        # Add a Geometry Nodes modifier
        geo_modifier = self.object.modifiers.new(
            name="GeometryNodes", type="NODES")

        # Create a new node tree for Geometry Nodes
        node_tree = bpy.data.node_groups.new(
            name="CurveConversion", type="GeometryNodeTree"
        )
        geo_modifier.node_group = node_tree

        # Set up the inputs and outputs for the node group
        node_tree.interface.new_socket(
            name="Mesh", in_out="INPUT", socket_type="NodeSocketGeometry"
        )
        node_tree.interface.new_socket(
            name="Curve", in_out="OUTPUT", socket_type="NodeSocketGeometry"
        )

        # Add nodes to the node tree
        nodes = node_tree.nodes
        group_input = nodes.new(type="NodeGroupInput")
        group_output = nodes.new(type="NodeGroupOutput")
        mesh_to_curve = nodes.new(type="GeometryNodeMeshToCurve")
        curve_circle = nodes.new(type="GeometryNodeCurvePrimitiveCircle")
        curve_to_mesh = nodes.new(type="GeometryNodeCurveToMesh")
        set_material = nodes.new(type="GeometryNodeSetMaterial")

        # Set the default properties for the nodes
        curve_circle.inputs["Radius"].default_value = 0.0025
        curve_circle.inputs["Resolution"].default_value = 16
        set_material.inputs["Material"].default_value = bpy.data.materials["MN Default"]

        # Set node positions
        group_input.location = (0, 0)
        mesh_to_curve.location = (200, 0)
        curve_circle.location = (200, -200)
        curve_to_mesh.location = (400, 0)
        set_material.location = (600, 0)
        group_output.location = (800, 0)

        # Set up node connections
        node_tree.links.new(
            group_input.outputs["Mesh"], mesh_to_curve.inputs["Mesh"])
        node_tree.links.new(
            mesh_to_curve.outputs["Curve"], curve_to_mesh.inputs["Curve"]
        )
        node_tree.links.new(
            curve_circle.outputs["Curve"], curve_to_mesh.inputs["Profile Curve"]
        )
        node_tree.links.new(
            curve_to_mesh.outputs["Mesh"], set_material.inputs["Geometry"]
        )
        node_tree.links.new(
            set_material.outputs["Geometry"], group_output.inputs["Curve"]
        )

    def _collect_interaction_geometry(
        self, frame: int = 0
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        vertices = []
        edge_types = []
        frame_str = str(frame)

        # Get the trajectory object
        blender_object = bpy.data.objects.get(self.trajectory_name)
        if not blender_object:
            return np.array([]), np.array([]), np.array([])

        # Collect vertices for all interactions
        for interaction_type, frames in self.frame_dict.items():
            frame_data = frames.get(frame_str, set())
            if not frame_data:
                continue

            for a, b in frame_data:
                interaction_type_index = self.INTERACTION_TYPES.index(
                    interaction_type)

                if interaction_type == "PiStacking":
                    # Both a and b are rings (multiple atoms)
                    pos_a = self.calculate_centroid(blender_object, a)
                    pos_b = self.calculate_centroid(blender_object, b)

                elif interaction_type == "CationPi":
                    # One group may have a single atom
                    if len(a) > 1 and len(b) == 1:
                        pos_a = self.calculate_centroid(blender_object, a)
                        pos_b = blender_object.data.vertices[int(b[0])].co
                    elif len(b) > 1 and len(a) == 1:
                        pos_a = blender_object.data.vertices[int(a[0])].co
                        pos_b = self.calculate_centroid(blender_object, b)
                    else:
                        # Default case if both have multiple atoms (unlikely but safe)
                        pos_a = self.calculate_centroid(blender_object, a)
                        pos_b = self.calculate_centroid(blender_object, b)

                else:
                    # For other interactions, a and b are individual atom indices
                    a_idx, b_idx = int(a), int(b)

                    # Get positions from named attributes
                    pos_a = databpy.named_attribute(
                        blender_object, "position")[a_idx]
                    pos_b = databpy.named_attribute(
                        blender_object, "position")[b_idx]

                vertices.extend([pos_a, pos_b])
                edge_types.append(interaction_type_index)

        # Convert to numpy arrays
        vertices_array = np.array(vertices)
        if len(vertices_array) > 0:
            edges_array = np.array(
                [(i, i + 1) for i in range(0, len(vertices_array), 2)]
            )
        else:
            edges_array = np.array([])

        return vertices_array, edges_array, np.array(edge_types)

    def compute_colors(self, edge_types: np.ndarray) -> np.ndarray:
        # Convert edge types to RGB colors
        return np.array(
            [self.EDGE_COLORS.get(et, [1.0, 1.0, 1.0, 1.0])
             for et in edge_types],
            dtype=np.float32,
        )

    def create_object(self, style: str = "vdw") -> bpy.types.Object:
        vertices, edges, edge_types = self._collect_interaction_geometry(0)

        self.object = databpy.create_object(
            vertices=vertices,
            edges=edges,
            name=f"{self.trajectory_name}_interactions",
            collection=bpy.data.collections["MolecularNodes"],
        )

        self.store_named_attribute(
            self.compute_colors(edge_types),
            "Color",
            atype=AttributeTypes.FLOAT_COLOR,
            domain=Domains.EDGE,
        )
        self.store_named_attribute(
            edge_types,
            "interaction_type",
            atype=AttributeTypes.INT,
            domain=Domains.EDGE,
        )

        return self.object

    def set_frame(self, frame: int) -> None:
        vertices, edges, edge_types = self._collect_interaction_geometry(frame)

        self.new_from_pydata(vertices, edges)
        self.store_named_attribute(
            self.compute_colors(edge_types),
            "Color",
            atype=AttributeTypes.FLOAT_COLOR,
            domain=Domains.EDGE,
        )
        self.store_named_attribute(
            edge_types,
            "interaction_type",
            atype=AttributeTypes.INT,
            domain=Domains.EDGE,
        )

    @staticmethod
    def calculate_centroid(trajectory_object, indices):
        vertices = [
            trajectory_object.data.vertices[int(idx)].co for idx in indices]
        return sum(vertices, Vector((0, 0, 0))) / len(vertices)

    def setup_from_json(self, json_file, object_name):
        self.trajectory_name = object_name

        with open(json_file, "r") as file:
            hb_data = json.load(file)

        for i_type, frames in hb_data.items():
            self.frame_dict[i_type] = {}
            for frame, pairs in frames.items():
                self.frame_dict[i_type][frame] = pairs


class MN_OT_interaction_visualiser(Operator):
    bl_idname = "mn.interaction_visualiser"
    bl_label = "Add Interactions"
    bl_description = "Loads a .json file with per-frame interaction information as an annotation object"
    bl_options = {"REGISTER", "UNDO"}

    filepath: StringProperty(
        name="JSON File",
        description="Path to the JSON file containing the interaction data",
        default="",
        subtype="FILE_PATH",
    )  # type: ignore

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {"RUNNING_MODAL"}

    def execute(self, context):
        # Create and setup interaction
        traj_name = context.active_object.name
        interaction = Interaction()
        interaction.setup_from_json(self.filepath, traj_name)
        interaction.create_object()
        interaction.setup_geometry_nodes()

        self.report({"INFO"}, "Interaction visualization setup complete")
        return {"FINISHED"}


CLASSES = [
    MN_OT_interaction_visualiser,
]
