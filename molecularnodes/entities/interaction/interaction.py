import bpy
from bpy.types import Operator, PropertyGroup
from bpy.props import StringProperty, FloatProperty
from mathutils import Vector
import json

from ..entity import EntityType, MolecularEntity


class Interaction(MolecularEntity):
    def __init__(self):
        super().__init__()
        self._entity_type = EntityType.INTERACTION
        self.frame_dict = {}
        self.bond_objects = {}
        self.trajectory_name = ""
        self.bond_width = 0.0025
        self._frame = None

    def set_frame(self, frame: int) -> None:
        # self._frame = self.frame_mapper(frame)
        self.update_positions(frame)

    @staticmethod
    def create_bond_material(interaction_type):
        material_name = f"{interaction_type}_Material"
        mat = bpy.data.materials.get(material_name)
        if mat is None:
            mat = bpy.data.materials.new(name=material_name)
            mat.use_nodes = True

            nodes = mat.node_tree.nodes
            nodes.clear()

            # Material nodes
            texture_coord = nodes.new(type="ShaderNodeTexCoord")
            color_ramp = nodes.new(type="ShaderNodeValToRGB")
            wave_texture = nodes.new(type="ShaderNodeTexWave")
            principled_bsdf = nodes.new(type="ShaderNodeBsdfPrincipled")
            output_node = nodes.new(type="ShaderNodeOutputMaterial")

            # Node positioning
            texture_coord.location = (-400, 0)
            wave_texture.location = (-200, 0)
            color_ramp.location = (0, 0)
            principled_bsdf.location = (200, 0)
            output_node.location = (400, 0)

            links = mat.node_tree.links
            links.new(texture_coord.outputs["UV"],
                      wave_texture.inputs["Vector"])
            links.new(wave_texture.outputs["Fac"], color_ramp.inputs["Fac"])
            links.new(color_ramp.outputs["Color"],
                      principled_bsdf.inputs["Alpha"])
            links.new(
                principled_bsdf.outputs["BSDF"], output_node.inputs["Surface"])

            # Customize based on interaction type
            if interaction_type == "Hydrogen":
                principled_bsdf.inputs["Base Color"].default_value = (
                    0, 1, 0, 1)  # Green
                # Tighter wave pattern
                wave_texture.inputs[1].default_value = 1
            elif interaction_type == "Salt_Bridge":
                principled_bsdf.inputs["Base Color"].default_value = (
                    1, 0.5, 0, 1)  # Orange
                # Moderate wave pattern
                wave_texture.inputs[1].default_value = 1
            elif interaction_type == "PiStacking":
                principled_bsdf.inputs["Base Color"].default_value = (
                    0, 0, 1, 1)  # Blue
                # Very tight wave pattern
                wave_texture.inputs[1].default_value = 1
            elif interaction_type == "CationPi":
                principled_bsdf.inputs["Base Color"].default_value = (
                    0.5, 0, 0.5, 1)  # Purple
                # Slightly denser wave pattern
                wave_texture.inputs[1].default_value = 1

            # Wave texture settings
            wave_texture.wave_type = "BANDS"
            wave_texture.bands_direction = "X"
            wave_texture.wave_profile = "SIN"

            # Wave transparency settings
            color_ramp.color_ramp.elements[0].position = 0.0
            color_ramp.color_ramp.elements[0].color = (
                1, 1, 1, 1)  # Opaque part
            color_ramp.color_ramp.interpolation = "CONSTANT"
            color_ramp.color_ramp.elements[1].position = 0.5
            color_ramp.color_ramp.elements[1].color = (
                0, 0, 0, 0)  # Transparent part

            principled_bsdf.inputs["Alpha"].default_value = 1

            mat.blend_method = "CLIP"
            mat.shadow_method = "CLIP"

        return mat

    @staticmethod
    def calculate_centroid(trajectory_object, indices):
        """Calculate the centroid of a set of vertex indices."""
        vertices = [
            trajectory_object.data.vertices[int(idx)].co for idx in indices]
        return sum(vertices, Vector((0, 0, 0))) / len(vertices)

    def create_bond_objects(self, all_bonds, bond_material, interaction_type):
        # Creates trajectory-children collections for each interaction
        interactions_collection = f"{self.trajectory_name}_Interactions"
        mn_collection = bpy.data.collections["MolecularNodes"]

        if interactions_collection not in mn_collection.children:
            bonds_collection = bpy.data.collections.new(
                interactions_collection)
            mn_collection.children.link(bonds_collection)

        interaction_collection = mn_collection.get(interaction_type)
        if not interaction_collection:
            interaction_collection = bpy.data.collections.new(interaction_type)
            bpy.data.collections[interactions_collection].children.link(
                interaction_collection)

        # Creates bond objects that will be given coordinates on update
        bond_objects = {}
        for couple in all_bonds:
            # bond_name = f"{interaction_type.lower()[:6]}_{str(couple[0])}_{str(couple[1])}"
            bond_name = f"{str(couple[0])}_{str(couple[1])}"
            curve_data = bpy.data.curves.new(name=bond_name, type="CURVE")
            curve_object = bpy.data.objects.new(bond_name, curve_data)
            interaction_collection.objects.link(curve_object)
            curve_data.dimensions = "3D"
            curve_data.bevel_depth = self.bond_width
            spline = curve_data.splines.new(type="POLY")
            spline.points.add(1)
            curve_object.data.materials.append(bond_material)
            bond_objects[couple] = curve_object

        return bond_objects

    def update_positions(self, current_frame):
        blender_object = bpy.data.objects.get(self.trajectory_name)
        if not blender_object or not blender_object.data:
            return

        frame = str(current_frame)
        for interaction_type, type_bonds in self.bond_objects.items():
            couples = self.frame_dict[interaction_type].get(frame, set())

            for bond, obj in type_bonds.items():
                if obj and obj.data:
                    visible = bond in couples
                    obj.hide_viewport = obj.hide_render = not visible

                    if visible:
                        spline = obj.data.splines[0]

                        if interaction_type == "PiStacking":
                            lig_centroid = self.calculate_centroid(
                                blender_object, bond[0])
                            prot_centroid = self.calculate_centroid(
                                blender_object, bond[1])
                            spline.points[0].co = (*lig_centroid, 1)
                            spline.points[1].co = (*prot_centroid, 1)
                        else:
                            spline.points[0].co = (
                                *blender_object.data.vertices[int(bond[0])].co, 1)
                            spline.points[1].co = (
                                *blender_object.data.vertices[int(bond[1])].co, 1)

    def setup_from_json(self, json_file, object_name):
        self.trajectory_name = object_name

        with open(json_file, "r") as file:
            hb_data = json.load(file)

        # Setup frame_dict
        for interaction_type, frames in hb_data.items():
            if interaction_type not in self.frame_dict:
                self.frame_dict[interaction_type] = {}

            for frame_key, interactions in frames.items():
                if frame_key not in self.frame_dict[interaction_type]:
                    self.frame_dict[interaction_type][frame_key] = set()

                for interaction in interactions:
                    if interaction_type == "PiStacking":
                        self.frame_dict[interaction_type][frame_key].add(
                            (tuple(interaction["Ligand"]),
                             tuple(interaction["Protein"]))
                        )
                    else:
                        self.frame_dict[interaction_type][frame_key].add(
                            (float(interaction["Ligand"]),
                             float(interaction["Protein"]))
                        )

        # Create bond objects
        for interaction_type, frames in self.frame_dict.items():
            all_bonds = set().union(
                *self.frame_dict[interaction_type].values())
            bond_material = self.create_bond_material(interaction_type)
            self.bond_objects[interaction_type] = self.create_bond_objects(
                all_bonds, bond_material, interaction_type
            )

    def update_bevel_depth(self, new_depth):
        self.bond_width = new_depth
        for interaction_type, bond_objects in self.bond_objects.items():
            for bond, obj in bond_objects.items():
                if obj and obj.type == 'CURVE':
                    obj.data.bevel_depth = new_depth*0.001


def update_bevel_depth(self, context):
    new_depth = self.bond_width  # Get the updated value
    scene = context.scene

    # Loop through the entities and update their bevel depth
    for obj in scene.MNSession.entities.values():
        if obj._entity_type.value == "interaction":
            obj.update_bevel_depth(new_depth)


class OBJECT_OT_interaction_visualiser(Operator):
    bl_idname = "object.interaction_visualiser"
    bl_label = "Visualize Interactions"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        scene = context.scene
        interaction_props = scene.interaction_visualiser
        json_file = interaction_props.json_file

        # Create and setup interaction
        interaction = Interaction()

        # Store the active object's name
        blender_object = context.active_object
        if blender_object:
            object_name = blender_object.name
            interaction.setup_from_json(json_file, object_name)

            self.report({"INFO"}, "Interaction visualization setup complete")
            return {"FINISHED"}
        else:
            self.report({"ERROR"}, "No active object selected")
            return {"CANCELLED"}


class InteractionVisualiserProperties(PropertyGroup):
    json_file: StringProperty(
        name="JSON File",
        description="Path to the JSON file containing data",
        default="",
        maxlen=1024,
        subtype="FILE_PATH",
    )
    object_name: StringProperty(
        name="Blender Object Name",
        description="Name of the Blender object representing the molecule",
        default="",
    )
    bond_width: FloatProperty(
        name="Bond Width",
        description="Update the bevel depth of the interaction objects",
        default=2.5,
        update=update_bevel_depth,
    )


CLASSES = [
    InteractionVisualiserProperties,
    OBJECT_OT_interaction_visualiser,
]
