import bpy

from .ui import (
    menu_item_interface, 
    button_custom_color, 
    button_custom_selection
)

class MN_MT_Node_Color(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_COLOR'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Set Color', 'MN_color_set', 
                            "Sets a new color for the selected atoms")
        layout.separator()
        menu_item_interface(layout, 'Goodsell Colors', 'MN_color_goodsell', 
                            "Adjusts the given colors to copy the 'Goodsell Style'.\n \
                            Darkens the non-carbon atoms and keeps the carbon atoms \
                            the same color. Highlights differences without being too \
                            visually busy")
        layout.separator()
        # menu_item_interface(layout, 'Color by B Factor', 'MN_color_map_attribute')
        menu_item_interface(layout, 'Attribute Map', 'MN_color_attribute_map')
        menu_item_interface(layout, 'Attribute Random', 'MN_color_attribute_random')
        layout.separator()
        # menu_item_color_chains(layout, 'Color by Chains')
        button_custom_color(layout, 'Chain', 'chain_id', 'Chain', 'chain_id_unique', 'chain')
        button_custom_color(layout, 'Entity', 'entity_id', '', 'entity_names', 'entity')
        button_custom_color(layout, 'Ligand', 'res_name', '', 'ligands', 'ligand', starting_value = 100)
        layout.separator()
        
        menu_item_interface(layout, 'Backbone', 'MN_color_backbone')
        menu_item_interface(layout, 'Secondary Structure', 'MN_color_sec_struct', 
                            "Specify colors based on the secondary structure")
        menu_item_interface(layout, 'Element', 'MN_color_element', 
                            "Choose a color for each of the first 20 elements")
        menu_item_interface(layout, 'Atomic Number', 'MN_color_atomic_number',
                            "Creates a color based on atomic_number field")
        menu_item_interface(layout, 'Res Name Peptide', 'MN_color_res_name_peptide')
        menu_item_interface(layout, 'Res Name Nucleic', 'MN_color_res_name_nucleic')
        menu_item_interface(layout, 'Element Common', 'MN_color_common', 
                            "Choose a color for the most common elements in PDB \
                            structures")

class MN_MT_Node_Bonds(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_BONDS'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Find Bonds', 'MN_bonds_find', 
                            "Finds bonds between atoms based on distance.\n\
                            Based on the vdw_radii for each point, finds other points \
                            within a certain radius to create a bond to. Does not \
                            preserve the index for the points. Does not detect bond type")
        menu_item_interface(layout, 'Break Bonds', 'MN_bonds_break', 
                            "Will delete a bond between atoms that already exists \
                            based on a distance cutoff")
        

class MN_MT_Node_Style(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_SYLE'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Presets', 'MN_style_presets')
        layout.separator()
        menu_item_interface(layout, 'Spheres', 'MN_style_spheres', 
                            'Create a sphere for each atom in the structure.')
        menu_item_interface(layout, 'Cartoon', 'MN_style_cartoon', 
                            'Create a cartoon representation, highlighting secondary \
                            structure through arrows and ribbons.')
        menu_item_interface(layout, 'Ribbon', 'MN_style_ribbon', 
                            'Create a "ribbon" or "licorice" style for peptide and nucleic acids.')
        menu_item_interface(layout, 'Surface', 'MN_style_surface', 
                            "Create a surface representation of the atoms.")
        menu_item_interface(layout, 'Ball and Stick', 'MN_style_ball_and_stick', 
                            "A style node to create ball and stick representation. \
                            Icospheres are instanced on atoms and cylinders for bonds. \
                            Bonds can be detected if they are not present in the \
                            structure")
        menu_item_interface(layout, 'Sticks', 'MN_style_sticks', 
                            "Turn each bond into a cylinder mesh")
        layout.separator()
        layout.label(text = 'Utilities')



class MN_MT_Node_Select(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_SELECT'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Separate Atoms', 'MN_select_separate_atoms', 
                            "Separate atoms based on a selection field.\n" +
                            "Takes atoms and splits them into the selected atoms the \
                            inverted atoms, based on a selection field")
        menu_item_interface(layout, 'Separate Polymers', 'MN_select_separate_polymers', 
                            "Separate the Geometry into the different polymers.\n" + 
                            "Outputs for protein, nucleic & sugars")
        layout.separator()
        button_custom_selection(layout, 'Chain', 'chain_id', 'Chain ', 'chain_id_unique', 'chain')
        button_custom_selection(layout, 'Entity', 'entity_id', '', 'entity_names', 'entity')
        button_custom_selection(layout, 'Ligands', 'res_name', '', 'ligands', 'ligand', starting_value = 100)
        layout.separator()
        menu_item_interface(layout, 'Secondary Structure', 'MN_select_sec_struct')
        menu_item_interface(layout, 'Backbone', 'MN_select_backbone', 
                            "Select atoms it they are part of the side chains or backbone.")
        menu_item_interface(layout, 'Atomic Number', 'MN_select_atomic_number', 
                            "Create a selection if input value equal to the \
                            atomic_number field.")
        menu_item_interface(layout, 'Element', 'MN_select_element', 
                            "Create a selection of particular elements by name. Only \
                            first 20 elements supported")
        menu_item_interface(layout, 'Attribute', 'MN_select_attribute')
        menu_item_interface(layout, 'Bonded Atoms', 'MN_select_bonded', 
                            "Based on an initial selection, finds atoms which are \
                            within a certain number of bonds away")
        layout.separator()
        menu_item_interface(layout, 'Proximity', 'MN_select_proximity', 
                            "Select atoms within a certain proximity of some target atoms.")
        menu_item_interface(layout, 'Cube', 'MN_select_cube', 
                            "Create a selection using an Empty Cube", 
                            node_link = False)
        menu_item_interface(layout, 'Sphere', 'MN_select_sphere', 
                            "Create a selection using an Empty Sphere", 
                            node_link = False)
        layout.separator()
        layout.operator('mn.residues_selection_custom', 
                        text = 'Res ID', 
                        emboss = True, 
                        depress = True)                        
        menu_item_interface(layout, 'Res ID Single', 'MN_select_res_id_single', 
                            "Create a selection if res_id matches input field")
        menu_item_interface(layout, 'Res ID Range', 'MN_select_res_id_range', 
                            "Create a selection if the res_id is within the given \
                            thresholds")
        menu_item_interface(layout, 'Res Name Peptide', 'MN_select_res_name_peptide', 
                            "Create a selection of particular amino acids by name")
        menu_item_interface(layout, 'Res Name Nucleic', 'MN_select_res_name_nucleic', 
                            "Create a selection of particular nucleic acids by name")
        menu_item_interface(layout, 'Res Whole', 'MN_select_whole_res', 
                            "Expand the selection to every atom in a residue, if any \
                            of those atoms are in the initial selection")

class MN_MT_Node_Assembly(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_ASSEMBLY'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        layout.operator("mn.assembly_bio", 
                        text = "Biological Assembly", 
                        emboss = True, 
                        depress=True
                        )
        menu_item_interface(layout, 'Center Assembly', 'MN_assembly_center', 
                            "Center the structure on the world origin based on \
                            bounding box")

class MN_MT_Node_Membranes(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_MEMBRANES'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Setup Atomic Properties', 'MN_prop_setup')

class MN_MT_Node_DNA(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_DNA'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Double Helix', 'MN_dna_double_helix', 
                            "Create a DNA double helix from an input curve.\n" + 
                            "Takes an input curve and instances for the bases, returns \
                            instances of the bases in a double helix formation")
        menu_item_interface(layout, 'Bases', 'MN_dna_bases', 
                            "Provide the DNA bases as instances to be styled and \
                            passed onto the Double Helix node")
        layout.separator()
        menu_item_interface(layout, 'Style Atoms Cyeles', 'MN_dna_style_atoms_cycles', 
                            "Style the DNA bases with spheres only visible in Cycles")
        menu_item_interface(layout, 'Style Spheres EEVEE', 'MN_dna_style_spheres_eevee', 
                            "Style the DNA bases with spheres visible in Cycles and \
                            EEVEE")
        menu_item_interface(layout, 'Style Surface', 'MN_dna_style_surface', 
                            "Style the DNA bases with surface representation")
        menu_item_interface(layout, 'Style Ball and Stick', 
                            'MN_dna_style_ball_and_stick', 
                            "Style the DNA bases with ball and stick representation")

class MN_MT_Node_Animate(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_ANIMATE'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Animate Frames', 'MN_animate_frames', 
                            "Interpolate between frames of a trajectory." + 
                            "Given a collection of frames for a trajectory, this node \
                            interpolates between them from start to finish based on \
                            the Animate field taking a value from 0 to 1. The \
                            positions of the Atoms are then moved based on this field")
        # menu_item_interface(layout, 'Animate Field', 'MN_animate_field')
        menu_item_interface(layout, 'Animate Value', 'MN_animate_value', 
                            "Animates between given start and end values, based on \
                            the input start and end frame of the timeline. Clamped \
                            will limit the output to the 'To Min' and 'To Max', while \
                            unclamped will continue to interpolate past these values. \
                            'Smoother Step' will ease in and out of these values, with \
                            default being linear interpolation")
        layout.separator()
        menu_item_interface(layout, 'Res Wiggle', "MN_animate_res_wiggle", 
                            "Wiggles the side chains of amino acids based on b_factor, \
                            adding movement to a structure.")
        menu_item_interface(layout, 'Res to Curve', "MN_animate_res_to_curve", 
                            "Takes atoms and maps them along a curve, as a single \
                            long peptide chain.")
        layout.separator()
        menu_item_interface(layout, 'Noise Position', 'MN_animate_noise_position', 
                            "Generate 3D noise field based on the position attribute")
        menu_item_interface(layout, 'Noise Field', 'MN_animate_noise_field', 
                            "Generate a 3D noise field based on the given field")
        menu_item_interface(layout, 'Noise Repeat', 'MN_animate_noise_repeat', 
                            "Generate a 3D noise field that repeats, based on the \
                            given field")

class MN_MT_Node_Utilities(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_UTILITIES'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Curve Resample', 'MN_utils_curve_resample')
        menu_item_interface(layout, 'Determine Secondary Structure', 'MN_utils_dssp')
        menu_item_interface(layout, 'Cartoon Utilities', 'MN_utils_style_cartoon')
        menu_item_interface(layout, 'Spheres Cycles', 'MN_utils_style_atoms_cycles', 
                            'A sphere atom representation, visible ONLY in Cycles. \
                            Based on point-cloud rendering')
        menu_item_interface(layout, 'Spheres EEVEE', 'MN_utils_style_spheres_eevee', 
                            'A sphere atom representation, visible in EEVEE and \
                            Cycles. Based on mesh instancing which slows down viewport \
                            performance')

class MN_MT_Node_CellPack(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_CELLPACK"
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        menu_item_interface(layout, 'Pack Instances', "MN_pack_instances")

class MN_MT_Node_Density(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_DENSITY'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Style Surface', 'MN_density_style_surface')
        menu_item_interface(layout, 'Style Wire', 'MN_density_style_wire')
        menu_item_interface(layout, 'Sample Nearest Attribute', 'MN_density_sample_searest')

class MN_MT_Node(bpy.types.Menu):
    bl_idname = "MN_MT_NODE"
    bl_label = "Menu for Adding Nodes in GN Tree"

    @classmethod
    def poll(cls, context):
        return not (False)

    def draw(self, context):
        layout = self.layout.column_flow(columns=1)
        layout.operator_context = "INVOKE_DEFAULT"
        layout.menu('MN_MT_NODE_SYLE', 
                    text='Style', icon_value=77)
        layout.menu('MN_MT_NODE_SELECT', 
                    text='Selection', icon_value=256)
        layout.menu('MN_MT_NODE_COLOR', 
                    text='Color', icon = 'COLORSET_07_VEC')
        layout.menu('MN_MT_NODE_ANIMATE', 
                    text='Animation', icon_value=409)
        layout.menu('MN_MT_NODE_ASSEMBLY', 
                    text='Assemblies', icon = 'GROUP_VERTEX')
        layout.menu('MN_MT_NODE_CELLPACK', 
                    text = 'CellPack', icon = 'PARTICLE_POINT')
        layout.menu('MN_MT_NODE_DENSITY', icon = "VOLUME_DATA", 
                    text = "Density")
        layout.menu('MN_MT_NODE_DNA', 
                    text='DNA', icon='GP_SELECT_BETWEEN_STROKES')
        # layout.menu('MN_MT_NODE_BONDS', 
        #             text='Bonds', icon = 'FIXED_SIZE')
        layout.menu('MN_MT_NODE_UTILITIES', 
                    text='Utilities', icon_value=92)

def MN_add_node_menu(self, context):
    if ('GeometryNodeTree' == bpy.context.area.spaces[0].tree_type):
        layout = self.layout
        layout.menu('MN_MT_NODE', text='Molecular Nodes', icon_value=88)
