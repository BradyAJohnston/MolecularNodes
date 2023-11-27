import bpy

from .ui import (
    menu_item_interface, 
    button_custom_color, 
    button_custom_selection
)
menu_items_color = [
    {
        'label': 'Set Color',
        'name': 'MN_color_set',
        'description': "Sets a new color for the selected atoms"
    },
    "break",
    {
        'label': 'custom',
        'function': button_custom_color, 
        'values': [
            {'label': 'Chain',  'field': 'chain_id', 'prefix': 'Chain', 'property_id': 'chain_id_unique'},
            {'label': 'Entity', 'field': 'entity_id', 'prefix': '',     'property_id': 'entity_names'},
            {'label': 'Ligand', 'field': 'res_name', 'prefix': '',      'property_id': 'ligands'}
        ]
    },
    "break",
    {
        'label': 'Goodsell Colors',
        'name': 'MN_color_goodsell',
        'description': "Adjusts the given colors to copy the 'Goodsell Style'.\n \
                        Darkens the non-carbon atoms and keeps the carbon atoms \
                        the same color. Highlights differences without being too \
                        visually busy"
    },
    {
        'label': 'Attribute Map',
        'name': 'MN_color_attribute_map',
        'description': ""
    },
    {
        'label': 'Attribute Random',
        'name': 'MN_color_attribute_random',
        'description': ""
    },
    {
        'label': 'Backbone',
        'name': 'MN_color_backbone',
        'description': ""
    },
    {
        'label': 'Secondary Structure',
        'name': 'MN_color_sec_struct',
        'description': "Specify colors based on the secondary structure"
    },
    "break",
    {
        'label': 'Element',
        'name': 'MN_color_element',
        'description': "Choose a color for each of the first 20 elements"
    },
    {
        'label': 'Atomic Number',
        'name': 'MN_color_atomic_number',
        'description': "Creates a color based on atomic_number field"
    },
    {
        'label': 'Res Name Peptide',
        'name': 'MN_color_res_name_peptide',
        'description': ""
    },
    {
        'label': 'Res Name Nucleic',
        'name': 'MN_color_res_name_nucleic',
        'description': ""
    },
    {
        'label': 'Element Common',
        'name': 'MN_color_common',
        'description': "Choose a color for the most common elements in PDB \
                        structures"
    }
]

menu_items_bonds = [
    {
        'label': 'Find Bonds',
        'name': 'MN_bonds_find',
        'description': "Finds bonds between atoms based on distance. Based on the vdw_radii for each point, finds other points within a certain radius to create a bond to. Does not preserve the index for the points. Does not detect bond type"
    },
    {
        'label': 'Break Bonds',
        'name': 'NB_bonds_break',
        'description': "Will delete a bond between atoms that already exists based on a distance cutoff"
    }
]

menu_items_style = [
    {
        'label': 'Presets',
        'name': 'MN_style_presets',
        'description': ''
    },
    {
        'label': 'Spheres',
        'name': 'MN_style_spheres',
        'description': 'Create a sphere for each atom in the structure.'
    },
    {
        'label': 'Cartoon',
        'name': 'MN_style_cartoon',
        'description': 'Create a cartoon representation, highlighting secondary structure through arrows and ribbons.'
    },
    {
        'label': 'Ribbon',
        'name': 'MN_style_ribbon',
        'description': 'Create a "ribbon" or "licorice" style for peptide and nucleic acids.'
    },
    {
        'label': 'Surface',
        'name': 'MN_style_surface',
        'description': 'Create a surface representation of the atoms.'
    },
    {
        'label': 'Ball and Stick',
        'name': 'MN_style_ball_and_stick',
        'description': 'A style node to create ball and stick representation. Icospheres are instanced on atoms and cylinders for bonds. Bonds can be detected if they are not present in the structure'
    },
    {
        'label': 'Sticks',
        'name': 'MN_style_sticks',
        'description': 'Turn each bond into a cylinder mesh'
    }
]

menu_item_selection = [
    {
        'label': 'Separate Atoms',
        'name': 'MN_select_separate_atoms',
        'description': "Separate atoms based on a selection field.\nTakes atoms and splits them into the selected atoms the inverted atoms, based on a selection field"
    },
    {
        'label': 'Separate Polymers',
        'name': 'MN_select_separate_polymers',
        'description': "Separate the Geometry into the different polymers.\nOutputs for protein, nucleic & sugars"
    },
    "break",
    {
        'label': 'custom',
        'function': button_custom_selection, 
        'values': [
            {'label': 'Chain',  'field': 'chain_id', 'prefix': 'Chain', 'property_id': 'chain_id_unique'},
            {'label': 'Entity', 'field': 'entity_id', 'prefix': '',     'property_id': 'entity_names'},
            {'label': 'Ligand', 'field': 'res_name', 'prefix': '',      'property_id': 'ligands'}
        ]
    },
    "break",
    {
        'label': 'Secondary Structure',
        'name': 'MN_select_sec_struct',
        'description': ''
    },
    {
        'label': 'Backbone',
        'name': 'MN_select_backbone',
        'description': "Select atoms it they are part of the side chains or backbone."
    },
    {
        'label': 'Atomic Number',
        'name': 'MN_select_atomic_number',
        'description': "Create a selection if input value equal to the atomic_number field."
    },
    {
        'label': 'Element',
        'name': 'MN_select_element',
        'description': "Create a selection of particular elements by name. Only first 20 elements supported"
    },
    {
        'label': 'Attribute',
        'name': 'MN_select_attribute',
        'description': ''
    },
    {
        'label': 'Bonded Atoms',
        'name': 'MN_select_bonded',
        'description': "Based on an initial selection, finds atoms which are within a certain number of bonds away"
    },
    "break",
    {
        'label': 'Proximity',
        'name': 'MN_select_proximity',
        'description': "Select atoms within a certain proximity of some target atoms."
    },
    {
        'label': 'Cube',
        'name': 'MN_select_cube',
        'description': "Create a selection using an Empty Cube"
    },
    {
        'label': 'Sphere',
        'name': 'MN_select_sphere',
        'description': "Create a selection using an Empty Sphere"
    },
    "break",
    {
        'label': 'Res ID',
        'name': 'mn.residues_selection_custom',
        'description': ''
    },
    {
        'label': 'Res ID Single',
        'name': 'MN_select_res_id_single',
        'description': "Create a selection if res_id matches input field"
    },
    {
        'label': 'Res ID Range',
        'name': 'MN_select_res_id_range',
        'description': "Create a selection if the res_id is within the given thresholds"
    },
    {
        'label': 'Res Name Peptide',
        'name': 'MN_select_res_name_peptide',
        'description': "Create a selection of particular amino acids by name"
    },
    {
        'label': 'Res Name Nucleic',
        'name': 'MN_select_res_name_nucleic',
        'description': "Create a selection of particular nucleic acids by name"
    },
    {
        'label': 'Res Whole',
        'name': 'MN_select_whole_res',
        'description': "Expand the selection to every atom in a residue, if any of those atoms are in the initial selection"
    }
]

menu_item_assembly = [
    {
        'label': 'Biological Assembly', 
        'name': 'mn.assembly_bio', 
        'description': ''
    }, 
    {
        'label': 'Center Assembly', 
        'name': 'MN_assembly_center', 
        'description': 'Center the structure on the world origin, based on the bounding box'
    }
]

menu_item_dna = [
    {
        'label': 'Double Helix',
        'name': 'MN_dna_double_helix',
        'description': "Create a DNA double helix from an input curve.\nTakes an input curve and instances for the bases, returns instances of the bases in a double helix formation"
    },
    {
        'label': 'Bases',
        'name': 'MN_dna_bases',
        'description': "Provide the DNA bases as instances to be styled and passed onto the Double Helix node"
    },
    "break",
    {
        'label': 'Style Atoms Cycles',
        'name': 'MN_dna_style_atoms_cycles',
        'description': "Style the DNA bases with spheres only visible in Cycles"
    },
    {
        'label': 'Style Spheres EEVEE',
        'name': 'MN_dna_style_spheres_eevee',
        'description': "Style the DNA bases with spheres visible in Cycles and EEVEE"
    },
    {
        'label': 'Style Surface',
        'name': 'MN_dna_style_surface',
        'description': "Style the DNA bases with surface representation"
    },
    {
        'label': 'Style Ball and Stick',
        'name': 'MN_dna_style_ball_and_stick',
        'description': "Style the DNA bases with ball and stick representation"
    }
]

menu_item_animate = [
    {
        'label': 'Animate Frames',
        'name': 'MN_animate_frames',
        'description': "Interpolate between frames of a trajectory. Given a collection of frames for a trajectory, this node interpolates between them from start to finish based on the Animate field taking a value from 0 to 1. The positions of the Atoms are then moved based on this field"
    },
    {
        'label': 'Animate Value',
        'name': 'MN_animate_value',
        'description': "Animates between given start and end values, based on the input start and end frame of the timeline. Clamped will limit the output to the 'To Min' and 'To Max', while unclamped will continue to interpolate past these values. 'Smoother Step' will ease in and out of these values, with default being linear interpolation"
    },
    "break",
    {
        'label': 'Res Wiggle',
        'name': 'MN_animate_res_wiggle',
        'description': "Wiggles the side chains of amino acids based on b_factor, adding movement to a structure."
    },
    {
        'label': 'Res to Curve',
        'name': 'MN_animate_res_to_curve',
        'description': "Takes atoms and maps them along a curve, as a single long peptide chain."
    },
    "break",
    {
        'label': 'Noise Position',
        'name': 'MN_animate_noise_position',
        'description': "Generate 3D noise field based on the position attribute"
    },
    {
        'label': 'Noise Field',
        'name': 'MN_animate_noise_field',
        'description': "Generate a 3D noise field based on the given field"
    },
    {
        'label': 'Noise Repeat',
        'name': 'MN_animate_noise_repeat',
        'description': "Generate a 3D noise field that repeats, based on the given field"
    }
]

menu_item_utils = [
    {'label': 'Curve Resample', 'name': 'MN_utils_curve_resample', 'description': ''},
    {'label': 'Determine Secondary Structure', 'name': 'MN_utils_dssp', 'description': ''},
    {'label': 'Cartoon Utilities', 'name': 'MN_utils_style_cartoon', 'description': ''},
    {'label': 'Spheres Cycles', 'name': 'MN_utils_style_atoms_cycles', 'description': 'A sphere atom representation, visible ONLY in Cycles. Based on point-cloud rendering'},
    {'label': 'Spheres EEVEE', 'name': 'MN_utils_style_spheres_eevee', 'description': 'A sphere atom representation, visible in EEVEE and Cycles. Based on mesh instancing which slows down viewport performance'}
]

menu_items_cellpack = [
    {
        'label' : 'Pack Instances',
        'name'  : 'MN_pack_instances', 
        'description': ''
    }
]

menu_items_density = [
        {
            'label': 'Style Surface',
            'name': 'MN_density_style_surface',
            'description': ''
        },
        {
            'label': 'Style Wire',
            'name': 'MN_density_style_wire',
            'description': ''
        },
        {
            'label': 'Sample Nearest Attribute',
            'name': 'MN_density_sample_nearest',
            'description': ''
        }
    ]


node_menu_menus = {
    'style': menu_items_style, 
    'selection': menu_item_selection, 
    'color': menu_items_color, 
    'animate': menu_item_animate, 
    'assemblies': menu_item_assembly, 
    'cellpack': menu_items_cellpack, 
    'density': menu_items_density, 
    'dna': menu_item_dna, 
    'utils': menu_item_utils
}

def build_menu(layout, items):
    for item in items:
        if item == "break":
            layout.separator()
        elif item['label'] == "custom":
            for button in item['values']:
                item['function'](layout, 
                    label = button['label'], 
                    field = button['field'], 
                    prefix = button['prefix'], 
                    property_id = button['property_id']
                    )
        elif item['name'].startswith("mn."):
            layout.operator(item['name'])
        else:
            menu_item_interface(layout, item['label'], item['name'], item['description'])

class MN_MT_Node_Color(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_COLOR'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        build_menu(layout, menu_items_color)


class MN_MT_Node_Bonds(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_BONDS'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        build_menu(layout, menu_items_bonds)


class MN_MT_Node_Style(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_STYLE'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        build_menu(layout, menu_items_style)

class MN_MT_Node_Select(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_SELECT'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        build_menu(layout, menu_item_selection)

class MN_MT_Node_Assembly(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_ASSEMBLY'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        build_menu(layout, menu_item_assembly)

class MN_MT_Node_DNA(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_DNA'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        build_menu(layout, menu_item_dna)

class MN_MT_Node_Animate(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_ANIMATE'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        build_menu(layout, menu_item_animate)

class MN_MT_Node_Utilities(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_UTILITIES'
    bl_label = ''
    
    def draw(self, context):
        build_menu(self.layout, menu_item_utils)

class MN_MT_Node_CellPack(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_CELLPACK"
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        build_menu(layout, menu_items_cellpack)

class MN_MT_Node_Density(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_DENSITY'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        build_menu(layout, menu_items_density)



class MN_MT_Node(bpy.types.Menu):
    bl_idname = "MN_MT_NODE"
    bl_label = "Menu for Adding Nodes in GN Tree"

    @classmethod
    def poll(cls, context):
        return True

    def draw(self, context):
        layout = self.layout.column_flow(columns=1)
        layout.operator_context = "INVOKE_DEFAULT"

        layout.menu('MN_MT_NODE_STYLE', text='Style', icon_value=77)
        layout.menu('MN_MT_NODE_SELECT', text='Selection', icon_value=256)
        layout.menu('MN_MT_NODE_COLOR', text='Color', icon='COLORSET_07_VEC')
        layout.menu('MN_MT_NODE_ANIMATE', text='Animation', icon_value=409)
        layout.menu('MN_MT_NODE_ASSEMBLY', text='Assemblies', icon='GROUP_VERTEX')
        layout.menu('MN_MT_NODE_CELLPACK', text='CellPack', icon='PARTICLE_POINT')
        layout.menu('MN_MT_NODE_DENSITY', text='Density', icon="VOLUME_DATA")
        layout.menu('MN_MT_NODE_DNA', text='DNA', icon='GP_SELECT_BETWEEN_STROKES')
        layout.menu('MN_MT_NODE_UTILITIES', text='Utilities', icon_value=92)

def MN_add_node_menu(self, context):
    if ('GeometryNodeTree' == bpy.context.area.spaces[0].tree_type):
        layout = self.layout
        layout.menu('MN_MT_NODE', text='Molecular Nodes', icon_value=88)
