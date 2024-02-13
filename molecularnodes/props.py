import bpy


bpy.types.Scene.MN_import_centre = bpy.props.BoolProperty(
    name="Centre Structure",
    description="Move the imported Molecule on the World Origin",
    default=False
)
bpy.types.Scene.MN_import_del_solvent = bpy.props.BoolProperty(
    name="Remove Solvent",
    description="Delete the solvent from the structure on import",
    default=True
)
bpy.types.Scene.MN_import_panel_selection = bpy.props.IntProperty(
    name="MN_import_panel_selection",
    description="Import Panel Selection",
    subtype='NONE',
    default=0
)
bpy.types.Scene.MN_import_build_assembly = bpy.props.BoolProperty(
    name='Build Assembly',
    default=False
)
bpy.types.Scene.MN_import_node_setup = bpy.props.BoolProperty(
    name="Setup Nodes",
    default=True,
    description='Create and set up a Geometry Nodes tree on import'
)


class MolecularNodesObjectProperties(bpy.types.PropertyGroup):
    subframes: bpy.props.IntProperty(
        name="Subframes",
        description="Number of subframes to interpolate for MD trajectories",
        default=0
    )
    molecule_type: bpy.props.StringProperty(
        name="Molecular Type",
        description="How the file was imported, dictating how MN interacts with it",
        default=""
    )
    pdb_code: bpy.props.StringProperty(
        name="PDB",
        description="PDB code used to download this structure",
        maxlen=4,
        options={'HIDDEN'}
    )
    star_type: bpy.props.StringProperty(
        name="Star Type",
        description="The type of star file that was imported",
        default=""
    )
