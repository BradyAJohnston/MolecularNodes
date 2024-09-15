from typing import List, Union
from abc import ABCMeta
import bpy


class Item:
    def __init__(self) -> None:
        self.is_break = False
        self.is_custom = False
        self.backup: str = None

    @classmethod
    def menu(
        self, layout: bpy.types.UILayout, context: bpy.types.Context = None
    ) -> None:
        pass


class MenuItem(Item):
    def __init__(
        self,
        name: str,
        label: str = None,
        description: str = None,
        video_url: str = None,
        backup: str = None,
    ) -> None:
        super().__init__()
        self.name = name
        self.label = label
        if self.label is None:
            self.label = self.name
        self.description = description
        self.video_url = video_url
        self.backup = backup

    def short_description(self):
        return self.description.split("\n")[0].removesuffix(".")

    def menu(
        self, layout: bpy.types.UILayout, context: bpy.types.Context = None
    ) -> None:
        if self.label is None:
            self.label = self.name

        if self.name.startswith("mn."):
            layout.operator(self.name)
            return None

        op = layout.operator("mn.add_custom_node_group", text=self.label)
        op.node_label = self.label
        op.node_name = self.name
        op.node_description = self.description
        op.node_link = False


class CustomItem(Item):
    def __init__(
        self,
        label: str,
        field: str,
        dtype: str,
        name: str,
        prefix: str,
        property_id: str,
        description: str,
        video_url: str = None,
    ) -> None:
        super().__init__()
        self.label = label
        self.field = field
        self.dtype = dtype
        self.name = name
        self.prefix = prefix
        self.property_id = property_id
        self.description = description
        self.video_url = video_url
        self.is_custom = True

    def menu(
        self, layout: bpy.types.UILayout, context: bpy.types.Context = None
    ) -> None:
        row = layout.row()
        op = row.operator("mn.iswitch_custom", text=self.label)
        op.field = self.field
        op.dtype = self.dtype
        op.prefix = self.prefix
        op.node_property = self.property_id
        op.node_name = self.label

        if self.dtype == "RGBA":
            op.description = f"Choose custom colors for {self.label}"
        elif self.dtype == "BOOLEAN":
            op.description = f"Choose custom selections for {self.label}"
        else:
            raise ValueError(f"Data type currently not supported: {self.dtype}")
        row.enabled = bool(context.active_object.get(self.property_id))


class Break:
    def __init__(self, text: str = None) -> None:
        super().__init__()
        self.is_break = True
        self.text = text

    def menu(
        self, layout: bpy.types.UILayout, context: bpy.types.Context = None
    ) -> None:
        layout.separator()
        # optionally we can add a subtitle for the next section of nodes
        if self.text and self.text.strip() != "":
            layout.label(text=self.text)


class Submenu:
    def __init__(self, name, items, title: str = None, description: str = None) -> None:
        self.name: str = name
        self.items: List[Union[MenuItem, Break, CustomItem]] = items
        self.title = title
        self.description = description

    def node_names(self):
        return [item.name for item in self.items if not isinstance(item, Break)]

    def menu(self, layout: bpy.types.UILayout, context: bpy.types.Context):
        for item in self.items:
            item.menu(layout=layout, context=context)


class Menu:
    def __init__(self, submenus: List[Submenu]) -> None:
        self.submenus = submenus

    def get_submenu(self, name: str) -> Submenu:
        for sub in self.submenus:
            if sub.name == name:
                return sub


menu_items = Menu(
    submenus=[
        Submenu(
            name="style",
            title="Style",
            items=[
                MenuItem(
                    label="Spheres",
                    name="Style Spheres",
                    description="Style to apply the traditional space-filling atomic representation of atoms. Spheres are scaled based on the `vdw_radii` attribute. By default the _Point Cloud_ rendering system is used, which is only visible inside of Cycles.",
                    video_url="https://imgur.com/KjKkF2u",
                ),
                MenuItem(
                    label="Cartoon",
                    name="Style Cartoon",
                    description="Style to apply the traditional cartoon representation of protein structures. This style highlights alpha-helices and beta-sheets with arrows and cylinders.",
                    video_url="https://imgur.com/1xmdfxZ",
                ),
                MenuItem(
                    label="Ribbon",
                    name="Style Ribbon",
                    description="Style that creates a continuous solid ribbon or licorice tube through the backbones of peptides and nucleic acids.",
                    video_url="https://imgur.com/iMxEJaH",
                ),
                MenuItem(
                    label="Surface",
                    name="Style Surface",
                    description="Style that creates a surface representation based on the proximity of atoms to a probe that is moved through the entire structure.",
                    video_url="https://imgur.com/ER8pcYf",
                ),
                MenuItem(
                    label="Ball and Stick",
                    name="Style Ball and Stick",
                    description="Style that creates cylinders for bonds and spheres for atoms. The atoms can be either Eevee or Cycles compatible, with customisation to resolution and radius possible.",
                    video_url="https://imgur.com/kuWuOsw",
                ),
                MenuItem(
                    label="Sticks",
                    name="Style Sticks",
                    description="Style that creates a cylinder for each bond. Cylindrical caps to the cylinders are currently not supported. Best to use [`Style Ball and Stick`](#style-ball-and-stick).",
                    video_url="https://imgur.com/4ZK1AMo",
                ),
                Break(),
                MenuItem(
                    label="Preset 1",
                    name="Style Preset 1",
                    description="Quickly switch between several different pre-made preset styles. Best used when using MolecularNodes via scripts, ensuring all atoms are displayed using a combination of cartoons and atoms.",
                    video_url="https://imgur.com/gCQRWBk.mp4",
                ),
                MenuItem(
                    label="Preset 2",
                    name="Style Preset 2",
                    description="Quickly switch between several different pre-made preset styles. Best used when using MolecularNodes via scripts, ensuring all atoms are displayed using a combination of cartoons and atoms.",
                    video_url="https://imgur.com/gCQRWBk.mp4",
                ),
                MenuItem(
                    label="Preset 3",
                    name="Style Preset 3",
                    description="Quickly switch between several different pre-made preset styles. Best used when using MolecularNodes via scripts, ensuring all atoms are displayed using a combination of cartoons and atoms.",
                    video_url="https://imgur.com/gCQRWBk.mp4",
                ),
                MenuItem(
                    label="Preset 4",
                    name="Style Preset 4",
                    description="Quickly switch between several different pre-made preset styles. Best used when using MolecularNodes via scripts, ensuring all atoms are displayed using a combination of cartoons and atoms.",
                    video_url="https://imgur.com/gCQRWBk.mp4",
                ),
            ],
        ),
        Submenu(
            name="select",
            title="Select",
            items=[
                MenuItem(
                    name="Separate Atoms",
                    description="Select only the desired input atoms. The output is bits of geometry, which include the selection and include the inverse of the selected atoms. You can expand the selection to include an entire residue if a single atom in that residue is selected, by setting `Whole Residue` to `True`.",
                    video_url="https://imgur.com/VsCW0HY",
                ),
                MenuItem(
                    name="Separate Polymers",
                    description="Separate the input atomic geometry into it's different polymers or `Protein`, `Nucleic Acid` and `other`.",
                    video_url="https://imgur.com/ICQZxxz",
                ),
                Break(),
                CustomItem(
                    label="Chain",
                    field="chain_id",
                    dtype="BOOLEAN",
                    name="Select Chain_",
                    prefix="",
                    property_id="chain_ids",
                    description="Select single or multiple of the different chains. Creates a selection based on the `chain_id` attribute.",
                    video_url="https://imgur.com/P9ZVT2Z",
                ),
                CustomItem(
                    label="Entity",
                    field="entity_id",
                    name="Select Entity_",
                    dtype="BOOLEAN",
                    prefix="",
                    property_id="entity_ids",
                    description="Select single or multiple of the different entities. Creates a selection based on the `entity_id` attribute.",
                    video_url="https://imgur.com/fKQIfGZ",
                ),
                CustomItem(
                    label="Ligand",
                    field="res_name",
                    name="Select Ligand_",
                    dtype="BOOLEAN",
                    prefix="",
                    property_id="ligands",
                    description="Select single or multiple of the different ligands.",
                    video_url="https://imgur.com/s2seWIw",
                ),
                CustomItem(
                    label="Segment",
                    field="segid",
                    name="Select Segment_",
                    dtype="BOOLEAN",
                    prefix="",
                    property_id="segments",
                    description="",
                ),
                MenuItem(
                    label="Atomic Number",
                    name="Select Atomic Number",
                    description="Select single elements, by matching to the `atomic_number` field. Useful for selecting single elements, or combining to select elements higher than 20 on the periodic table.",
                    video_url="https://imgur.com/Bxn33YK",
                ),
                MenuItem(
                    label="Element",
                    name="Select Element",
                    description="Select individual elements, for the first 20 elements on the periodic table. For selections of higher elements, use [`MN_select_atomic_number`](#select-atomic-number). Creating a node which includes more elements becomes too large to be practical.",
                    video_url="https://imgur.com/nRQwamG",
                ),
                MenuItem(
                    label="Res Name",
                    name="Select Res Name",
                    description="Select protein or nucleic acids based on their residue name.",
                    video_url="https://imgur.com/kjzH9Rs",
                ),
                MenuItem(
                    label="Res ID Single",
                    name="Select Res ID",
                    description="Select a atoms based on their `res_id` number.",
                    video_url="https://imgur.com/BL6AOP4",
                ),
                MenuItem(
                    label="Res ID Range",
                    name="Select Res ID Range",
                    description="Select multiple residues by specifying a _Min_ and a _Max_, defining a range that includes or excludes based on the `res_id` number.",
                    video_url="https://imgur.com/NdoQcdE",
                ),
                MenuItem(
                    label="Res ID",
                    name="mn.residues_selection_custom",
                    backup="Select Res ID_",
                    description="Create a more complex selection for the `res_id` field, by specifying multiple ranges and potential single `res_id` numbers. This node is built uniquely each time, to the inputs will look different for each user.\nIn the example below, residues 10 & 15 are selected, as well as residues between and including 20-100.\nThe node was created by inputting `10, 15, 20-100` into the node creation field.",
                    video_url="https://imgur.com/OwAXsbG",
                ),
                Break(),
                MenuItem(
                    label="Attribute",
                    name="Select Attribute",
                    description="Select atoms that have true for the given attribute name.",
                ),
                MenuItem(
                    name="Is Peptide",
                    description="Select the atoms involved in a peptide chain.",
                ),
                MenuItem(
                    name="Is Nucleic",
                    description="Select the atoms involved in nucleic acid polymer.",
                ),
                MenuItem(
                    name="Is Lipid",
                    description="Select the atoms involved in lipid molecules.",
                ),
                MenuItem(
                    name="Is Solvent",
                    description="Select the atoms that are part of the solvent.",
                ),
                MenuItem(
                    name="Is Alpha Carbon",
                    description="Select the alpha carbons of a peptide.",
                ),
                MenuItem(
                    name="Is Backbone",
                    description="Select the backbone atoms of a peptide or nucleic acid polymer.",
                ),
                MenuItem(
                    name="Is Side Chain",
                    description="Select the side chain atoms of a peptide or nucleic acid polymer.",
                ),
                MenuItem(
                    name="Is Helix",
                    description="Select the atoms in a alpha-helix or similar.",
                ),
                MenuItem(
                    name="Is Sheet",
                    description="Select the atoms in a beta-sheet or similar.",
                ),
                MenuItem(
                    name="Is Loop",
                    description="Select the atoms not in a sheet or helix.",
                ),
                Break(),
                MenuItem(
                    label="Bonded",
                    name="Select Bonded",
                    description="Based on an initial selection, finds atoms which are within a certain number of bonds of this selection. Output can include or excluded the original selection.",
                    video_url="https://imgur.com/g8hgXup",
                ),
                MenuItem(
                    label="Res Whole",
                    name="Select Res Whole",
                    description="Expand the given selection to include a whole residue, if a single atom in that residue is selected. Useful for when a distance or proximity selection includes some of the residue and you wish to include all of the residue.",
                    video_url="https://imgur.com/JFzwE0i",
                ),
                Break(),
                MenuItem(
                    label="Cube",
                    name="Select Cube",
                    description="Create a selection that is inside the `Empty_Cube` object. When this node is first created, an _empty_ object called `Empty_Cube` should be created. You can always create additional empty objects through the add menu, to use a different object. The rotation and scale of the object will be taken into account for the selection.",
                    video_url="https://imgur.com/P4GZ7vq",
                ),
                MenuItem(
                    label="Sphere",
                    name="Select Sphere",
                    description="Create a selection that is within a spherical radius of an object, based on that object's scale. By default an _empty_ object called `Empty_Sphere` is created. You can use other objects or create a new empty to use. The origin point for the object will be used, which should be taken in to account when using molecules. Use [`MN_select_proximity`](#select-proximity) for selections which are within a certain distance of a selection of atoms instead of a single origin point.",
                    video_url="https://imgur.com/xdeTZR7",
                ),
                MenuItem(
                    label="Proximity",
                    name="Select Proximity",
                    description="Create a selection based on the proximity to the Target Atoms of the input. A sub-selection of the Target atoms can be used if the `Selection` input is used. You can expand the selection to include an entire residue if a single atom in that residue is selected, by setting `Whole Residue` to `True`.\nIn the example below, the `Style Atoms` is being applied to a selection, which is being calculated from the proximity of atoms to specific chains. As the cutoff for the selection is changed, it includes or excludes more atoms. The `Whole Residue` option also ensures that entire residues are shown.",
                    video_url="https://imgur.com/RI80CRY",
                ),
            ],
        ),
        Submenu(
            name="color",
            title="Color",
            items=[
                MenuItem(
                    label="Set Color",
                    name="Set Color",
                    description="The is the primary way to change the color of structures in Molecular Nodes. Colors for cartoon and ribbon are taken from the _alpha-carbons_ of the structures. Change the color of the input atoms, based on a selection and a color field. The color field can be as complex of a calculation as you wish. In the example below the color for the whole structure can be set, or the color can be based on a color for each chain, or the result of mapping a color to an attribute such as `b_factor`.",
                    video_url="https://imgur.com/667jf0O",
                ),
                Break(),
                CustomItem(
                    label="Chain",
                    field="chain_id",
                    dtype="RGBA",
                    name="Color Chain_",
                    prefix="",
                    property_id="chain_ids",
                    description="Choose the colors for individual chains in the structure. This node is generated for each particular molecule, so the inputs will look different based on the imported structure. For larger structures with many chains this node may become too large to be practical, in which case you might better use [`Color Entity ID`](#color-entity-id).",
                    video_url="https://imgur.com/9oM24vB",
                ),
                CustomItem(
                    label="Segment",
                    field="segid",
                    name="Color Segment_",
                    dtype="RGBA",
                    prefix="",
                    property_id="segments",
                    description="",
                ),
                CustomItem(
                    label="Entity",
                    field="entity_id",
                    name="Color Entity_",
                    dtype="RGBA",
                    prefix="",
                    property_id="entity_ids",
                    description="Choose the colors for individual entities in the structure. Multiple chains may be classified as the same entity, if they are copies of the same chain but in different conformations or positions and rotations. The nodes is generated for each individual structure, if `entity_id` is available.",
                    video_url="https://imgur.com/kEvj5Jk",
                ),
                CustomItem(
                    label="Ligand",
                    field="res_name",
                    name="Color Ligand_",
                    dtype="RGBA",
                    prefix="",
                    property_id="ligands",
                    description="Choose the colors for individual ligands in the structure.",
                    video_url="https://imgur.com/bQh8Fd9",
                ),
                MenuItem(
                    label="Element",
                    name="Color Element",
                    description="Choose a color for each of the first 20 elements on the periodic table. For higher atomic number elements use [`Color Atomic Number`](#color-atomic-number).",
                    video_url="https://imgur.com/iMGZKCx",
                ),
                MenuItem(
                    label="Atomic Number",
                    name="Color Atomic Number",
                    description="Choose a color for an individual element. Select the element based on `atomic_number`. Useful for higher atomic number elements which are less commonly found in structures.",
                    video_url="https://imgur.com/pAloaAF",
                ),
                MenuItem(
                    label="Res Name",
                    name="Color Res Name",
                    description="Choose a color for each of the 20 naturally occurring amino acids or the 4 base nucleic acids (DNA / RNA)",
                    video_url="https://imgur.com/1yhSVsW",
                ),
                MenuItem(
                    label="Common Elements",
                    name="Color Common",
                    description="Choose a color for each of the common elements. This is a smaller convenience node for elements which commonly appear in macromolecular structures",
                    video_url="https://imgur.com/GhLdNwy",
                ),
                Break(),
                MenuItem(
                    label="Goodsell",
                    name="Color Goodsell",
                    description="Change the inputted color to be darker for non-carbon atoms. Creates a _Goodsell Style_ color scheme for individual chains.",
                    video_url="https://imgur.com/gPgMSRa",
                ),
                MenuItem(
                    label="Rainbow",
                    name="Color Rainbow",
                    description="Generate a rainbow color palette, that changes over from start to finish along a peptide chain. Can be one rainbow over the entire structure, or create a rainbow of a per-chani basis.",
                ),
                MenuItem(
                    label="Attribute Map",
                    name="Color Attribute Map",
                    description="Interpolate between two or three colors, based on the value of an attribute field such as `b_factor`. Choosing the minimum and maximum values with the inputs.",
                    video_url="https://imgur.com/lc2o6e1",
                ),
                MenuItem(
                    label="Attribute Random",
                    name="Color Attribute Random",
                    description="Generate a random color, based on the given attribute. Control the lightness and saturation of the color with the inputs.",
                    video_url="https://imgur.com/5sMcpAu",
                ),
                MenuItem(
                    label="pLDDT",
                    name="Color pLDDT",
                    description="Assigns colors using the `b_factor` attribute, which contains the `pLDDT` attribute for models that come from AlphaFold.",
                ),
                MenuItem(
                    label="Backbone",
                    name="Color Backbone",
                    description="Color atoms by whether or not they form part of a peptide or nucleic backbone",
                ),
                MenuItem(
                    label="Secondary Structure",
                    name="Color Sec Struct",
                    description="Choose a color for the different secondary structures, based on the `sec_struct` attribute.",
                    video_url="https://imgur.com/wcJAUp9",
                ),
                Break(),
            ],
        ),
        Submenu(
            name="topology",
            title="Topology",
            items=[
                MenuItem(
                    label="DSSP",
                    name="Topology DSSP",
                    description="Calculate the secondary structure of a structure, storing it on the `sec_struct` attribute.",
                ),
                MenuItem(
                    name="Residue Mask",
                    description="Returns the index for the atom for each unique group (from res_id) for each point in that group. Allows for example, all atoms in a group to be rotated around the position of the selected atom.\n\nIn the video example, the `atom_name` is used to select an atom within the groups. Each atom's position is then offset to that position, showing the group-wise selection.",
                    video_url="https://imgur.com/sD3jRTR",
                ),
                MenuItem(
                    name="Backbone Positions",
                    description='If the atoms have been through the "Compute Backbone" node, then the backbone atom positions will be available as attributes through this node.\n\nIn the video example, the `Alpha Carbons` output is styled as spheres, where the position is mixed with some of the backbone posiitons. The backbone positions can also be selected from the AA residue higher or lower with the specified offset.',
                    video_url="https://imgur.com/6X2wnpY",
                ),
                MenuItem(
                    name="Dihedral Phi",
                    description="",
                    video_url="",
                ),
                MenuItem(
                    name="Dihedral Psi",
                    description="",
                    video_url="",
                ),
                MenuItem(
                    name="Backbone Position",
                    description="Return the backbone position for the peptide residue, and recalculate if the attribute doesn't exist",
                    video_url="",
                ),
                MenuItem(
                    name="Backbone Vectors",
                    description="Calculate `Normal`, `Tangent` and `Bitangent` values from protein backbone atom positions",
                    video_url="",
                ),
                Break(),
                MenuItem(
                    name="Atom ID",
                    description="The `atom_id` attribute which is read from the file. Will be increasing linearly for each subsequent atom",
                ),
                MenuItem(
                    name="Chain ID",
                    description="The 'chain_id' attribute, an integer representation of the Chain IDs from the structure. Chains are sorted alphabetically then assigned an ID startin at `0` and increasing.",
                ),
                MenuItem(
                    name="Atom Name",
                    description="The `atom_name` attribute, an integer representation of the atom names such as C for carbon, CA for alpha carbon",
                ),
                MenuItem(
                    name="Residue Name",
                    description="The `res_name` attribute, an integer representation of the residue names. Amino acids are sorted alphabetically and assigned a value starting at `0`. DNA starts at 40 and DNA starts at 30",
                ),
                Break(),
                MenuItem(
                    name="Res Info",
                    description="Read information about the atoms with the context of each residue the atom is in",
                ),
                MenuItem(
                    name="Chain Info",
                    description="Read information about the residues within the context of each chain",
                ),
                MenuItem(
                    name="Res Group ID",
                    description="A unique Group ID that is calculated for every residue in the structure",
                ),
                Break(),
                MenuItem(
                    label="Find Bonds",
                    name="Topology Find Bonds",
                    description="Finds bonds between atoms based on distance. Based on the vdw_radii for each point, finds other points within a certain radius to create a bond to. Does not preserve the index for the points, detect bond type, or transfer all attributes",
                    video_url="https://imgur.com/oUo5TsM",
                ),
                MenuItem(
                    label="Break Bonds",
                    name="Topology Break Bonds",
                    description="Will delete a bond between atoms that already exists based on a distance cutoff, or is selected in the `Selection` input. Leaves the atoms unaffected",
                    video_url="https://imgur.com/n8cTN0k",
                ),
                MenuItem(
                    label="Bond Count",
                    name="Bond Count",
                    description="The number of bonds for an atom",
                ),
                MenuItem(
                    label="Edge Info",
                    name="Edge Info",
                    description='Get information for the selected edge, evaluated on the point domain. The "Edge Index" selects the edge from all possible connected edges. Edges are unfortunately stored somewhat randomly. The resulting information is between the evaluating point and the point that the edge is between. Point Index returns -1 if not connected.\n\nIn the video example, cones are instanced on each point where the Edge Index returns a valid connection. The Edge Vector can be used to align the instanced cone along that edge. The length of the edge can be used to scale the cone to the other point. As the "Edge Index" is changed, the selected edge changes. When "Edge Index" == 3, only the atoms with 4 connections are selected, which in this model (1BNA) are just the phosphates.',
                    video_url="https://imgur.com/Ykyis3e",
                ),
                MenuItem(
                    label="Edge Angle",
                    name="Point Edge Angle",
                    description=' Calculate the angle between two edges, selected with the edge indices. For molecule bonds, combinations of [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] will select all possible bond angles.\n\nIn the video example, two edges are selected with their "Edge Index" values. Those atoms which aren\'t valid return false and do not get instanced. The two edge vectors are used to calculate the perpendicular vector through cross product, around which the rotation for the cone is rotated. This demonstrates the ability to calculate the edge angle between the two selected edges.',
                    video_url="https://imgur.com/oQP6Cv8",
                ),
                MenuItem(
                    label="Points of Edge",
                    name="Points of Edge",
                    description='Finds the conntected point for the selected "Edge Index", and returns each point index for all of the points connected to that point. If the connection doesn\'t exist, or the connection is back to the original point, -1 is returned.\n\nIn the video example, a new point is selected based on the "Edge Index". At that point, all of the connecting points are exposed as indices `0, 1, 2, 3`. If that index is not a valid point or connection, or the point is the same as the original point that is being evaluated, then -1 is returned. \n\nThis is one of the more complicated topology nodes, but allows indexing of the atoms that are bonded to a bonded atom. This helps with doing calculations for planar molecules.',
                    video_url="https://imgur.com/fZ6srIS",
                ),
            ],
        ),
        Submenu(
            name="curves",
            title="Curves",
            items=[
                MenuItem(
                    name="Curve Vectors",
                    description="Vectors relevant to orientations on a curve: Normal, Tangent and Bitangent",
                ),
                MenuItem(
                    name="Curve Offset Dihedral",
                    description="Calculate the Dihedral between normals along a curve, mixing the resulting angle",
                ),
                MenuItem(
                    name="Curve Transform",
                    description="The transform for the point, combining the `Position` and `Curve Rotation` using `Normal` and `Tangent`",
                ),
                MenuItem(
                    name="Curve Rotation",
                    description="The `Rotation` for the point,  using the `Normal` and `Tangent` of the point",
                ),
                MenuItem(
                    name="Curve Offset Dot",
                    description="Calculate the dot product of the current point and another, using their `Normal` attributes. Also returns the dot producted thresholded for a particular cutoff",
                ),
                MenuItem(
                    name="Offset Point Along Curve",
                    description="Return the `Factor` or `Length` of a point, when offset an amount of points along their curve. A value of 1.5 returns the `Factor` / `Length` that is half way between the two points that are +1 and +2 from the current `Point Index`",
                ),
                MenuItem(
                    name="Curve Endpoint Values",
                    description="Integer values for the ends of splines, specifying the start and end sizes",
                ),
                Break(),
                MenuItem(
                    name="Curve Visualize",
                    description="Visualize curves, instancing Gimbals with the resulting curve rotation and positions",
                ),
                MenuItem(
                    name="Curve Custom Profile",
                    description="`Curve to Mesh` but with the potential for a custom profile, with fields for custom rotations and scaling along the curve while generating the geometry",
                ),
            ],
        ),
        Submenu(
            name="geometry",
            title="Geometry",
            items=[
                MenuItem(
                    name="Split to Centred Instances",
                    description="Split points to instances, with the origin points being the `Group ID` calculated centroids",
                ),
                MenuItem(
                    name="Centre on Selection",
                    description="Move the input points to be centred on their calculated cnetroid point, which is based on the selection. The optional `Group ID` value applies this transformation on a per-group basis",
                ),
                MenuItem(
                    name="Primitive Arrow",
                    description="A simple arrow geometry",
                ),
                MenuItem(
                    name="Primitive Gimbal",
                    description="A 3-axis gimbal made of `Primitive Arrow`s, useful for visualisation and debugging",
                ),
            ],
        ),
        Submenu(
            name="fields",
            title="Fields",
            description="For working with and manipulating fields. Evaluating at specific indices, offsetting, fallback values and picking from groups.",
            items=[
                MenuItem(
                    label="Index Mixed",
                    name="Index Mixed",
                    description="Offset the current point's `Index` by a given amount, returning the fractional Index value",
                ),
                MenuItem(
                    label="Index Mix Float",
                    name="Index Mix Float",
                    description="Sample and interpolate the `Float` value between the floor and ceiling of the given `Index`, by the fractional amount",
                ),
                MenuItem(
                    label="Index Mix Vector",
                    name="Index Mix Vector",
                    description="Sample and interpolate the `Vector` value between the floor and ceiling of the given `Index`, by the fractional amount",
                ),
                MenuItem(
                    label="Index Mix Rotation",
                    name="Index Mix Rotation",
                    description="Sample and interpolate the `Rotation` value between the floor and ceiling of the given `Index`, by the fractional amount",
                ),
                MenuItem(
                    label="Index Mix Color",
                    name="Index Mix Color",
                    description="Sample and interpolate the `Color` value between the floor and ceiling of the given `Index`, by the fractional amount",
                ),
                Break(),
                MenuItem(
                    label="Sample Position",
                    name="Sample Position",
                    description="Sample the `Position` attribute from a point at the given `Index`",
                ),
                MenuItem(
                    label="Sample Mixed Float",
                    name="Sample Mixed Float",
                    description="Sample the `Float` attribute, mixed between the floor and ceiling of the given `Index`",
                ),
                MenuItem(
                    label="Sample Mixed Vector",
                    name="Sample Mixed Vector",
                    description="Sample the `Vector` attribute, mixed between the floor and ceiling of the given `Index`",
                ),
                MenuItem(
                    label="Sample Mixed Rotation",
                    name="Sample Mixed Rotation",
                    description="Sample the `Rotation` attribute, mixed between the floor and ceiling of the given `Index`",
                ),
                MenuItem(
                    label="Sample Mixed Color",
                    name="Sample Mixed Color",
                    description="Sample the `Color` attribute, mixed between the floor and ceiling of the given `Index`",
                ),
                Break(),
                MenuItem(
                    label="Offset Integer",
                    name="Offset Integer",
                    description="Evaluate an `Integer` at an index that is offset by the specified amount",
                ),
                MenuItem(
                    label="Offset Float",
                    name="Offset Float",
                    description="Evaluate a `Float` value at an index that is offset by the specified amount",
                ),
                MenuItem(
                    label="Offset Vector",
                    name="Offset Vector",
                    description="Evaluate a `Vector` at an index that is offset by the specified amount",
                ),
                MenuItem(
                    label="Offset Boolean",
                    name="Offset Boolean",
                    description="Evaluate a `Boolean` at an index that is offset by the specified amount",
                ),
                MenuItem(
                    label="Offset Rotation",
                    name="Offset Rotation",
                    description="Evaluate a `Rotation` at an index that is offset by the specified amount",
                ),
                MenuItem(
                    name="Offset Matrix",
                    description="Evaluate a `4X4 Matrix` at an index that is offset by the specified amount",
                ),
                Break(),
                MenuItem(
                    name="Fallback Float",
                    description="Use the `Float` attribute specified by name. If the attribute doesn't exist, use the `Fallback` value instead",
                ),
                MenuItem(
                    name="Fallback Vector",
                    description="Use the `Vector` attribute specified by name. If the attribute doesn't exist, use the `Fallback` value instead",
                ),
                MenuItem(
                    name="Fallback Integer",
                    description="Use the `Integer` attribute specified by name. If the attribute doesn't exist, use the `Fallback` value instead",
                ),
                MenuItem(
                    name="Fallback Boolean",
                    description="Use the `Boolean` attribute specified by name. If the attribute doesn't exist, use the `Fallback` value instead",
                ),
                MenuItem(
                    name="Fallback Color",
                    description="Use the `Color` attribute specified by name. If the attribute doesn't exist, use the `Fallback` value instead",
                ),
                MenuItem(
                    name="Fallback Rotation",
                    description="Use the `Rotation` attribute specified by name. If the attribute doesn't exist, use the `Fallback` value instead",
                ),
                MenuItem(
                    name="Fallback Matrix",
                    description="Use the `Matrix` attribute specified by name. If the attribute doesn't exist, use the `Fallback` value instead",
                ),
                Break(),
                MenuItem(
                    label="Group Pick",
                    name="Group Pick",
                    description="For each group, return the index of the point for which the Selection is true. Only valid if there is a single true in the group. If not lvalid, returns -1",
                ),
                MenuItem(
                    label="Group Pick Vector",
                    name="Group Pick Vector",
                    description="For each group, return the Position of the point at which the selection is true. If there is more than one true for the group the pick is not valid and (0, 0, 0) is returned",
                ),
                MenuItem(
                    label="Group Info",
                    name="Group Info",
                    description="Based on the Group ID input, return the size of the group and the indices of the first and last items of the group",
                ),
                Break(),
                MenuItem(
                    label="Attribute Run",
                    name="Attribute Run",
                    description="Fill in gaps in a set of continuous boolean True values, up to a specific size",
                ),
                MenuItem(
                    label="Integer Run",
                    name="Integer Run",
                    description="Group mask output increments whenever the input value or the Group ID changes",
                ),
                MenuItem(
                    label="Boolean Run Fill",
                    name="Boolean Run Fill",
                    description="Fill in gaps in a set of continuous boolean True values, up to a specific size",
                ),
                MenuItem(
                    label="Boolean Run Mask",
                    name="Boolean Run Mask",
                    description="Mask a run of boolean values. Potentially trim the start or ending values and specifying a minimum length under which they are considered false",
                ),
            ],
        ),
        Submenu(
            name="ensemble",
            title="Ensemble",
            items=[
                MenuItem(
                    label="Biological Assembly",
                    name="mn.assembly_bio",
                    backup="MN_assembly_",
                    description="Creates a biological assembly by applying rotation and translation matrices to individual chains in the structure. It is created on an individual molecule basis, if assembly instructions are detected when imported.\n\n::: callout-caution \n\nStyle Spheres requires the material to be set again after the assembly node, as the material is currently lost when joining multiple point clouds.\n\n:::",
                    video_url="https://imgur.com/TNc102v",
                ),
                MenuItem(
                    label="Center Assembly",
                    name="MN_assembly_center",
                    description="Move an instanced assembly to the world origin. Some structures are not centred on the world origin, so this node can reset them to the world origin for convenient rotation and translation and animation.",
                    video_url="https://imgur.com/pgFTmgC",
                ),
                MenuItem(
                    label="Instance",
                    name="Ensemble Instance",
                    description="Instance the items of an ensemble onto the given points",
                ),
            ],
        ),
        Submenu(
            name="DNA",
            title="DNA",
            items=[
                MenuItem(
                    label="Double Helix",
                    name="MN_dna_double_helix",
                    description="Create a DNA double helix from an input curve.\nTakes an input curve and instances for the bases, returns instances of the bases in a double helix formation",
                ),
                MenuItem(
                    label="Bases",
                    name="MN_dna_bases",
                    description="Provide the DNA bases as instances to be styled and passed onto the Double Helix node",
                ),
                Break(),
                MenuItem(
                    label="Style Spheres Cycles",
                    name="MN_dna_style_spheres_cycles",
                    description="Style the DNA bases with spheres only visible in Cycles",
                ),
                MenuItem(
                    label="Style Spheres EEVEE",
                    name="MN_dna_style_spheres_eevee",
                    description="Style the DNA bases with spheres visible in Cycles and EEVEE",
                ),
                MenuItem(
                    label="Style Surface",
                    name="MN_dna_style_surface",
                    description="Style the DNA bases with surface representation",
                ),
                MenuItem(
                    label="Style Ball and Stick",
                    name="MN_dna_style_ball_and_stick",
                    description="Style the DNA bases with ball and stick representation",
                ),
            ],
        ),
        Submenu(
            name="animate",
            title="Animate",
            items=[
                MenuItem(
                    label="Animate Frames",
                    name="Animate Frames",
                    description="Animate the atoms of a structure, based on the frames of a trajectory from the `Frames` collection in the input. The structure animates through the trajectory from the given start frame to the given end frame, as the `Animate 0..1` value moves from `0` to `1`. Values higher than `1` start at the beginning again and the trajectory will loop repeating every `1.00`.\nPosition and `b_factor` are interpolated if available. By default linear interpolation is used. Smoothing in and out of each frame can be applied with the `Smoother Step`, or no interpolation at all.",
                    video_url="https://imgur.com/m3BPUxh",
                ),
                MenuItem(
                    label="Animate Value",
                    name="Animate Value",
                    description="Animate a float value between the specified min and max values, over specified range of frames. If clamped, frames above and below the start and end will result in the min and max output values, otherwise it will continue to linearly interpolate the value beyond the min and max values.",
                    video_url="https://imgur.com/2oOnwRm",
                ),
                Break(),
                MenuItem(
                    label="Animate Trails",
                    name="Animate Trails",
                    description="Add trails to the atoms as they are animated, which trail the specified number of frames behind the atoms",
                ),
                Break(),
                MenuItem(
                    label="Res Wiggle",
                    name="Animate Wiggle",
                    description="Create a procedural animation of side-chain movement. 'Wiggles' the side-chains of peptide amino acids based on the `b_factor` attribute. Wiggle is currently only supported for protein side-chains and does not check for steric clashing so higher amplitudes will result in strange results. The animation should seamlessly loop every `1.00` of the `Animate 0..1` input.",
                    video_url="https://imgur.com/GK1nyUz",
                ),
                MenuItem(
                    label="Peptide to Curve",
                    name="Animate Peptide to Curve",
                    description="Take the protein residues from a structure and align then along an input curve. Editing the curve will change how the atoms are arranged. The output atoms can be styled as normal.",
                    video_url="https://imgur.com/FcEXSZx",
                ),
                Break(),
                MenuItem(
                    label="Noise Position",
                    name="MN_animate_noise_position",
                    description="Create 3D noise vector based on the position of points in 3D space. Evolve the noise function with the `Animate` input, and change the characteristics of the noise function with the other inputs such as scale and detail. There is also a 1-dimensional noise output called `Fac`.\n\nAn example of using this noise is to offset the positions of atoms with the `Set Position` node.",
                    video_url="https://imgur.com/B8frW1C",
                ),
                MenuItem(
                    label="Noise Field",
                    name="MN_animate_noise_field",
                    description="Create a 3D noise vector based on the input field. Evolve the noise function with the `Animate` input, and change the characteristics of the noise function with the other inputs such as scale and detail. There is also a 1-dimensional noise output called `Fac`.\n\nAn example of using this noise is to offset the positions of atoms with the `Set Position` node. Different field inputs result in different noise being applied. Using the `chain_id` results in the same noise being generated for each atom in each chain, but different between chains.",
                    video_url="https://imgur.com/hqemVQy",
                ),
                MenuItem(
                    label="Noise Repeat",
                    name="MN_animate_noise_repeat",
                    description="Create a 3D noise vector based on the input field, that repeats every `1.00` for the `Animate 0..1` input. Evolve the noise function with the `Animate` input, and change the characteristics of the noise function with the other inputs such as scale and detail. There is also a 1-dimensional noise output called `Fac`.\n\nAn example of using this noise is to offset the positions of atoms with the `Set Position` node. Different field inputs result in different noise being applied. Using the `chain_id` results in the same noise being generated for each atom in each chain, but different between chains.",
                    video_url="https://imgur.com/GNQcIlx",
                ),
            ],
        ),
        Submenu(
            name="utils",
            title="Utilities",
            items=[
                MenuItem(
                    label="Curve Resample",
                    name="MN_utils_curve_resample",
                    description="",
                ),
                MenuItem(
                    name="Attribute Remap",
                    description="Sample an attribute from the mesh and remap from the minimum to the maximum to the specified values",
                ),
                MenuItem(
                    name="Field Remap",
                    description="Sample a field from the mesh and remap from the minimum to the maximum to the specified values",
                ),
                Break(),
                MenuItem(
                    label="Between Integer",
                    name="Between Integer",
                    description="Test if an integer is between (and including) the upper and lower bounds",
                ),
                MenuItem(
                    label="Between Float",
                    name="Between Float",
                    description="Test if a float is between the upper and lower bounds",
                ),
                MenuItem(
                    label="Between Vector",
                    name="Between Vector",
                    description="Test if a vector is element-wise between the upper and lower bounds.",
                ),
                Break(),
                MenuItem(
                    name="Vector from Point",
                    description="Calculate the vector from the current point's position to the input vector",
                ),
                MenuItem(
                    name="Mix Position",
                    description="Mix the current point's Position with the input vector",
                ),
                Break(),
                MenuItem(
                    name="Fractionate Float",
                    description="Test if a vector is element-wise between the upper and lower bounds.",
                ),
                Break(),
                MenuItem(
                    label="Centroid",
                    name="Centroid",
                    description="Calculate the centroid point for the selection for each group in the `Group ID`",
                ),
                MenuItem(
                    label="Vector Angle",
                    name="Vector Angle",
                    description="Compute the angle in radians between two vectors.",
                ),
                MenuItem(
                    label="Dihedral Angle",
                    name="Dihedral Angle",
                    description='Computes the angle between two vectors, AB & CD around around the axis of BC. The first vector AB is treated as the "12 O\'clock" up position, looking down the axis towards C, with angles being return in the range of (-Pi, Pi). Clockwise angles are positive and anti-clockwise angles are negative.',
                    video_url="",
                ),
                MenuItem(
                    label="3 Point Angle",
                    name="3 Point Angle",
                    description="Calculate the angle between 3 different points. These points are selected based on their index in the point domain, with Index B being the centre of the calculation.\n\nIn the video example, the same calculation that is occurring internally inside of the `MN_topo_edge_angle` node, is being handled explicity by this node. If the `Index` is being used as `Index B` then the current point that is being evaluated is the centre of the angle calculation. If this value is changed, then the point at the corresponding index is used, which results in a smaller angle in the example video.",
                    video_url="https://imgur.com/qXyy2ln",
                ),
                MenuItem(
                    label="2 Point Angle",
                    name="2 Point Angle",
                    description="Calculate the angle that two points make, relative to the current point being evaluated. Points are selected based on their index, with the centre of the angle calculation being the current point's position. Equivalent to using 3-Point angle and using `Index` as the `Index B`.\n\nIn the example video, the angle calculation is similar to that of the 3-Point Angle node, but the middle point is always the current point.",
                    video_url="https://imgur.com/xp7Vbaj",
                ),
                MenuItem(
                    label="Point Distance",
                    name="Point Distance",
                    description="Calculate the distance and the vector between the evaluating point and the point selected via the Index.\n\nIn the example video, each point is calculating a vector and a distance between itself and the indexed point. When the Point Mask node is used, this index is then on a per-group basis, so each point in the group points to just the group's corresponding point.",
                    video_url="https://imgur.com/AykNvDz",
                ),
                MenuItem(
                    label="Cartoon Utilities",
                    name=".MN_utils_style_cartoon",
                    description="The underlying node group which powers the cartoon style",
                ),
            ],
        ),
        Submenu(
            name="density",
            title="Density",
            items=[
                MenuItem(
                    label="Style Surface",
                    name="Style Density Surface",
                    description="A surface made from the electron density given a certain threshold value.",
                    video_url="https://imgur.com/jGgMSd4",
                ),
                MenuItem(
                    label="Style Wire",
                    name="Style Density Wire",
                    description="A wire surface made from the electron density given a certain threshold value.",
                    video_url="https://imgur.com/jGgMSd4",
                ),
                MenuItem(
                    label="Sample Nearest Attribute",
                    name="Sample Nearest Atoms",
                    description="Sample the nearest atoms from another object, to get the colors or other attributes and apply them to a volume mesh.",
                    video_url="https://imgur.com/UzNwLv2",
                ),
            ],
        ),
    ]
)
