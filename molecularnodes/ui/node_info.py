from .menu import Menu, Submenu, MenuItem, CustomItem, Break

menu_items = Menu(
    submenus=[
        Submenu(
            name="style",
            title="Style",
            items=[
                MenuItem(
                    label="Spheres",
                    name="Style Spheres",
                    description="Space filling atomic spheres scaled based on their `vdw_radii` attribute. By default they are _Point Clouds_ which don't render well inside of EEVEE, but are extremely fast to render and manipulate inside of Cycles. To show the atoms inside of EEVEE you can enable the `As Mesh` input which realises them into real geometry",
                    videos="https://imgur.com/eRDo2As",
                ),
                MenuItem(
                    label="Cartoon",
                    name="Style Cartoon",
                    description="Peptide helices and sheets are emphasized with the 'traditional' cartoon style. Secondary structure becomes distinct and nucleic acids become a simplified representation. The node also includes an option to calculate secondary structure by enabling the `DSSP` input, but this will differ to secondary structure imported from structures from the PDB.",
                    videos="https://imgur.com/jFdzd5J",
                ),
                MenuItem(
                    label="Ribbon",
                    name="Style Ribbon",
                    description="A simplified tube that through the alpha carbons of the protein. Controls for the raadius and smoothing of sheets are available.",
                    videos="https://imgur.com/LROZR8i",
                ),
                MenuItem(
                    label="Surface",
                    name="Style Surface",
                    description="The `vdw_radii` and a probe are used to calculate the surface of the molecule. This is still a close approximation of other surface generation algorithms and won't match those made by PyMol and ChimeraX exactly. Coloring of the surface is by default done by sampling the closest alpha carbon, but this can be disabled to use the closest atom for coloring the mesh",
                    videos="https://imgur.com/GM8TYZm",
                ),
                MenuItem(
                    label="Ball and Stick",
                    name="Style Ball and Stick",
                    description="Shows the atoms and bonds. The bonds are cylinders between the atoms, which can be split apart for double and triple bonds if the information is within the structure on import. Spheres can be displayed as mesh of as _Point Cloud_ as with the `Style Sheres` node. Bonds can also be calculated if none are present, but this is approximate and purely distance based so should not be relied upon",
                    videos="https://imgur.com/nXPQN7W",
                ),
                MenuItem(
                    label="Sticks",
                    name="Style Sticks",
                    description="Each bond terminates cleanly in a half sphere for a typical stick representation of the atoms and bonds",
                    videos="https://imgur.com/iNqE87M",
                ),
                Break(),
                MenuItem(
                    label="Preset 1",
                    name="Style Preset 1",
                    description="Quickly switch between several different pre-made preset styles. Best used when using MolecularNodes via scripts, ensuring all atoms are displayed using a combination of cartoons and atoms.",
                    videos="https://imgur.com/gCQRWBk.mp4",
                ),
                MenuItem(
                    label="Preset 2",
                    name="Style Preset 2",
                    description="Quickly switch between several different pre-made preset styles. Best used when using MolecularNodes via scripts, ensuring all atoms are displayed using a combination of cartoons and atoms.",
                    videos="https://imgur.com/gCQRWBk.mp4",
                ),
                MenuItem(
                    label="Preset 3",
                    name="Style Preset 3",
                    description="Quickly switch between several different pre-made preset styles. Best used when using MolecularNodes via scripts, ensuring all atoms are displayed using a combination of cartoons and atoms.",
                    videos="https://imgur.com/gCQRWBk.mp4",
                ),
                MenuItem(
                    label="Preset 4",
                    name="Style Preset 4",
                    description="Quickly switch between several different pre-made preset styles. Best used when using MolecularNodes via scripts, ensuring all atoms are displayed using a combination of cartoons and atoms.",
                    videos="https://imgur.com/gCQRWBk.mp4",
                ),
            ],
        ),
        Submenu(
            name="select",
            title="Select",
            items=[
                MenuItem(
                    name="Separate Atoms",
                    description="Separate the input atoms intwo the `Atoms` and `Inverted` based on the input selection. For simply styling this can the same as inputting the selecting directly into the Style node. The one additional output for this node is the `Index` field, which is the index of the atom in the original structure before the selection happened and potentially changed it's index",
                    videos="https://imgur.com/OAcekhf",
                ),
                MenuItem(
                    name="Separate Polymers",
                    description="Separate the input atomic geometry into it's different polymers or `Protein`, `Nucleic Acid` and `other`.",
                    videos="https://imgur.com/trg0voP",
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
                    videos="https://imgur.com/U7HqDct",
                ),
                CustomItem(
                    label="Entity",
                    field="entity_id",
                    name="Select Entity_",
                    dtype="BOOLEAN",
                    prefix="",
                    property_id="entity_ids",
                    description="Select single or multiple of the different entities. Creates a selection based on the `entity_id` attribute.",
                    videos="https://imgur.com/h5KZBXt",
                ),
                CustomItem(
                    label="Ligand",
                    field="res_name",
                    name="Select Ligand_",
                    dtype="BOOLEAN",
                    prefix="",
                    property_id="ligands",
                    description="Select single or multiple of the different ligands.",
                    videos="https://imgur.com/s2seWIw",
                ),
                CustomItem(
                    label="Segment",
                    field="segid",
                    name="Select Segment_",
                    dtype="BOOLEAN",
                    prefix="",
                    property_id="segments",
                    description="Output a selection based on the `segment_id` that is present in some MD trajectories",
                ),
                MenuItem(
                    label="Atomic Number",
                    name="Select Atomic Number",
                    description="Select points based on their `atomic_number` attribute, corresponding to the element's atomic number. Useful for selecting single elements quickly",
                    videos="https://imgur.com/B6V9W3F",
                ),
                MenuItem(
                    label="Element",
                    name="Select Element",
                    description="Select points for the first 80 elemnts of the preiodic table, via a boolean input. Elements are grouped into panels of 20 each for orgnaisation and convenience",
                    videos="https://imgur.com/d6Q3T7D",
                ),
                MenuItem(
                    label="Res Name",
                    name="Select Res Name",
                    description="Select points based on their `res_id` attribute for different peptide or nucleic acid residue names. Inputs are arranged alphabetically and in panels of _Protein_, _RNA_ and _DNA_ for layout",
                    videos="https://imgur.com/smwnKsL",
                ),
                MenuItem(
                    label="Res ID Single",
                    name="Select Res ID",
                    description="Select a atoms based on their `res_id` number.",
                    videos="https://imgur.com/BL6AOP4",
                ),
                MenuItem(
                    label="Res ID Range",
                    name="Select Res ID Range",
                    description="Select multiple residues by specifying a _Min_ and a _Max_, defining a range that includes or excludes based on the `res_id` number.",
                    videos="https://imgur.com/NdoQcdE",
                ),
                MenuItem(
                    label="Res ID",
                    name="mn.residues_selection_custom",
                    backup="Select Res ID_",
                    description="Create a more complex selection for the `res_id` field, by specifying multiple ranges and potential single `res_id` numbers. This node is built uniquely each time, to the inputs will look different for each user.\nIn the example below, residues 10 & 15 are selected, as well as residues between and including 20-100.\nThe node was created by inputting `10, 15, 20-100` into the node creation field.",
                    videos="https://imgur.com/OwAXsbG",
                ),
                Break(),
                MenuItem(
                    label="Attribute",
                    name="Select Attribute",
                    description="Select atoms that have true for the given attribute name.",
                ),
                MenuItem(
                    name="Is Peptide",
                    description="Outputs a selection for all of the points in a peptide",
                    videos="https://imgur.com/KGU2ElV",
                ),
                MenuItem(
                    name="Is Nucleic",
                    description="Outputs a selection for all of the points in a nucleic acid molecule",
                    videos="https://imgur.com/UBVBPnI",
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
                    description="Outputs a selection for all of the points that are alpha carbons (CA) in the structure",
                    videos="https://imgur.com/JxGa9Ou",
                ),
                MenuItem(
                    name="Is Backbone",
                    description="Outputs a selection for the backbone atoms of a peptide or nucleic acid polymer. Peptide backbone includes the peptide oxygen and the alpha carbon as part of the backbone",
                    videos="https://imgur.com/MBpgktt",
                ),
                MenuItem(
                    name="Is Side Chain",
                    description="Outputs a selection for the side chain atoms of a peptide or nucleic acid polymer. Peptide side chain can optionally include or exclude the alpha carbon CA of the side chain",
                    videos="https://imgur.com/Hpy5AQc",
                ),
                MenuItem(
                    name="Is Helix",
                    description="Outputs a selection for points that are part of a helix. The `sec_struct` attribute must exist, either imported from the file or computed with the `DSSP` node",
                    videos="https://imgur.com/Rp2CPvq",
                ),
                MenuItem(
                    name="Is Sheet",
                    description="Outputs a selection for points that are part of a sheet. The `sec_struct` attribute must exist, either import from the file or computed with the `DSSP` node",
                    videos="https://imgur.com/OObFAzq",
                ),
                MenuItem(
                    name="Is Loop",
                    description="Outputs a selection for those points in a peptide which are not part of any secondary structure (loop or helix)",
                    videos="https://imgur.com/Buy0OEu",
                ),
                Break(),
                MenuItem(
                    label="Bonded",
                    name="Select Bonded",
                    description="Based on an initial selection, finds atoms which are within a certain number of bonds of this selection. Output can include or excluded the original selection.",
                    videos="https://imgur.com/g8hgXup",
                ),
                MenuItem(
                    label="Res Whole",
                    name="Select Res Whole",
                    description="Expand the given selection to include a whole residue, if a single atom in that residue is selected. Useful for when a distance or proximity selection includes some of the residue and you wish to include all of the residue.",
                    videos="https://imgur.com/JFzwE0i",
                ),
                Break(),
                MenuItem(
                    label="Cube",
                    name="Select Cube",
                    description="Create a selection that is inside the `Empty_Cube` object. When this node is first created, an _empty_ object called `Empty_Cube` should be created. You can always create additional empty objects through the add menu, to use a different object. The rotation and scale of the object will be taken into account for the selection.",
                    videos="https://imgur.com/P4GZ7vq",
                ),
                MenuItem(
                    label="Sphere",
                    name="Select Sphere",
                    description="Create a selection that is within a spherical radius of an object, based on that object's scale. By default an _empty_ object called `Empty_Sphere` is created. You can use other objects or create a new empty to use. The origin point for the object will be used, which should be taken in to account when using molecules. Use [`MN_select_proximity`](#select-proximity) for selections which are within a certain distance of a selection of atoms instead of a single origin point.",
                    videos="https://imgur.com/xdeTZR7",
                ),
                MenuItem(
                    label="Proximity",
                    name="Select Proximity",
                    description="Create a selection based on the proximity to the Target Atoms of the input. A sub-selection of the Target atoms can be used if the `Subset` input is used. You can expand the selection to include an entire residue if a single atom in that residue is selected, by setting `Expand` to `True`.\nIn the example below, the `Style Atoms` is being applied to a selection, which is being calculated from the proximity of atoms to specific chains. As the cutoff for the selection is changed, it includes or excludes more atoms. The `Whole Residue` option also ensures that entire residues are shown.",
                    videos="https://imgur.com/QmvDZn2",
                ),
            ],
        ),
        Submenu(
            name="color",
            title="Color",
            description="Customise the colors of the structure using these nodes. The color for atoms are always set before using the `Style` nodes. The `Style Cartoon`, `Style Surface` and `Style Ribbon` nodes take their color information from the atoms that are used to generate them. There are nodes to generate colors for each chain, each residue. The setting of `Color` attribute is done using the `Set Color` attribute.",
            items=[
                MenuItem(
                    name="Set Color",
                    description="Set the `Color` attribute on the point domain, based on the input `Color` field. The `Selection` field limits the setting of the color to the input `Selection`, with the unselected points keeping their current `Color` attribute",
                    videos="https://imgur.com/1vTKUom",
                ),
                Break(),
                CustomItem(
                    label="Chain",
                    field="chain_id",
                    dtype="RGBA",
                    name="Color Chain_",
                    prefix="",
                    property_id="chain_ids",
                    description="Pick colors that are used for each `chain_id` that is present in the structure. The inputs can also accept fields to generate colors for that chain based on upstream field evaluation. This node is generated on a per-structure basis, so will only be useful for the corresponding structure.",
                    videos="https://imgur.com/F6f8MQW",
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
                    description="Choose the colors for the individual entities in the structure. If there are multiple copies of the same protein or molecule in the structure but in different conformations, it can helpt to color them based on their `entity_id` rather than individual `chain_id` values",
                    videos=[
                        "https://imgur.com/dqa9Jro",
                        "https://imgur.com/EMmQsEP",
                    ],
                ),
                CustomItem(
                    label="Ligand",
                    field="res_name",
                    name="Color Ligand_",
                    dtype="RGBA",
                    prefix="",
                    property_id="ligands",
                    description="Choose the colors for individual ligands in the structure.",
                    videos="https://imgur.com/bQh8Fd9",
                ),
                MenuItem(
                    label="Element",
                    name="Color Element",
                    description="Choose colors for the individual elements of the periodic table. The first 80 elements are listed, grouped into panels in chunks of 20 elements per panel for better organisation.",
                    videos="https://imgur.com/MHjCHJd",
                ),
                MenuItem(
                    label="Atomic Number",
                    name="Color Atomic Number",
                    description="Choose the color for an individual `atomic_number`. Useful for joining with fields and selecting individual elements for coloring",
                    videos="https://imgur.com/tErqBNb",
                ),
                MenuItem(
                    label="Res Name",
                    name="Color Res Name",
                    description="Choose colors for the individual amino acid and nucleic acid residues",
                    videos="https://imgur.com/7bEmJGU",
                ),
                MenuItem(
                    label="Common Elements",
                    name="Color Common",
                    description="For convenience, choose colors for the common 6 elements found protein structures",
                    videos="https://imgur.com/JZD0IvD",
                ),
                Break(),
                MenuItem(
                    label="Goodsell",
                    name="Color Goodsell",
                    description="Darken the color of non-carbon atoms, a similar strategy to the style developed by David Goodsell",
                    videos="https://imgur.com/RDykyT3",
                ),
                MenuItem(
                    label="Rainbow",
                    name="Color Rainbow",
                    description="Output a rainbow color spectrum along the structure. The rainbow can be generated along each chain in the structure, or along the entire structure. You can offset the strat of the rainbow and chain the color & satuation of the generated colors",
                    videos="https://imgur.com/8fYEqQn",
                ),
                MenuItem(
                    label="Attribute Map",
                    name="Color Attribute Map",
                    description="Map a value to interpolate between two colors, optionally using an intermediate color. User can specify the _Attribute Name_ that is used for the inerpolation, specifying the start and end values of the interpolation",
                    videos="https://imgur.com/ZWCSx19",
                ),
                MenuItem(
                    label="Attribute Random",
                    name="Color Attribute Random",
                    description="Generate a random color, based on the given attribute. Control the lightness and saturation of the color with the inputs.",
                    videos="https://imgur.com/5sMcpAu",
                ),
                MenuItem(
                    label="pLDDT",
                    name="Color pLDDT",
                    description="Uses the `b_factor` attribute to assign the typical AlphaFold `pLDDT` color palette. Node is automatically used when fetching directly from the AlphaFold database",
                    videos="https://imgur.com/3RAi1rm",
                ),
                MenuItem(
                    label="Backbone",
                    name="Color Backbone",
                    description="Select colors for atoms that are part of the backbone or the side chain of a peptide or nucleic acid polymer. Atoms that are neither maintain their current `Color` attribute",
                    videos="https://imgur.com/nadPkJd",
                ),
                MenuItem(
                    label="Secondary Structure",
                    name="Color Sec Struct",
                    description="Choose a color for the different secondary structures, based on the `sec_struct` attribute.",
                    videos="https://imgur.com/wcJAUp9",
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
                    videos="https://imgur.com/sD3jRTR",
                ),
                MenuItem(
                    name="Backbone Positions",
                    description='If the atoms have been through the "Compute Backbone" node, then the backbone atom positions will be available as attributes through this node.\n\nIn the video example, the `Alpha Carbons` output is styled as spheres, where the position is mixed with some of the backbone posiitons. The backbone positions can also be selected from the AA residue higher or lower with the specified offset.',
                    videos="https://imgur.com/6X2wnpY",
                ),
                MenuItem(
                    name="Dihedral Phi",
                    description="",
                    videos="",
                ),
                MenuItem(
                    name="Dihedral Psi",
                    description="",
                    videos="",
                ),
                MenuItem(
                    name="Rotate Backbone",
                    description="Rotate the atoms cumulatively for each residue, adjusting the `phi` and `psi` angles for the selected residues",
                    videos="",
                ),
                # MenuItem(
                #     name="Backbone Position",
                #     description="Return the backbone position for the peptide residue, and recalculate if the attribute doesn't exist",
                #     videos="",
                # ),
                MenuItem(name="Backbone N"),
                MenuItem(name="Backbone CA"),
                MenuItem(name="Backbone C"),
                MenuItem(name="Backbone O"),
                MenuItem(name="Backbone NH"),
                MenuItem(
                    name="Backbone Vectors",
                    description="Calculate `Normal`, `Tangent` and `Bitangent` values from protein backbone atom positions",
                    videos="",
                ),
                Break(),
                MenuItem(
                    name="Chain Group ID",
                    description="Assumes only CA points in the geometry. Unique Group ID for each chain, incrementing if distance between CA points are greater than threshold",
                ),
                MenuItem(
                    name="Chain ID",
                    description="The 'chain_id' attribute, an integer representation of the Chain IDs from the structure. Chains are sorted alphabetically then assigned an ID startin at `0` and increasing.",
                ),
                MenuItem(
                    name="Atom ID",
                    description="The `atom_id` attribute which is read from the file. Will be increasing linearly for each subsequent atom",
                ),
                MenuItem(
                    name="Entity ID",
                    description="The `entity_id` attribute of the point",
                ),
                MenuItem(
                    name="Residue ID",
                    description="The `residue_id` attribute of the point, which is the assigned number of the residue inside the chain of the structure",
                ),
                MenuItem(
                    name="Atomic Number",
                    description="The `atomic_number` attribute, ascending from `1` for each element on the periodic table",
                ),
                MenuItem(
                    name="Atom Name",
                    description="The `atom_name` attribute, an integer representation of the atom names such as C for carbon, CA for alpha carbon",
                ),
                MenuItem(
                    name="Residue Name",
                    description="The `res_name` attribute, an integer representation of the residue names. Amino acids are sorted alphabetically and assigned a value starting at `0`. DNA starts at 40 and DNA starts at 30",
                ),
                MenuItem(
                    name="Secondary Structure",
                    description="The `sec_struct` attribute of the point. `0` is non-peptide, 1 is helix, 2 is sheet and 3 is loop",
                ),
                MenuItem(
                    name="VDW Radii",
                    description="The `vdw_radii` attribute of the point, corresponding the radius of the element as defined in the data dictionary. Used for scaling sphers for `Style Spheres` and `Style Ball and Stick`",
                ),
                MenuItem(
                    name="Mass",
                    description="The `mass` attribute of the point, used for centre of mass calculations",
                ),
                MenuItem(
                    name="B Factor",
                    description="The `b_factor` attribute of the point, representing the 'temperature' of the point",
                ),
                MenuItem(
                    name="Color",
                    description="The `Color` attribute of the point, used for coloring the final generated geometry inside of the materials",
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
                    videos="https://imgur.com/oUo5TsM",
                ),
                MenuItem(
                    label="Break Bonds",
                    name="Topology Break Bonds",
                    description="Will delete a bond between atoms that already exists based on a distance cutoff, or is selected in the `Selection` input. Leaves the atoms unaffected",
                    videos="https://imgur.com/n8cTN0k",
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
                    videos="https://imgur.com/Ykyis3e",
                ),
                MenuItem(
                    label="Edge Angle",
                    name="Point Edge Angle",
                    description=' Calculate the angle between two edges, selected with the edge indices. For molecule bonds, combinations of [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] will select all possible bond angles.\n\nIn the video example, two edges are selected with their "Edge Index" values. Those atoms which aren\'t valid return false and do not get instanced. The two edge vectors are used to calculate the perpendicular vector through cross product, around which the rotation for the cone is rotated. This demonstrates the ability to calculate the edge angle between the two selected edges.',
                    videos="https://imgur.com/oQP6Cv8",
                ),
                MenuItem(
                    label="Points of Edge",
                    name="Points of Edge",
                    description='Finds the conntected point for the selected "Edge Index", and returns each point index for all of the points connected to that point. If the connection doesn\'t exist, or the connection is back to the original point, -1 is returned.\n\nIn the video example, a new point is selected based on the "Edge Index". At that point, all of the connecting points are exposed as indices `0, 1, 2, 3`. If that index is not a valid point or connection, or the point is the same as the original point that is being evaluated, then -1 is returned. \n\nThis is one of the more complicated topology nodes, but allows indexing of the atoms that are bonded to a bonded atom. This helps with doing calculations for planar molecules.',
                    videos="https://imgur.com/fZ6srIS",
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
                    description="Calculate the Dihredral angle between the current point and another point offset along the curve, using their `Normal`s",
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
                    name="Cumulative Length",
                    description="The length of the current point, added onto the cumulative length with all previous splines in the curve",
                ),
                MenuItem(
                    name="Curve Offset Dot",
                    description="Calculate the dot product of the current point and another, using their `Normal` attributes. Also returns the dot producted thresholded for a particular cutoff",
                ),
                MenuItem(
                    name="Offset Point Along Curve",
                    description="Return the `Factor` or `Length` of a point, when offset an amount of points along their curve. A value of 1.5 returns the `Factor` / `Length` that is half way between the two points that are +1 and +2 from the current `Point Index`,\n\nUseful for sampling the same curve, at a different point along the curve. Helps with interpolating along a bezier curve without having to evaluate to the intermediate points first",
                    videos=[
                        "https://imgur.com/NX6dZXR",
                        "https://imgur.com/05AmdjU",
                    ],
                ),
                MenuItem(
                    name="Curve Endpoint Values",
                    description="Integer values for the ends of splines, specifying the start and end sizes",
                ),
                Break(),
                MenuItem(
                    name="Curve Visualize",
                    description="Visualize curves, instancing Gimbals with the resulting curve rotation and positions",
                    videos="https://imgur.com/KKY7v12",
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
                    description="Same as `Split to Instances`, but with the origin point of the instance being the calculated `Centroid` for each `Group ID`. The `Selection` determines the points which contribute to the `Centroid` calculation, but all points are still offset and split into their instances",
                    videos="https://imgur.com/NcHuxsz",
                ),
                MenuItem(
                    name="Centre on Selection",
                    description="Offsets the input points so that their calculated `Centroid` is on the world origin. If the `Group ID` input is used this offset is applied on a per-group basis. If the selection is used, only the selected points contribute towards the calculation of the centroid, but all points are still moved",
                    videos="https://imgur.com/xSOH4Tr",
                ),
                MenuItem(
                    name="Separate First Point",
                    description="Separate out the first point for each group in the `Group ID`",
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
                    name="Index Mixed",
                    description="Offset the current point's `Index` by a given amount, returning the fractional Index value",
                ),
                MenuItem(
                    name="Index Mix Float",
                    description="Sample and interpolate the `Float` value between the floor and ceiling of the given `Index`, by the fractional amount",
                ),
                MenuItem(
                    name="Index Mix Vector",
                    description="Sample and interpolate the `Vector` value between the floor and ceiling of the given `Index`, by the fractional amount",
                ),
                MenuItem(
                    name="Index Mix Rotation",
                    description="Sample and interpolate the `Rotation` value between the floor and ceiling of the given `Index`, by the fractional amount",
                ),
                MenuItem(
                    name="Index Mix Color",
                    description="Sample and interpolate the `Color` value between the floor and ceiling of the given `Index`, by the fractional amount",
                ),
                Break(),
                MenuItem(
                    name="Sample Position",
                    description="Sample the `Position` attribute from a point at the given `Index`",
                ),
                MenuItem(
                    name="Sample Mixed Float",
                    description="Sample the `Float` attribute, mixed between the floor and ceiling of the given `Index`",
                ),
                MenuItem(
                    name="Sample Mixed Vector",
                    description="Sample the `Vector` attribute, mixed between the floor and ceiling of the given `Index`",
                ),
                MenuItem(
                    name="Sample Mixed Rotation",
                    description="Sample the `Rotation` attribute, mixed between the floor and ceiling of the given `Index`",
                ),
                MenuItem(
                    name="Sample Mixed Color",
                    description="Sample the `Color` attribute, mixed between the floor and ceiling of the given `Index`",
                ),
                Break(),
                MenuItem(
                    name="Offset Index",
                    description="Add an integer offset to the point's `Index`",
                ),
                MenuItem(
                    name="Offset Integer",
                    description="Evaluate an `Integer` at an index that is offset by the specified amount",
                ),
                MenuItem(
                    name="Offset Float",
                    description="Evaluate a `Float` value at an index that is offset by the specified amount",
                ),
                MenuItem(
                    name="Offset Vector",
                    description="Evaluate a `Vector` at an index that is offset by the specified amount",
                ),
                MenuItem(
                    name="Offset Boolean",
                    description="Evaluate a `Boolean` at an index that is offset by the specified amount",
                ),
                MenuItem(
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
                    name="Group Pick",
                    description="For each group, return the index of the point for which the Selection is true. Only valid if there is a single true in the group. If not lvalid, returns -1",
                ),
                MenuItem(
                    name="Group Pick First",
                    description="Similar to `Group Pick`, but will always return `Index` of first true element, even if there are multiple true values. If nothing is true then will return `-1`",
                ),
                MenuItem(
                    name="Group Pick Vector",
                    description="For each group, return the Position of the point at which the selection is true. If there is more than one true for the group the pick is not valid and (0, 0, 0) is returned",
                ),
                MenuItem(
                    name="Relative Index",
                    description="Based on the Group ID input, return the size of the group and the indices of the first and last items of the group",
                ),
                Break(),
                MenuItem(
                    name="Attribute Run",
                    description="Fill in gaps in a set of continuous boolean True values, up to a specific size",
                ),
                MenuItem(
                    name="Boolean First",
                    description="For each `Group ID`, every value becomes `False` except the first `True` value",
                ),
                MenuItem(
                    name="Integer Run",
                    description="Group mask output increments whenever the input value or the Group ID changes",
                ),
                MenuItem(
                    name="Boolean Run Fill",
                    description="Fill in gaps in a set of continuous boolean True values, up to a specific size",
                ),
                MenuItem(
                    name="Boolean Run Trim",
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
                    videos="https://imgur.com/TNc102v",
                ),
                MenuItem(
                    label="Center Assembly",
                    name="MN_assembly_center",
                    description="Move an instanced assembly to the world origin. Some structures are not centred on the world origin, so this node can reset them to the world origin for convenient rotation and translation and animation.",
                    videos="https://imgur.com/pgFTmgC",
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
                    videos="https://imgur.com/m3BPUxh",
                ),
                MenuItem(
                    label="Animate Value",
                    name="Animate Value",
                    description="Animate a float value between the specified min and max values, over specified range of frames. If clamped, frames above and below the start and end will result in the min and max output values, otherwise it will continue to linearly interpolate the value beyond the min and max values.",
                    videos="https://imgur.com/2oOnwRm",
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
                    videos="https://imgur.com/GK1nyUz",
                ),
                MenuItem(
                    label="Peptide to Curve",
                    name="Animate Peptide to Curve",
                    description="Take the protein residues from a structure and align then along an input curve. Editing the curve will change how the atoms are arranged. The output atoms can be styled as normal.",
                    videos="https://imgur.com/FcEXSZx",
                ),
                Break(),
                MenuItem(
                    label="Noise Position",
                    name="MN_animate_noise_position",
                    description="Create 3D noise vector based on the position of points in 3D space. Evolve the noise function with the `Animate` input, and change the characteristics of the noise function with the other inputs such as scale and detail. There is also a 1-dimensional noise output called `Fac`.\n\nAn example of using this noise is to offset the positions of atoms with the `Set Position` node.",
                    videos="https://imgur.com/B8frW1C",
                ),
                MenuItem(
                    label="Noise Field",
                    name="MN_animate_noise_field",
                    description="Create a 3D noise vector based on the input field. Evolve the noise function with the `Animate` input, and change the characteristics of the noise function with the other inputs such as scale and detail. There is also a 1-dimensional noise output called `Fac`.\n\nAn example of using this noise is to offset the positions of atoms with the `Set Position` node. Different field inputs result in different noise being applied. Using the `chain_id` results in the same noise being generated for each atom in each chain, but different between chains.",
                    videos="https://imgur.com/hqemVQy",
                ),
                MenuItem(
                    label="Noise Repeat",
                    name="MN_animate_noise_repeat",
                    description="Create a 3D noise vector based on the input field, that repeats every `1.00` for the `Animate 0..1` input. Evolve the noise function with the `Animate` input, and change the characteristics of the noise function with the other inputs such as scale and detail. There is also a 1-dimensional noise output called `Fac`.\n\nAn example of using this noise is to offset the positions of atoms with the `Set Position` node. Different field inputs result in different noise being applied. Using the `chain_id` results in the same noise being generated for each atom in each chain, but different between chains.",
                    videos="https://imgur.com/GNQcIlx",
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
                MenuItem(name="Angstrom to World"),
                MenuItem(name="Nanometre to World"),
                MenuItem(name="World to Angstrom"),
                MenuItem(name="World to Nanometre"),
                Break(),
                MenuItem(
                    name="Between Integer",
                    description="Test if an integer is between (and including) the upper and lower bounds",
                ),
                MenuItem(
                    name="Between Float",
                    description="Test if a float is between the upper and lower bounds",
                ),
                MenuItem(
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
                    description="Mix the current point's `Position` with the input vector. A convenience wrapper around a `Mix Vector` using `Position` as the first attribute",
                    videos="https://imgur.com/Cc538lr",
                ),
                Break(),
                MenuItem(
                    name="Fractionate Float",
                    description="Test if a vector is element-wise between the upper and lower bounds.",
                ),
                Break(),
                MenuItem(
                    name="Rotation from ZYZ",
                    description="Combine a rotation defined as `ZYZ`,common in electron tomography",
                ),
                Break(),
                MenuItem(
                    name="Transform Scale",
                    description="Scale the components of a transform individually with between 0 and their value",
                ),
                MenuItem(
                    name="Transform Relative",
                    description="The transform to get from B to A, relative to the CB axis",
                ),
                MenuItem(
                    name="Transform Mix",
                    description="Mix between two transforms, controlling the translation, rotation and scale independently",
                ),
                MenuItem(
                    name="Transform Local",
                    description="",
                ),
                MenuItem(
                    name="Transform Local Axis",
                    description="Create a transform that is rotation around an axis in local space, with local space defined by the origin point which defaults to Position",
                ),
                MenuItem(
                    name="Transform Accumulate",
                    description="Accumulate transforms on the domain if selected",
                ),
                MenuItem(
                    name="Transform Accumulate Point",
                    description="Accumulate transforms on a domain, applying these transforms to the Position",
                ),
                Break(),
                MenuItem(
                    name="Centroid",
                    description="Calculate the centroid point for the selection for each group in the `Group ID`. The centroid is the average position of all points in each `Group ID`. If a selection is given, then only the selected points contribute towards the overall centroid calculation, but the result is still available on the other points in the `Group ID`",
                    videos="https://imgur.com/Cc538lr",
                ),
                MenuItem(
                    name="Vector Direction",
                    description="Direction from one point to another, optionally normalized",
                ),
                MenuItem(
                    name="Vector Angle",
                    description="Compute the angle in radians between two vectors.",
                ),
                MenuItem(
                    name="Dihedral Angle",
                    description='Computes the angle between two vectors, AB & CD around around the axis of BC. The first vector AB is treated as the "12 O\'clock" up position, looking down the axis towards C, with angles being return in the range of (-Pi, Pi). Clockwise angles are positive and anti-clockwise angles are negative.',
                    videos="",
                ),
                MenuItem(
                    name="3 Point Angle",
                    description="Calculate the angle between 3 different points. These points are selected based on their index in the point domain, with Index B being the centre of the calculation.\n\nIn the video example, the same calculation that is occurring internally inside of the `MN_topo_edge_angle` node, is being handled explicity by this node. If the `Index` is being used as `Index B` then the current point that is being evaluated is the centre of the angle calculation. If this value is changed, then the point at the corresponding index is used, which results in a smaller angle in the example video.",
                    videos="https://imgur.com/qXyy2ln",
                ),
                MenuItem(
                    name="2 Point Angle",
                    description="Calculate the angle that two points make, relative to the current point being evaluated. Points are selected based on their index, with the centre of the angle calculation being the current point's position. Equivalent to using 3-Point angle and using `Index` as the `Index B`.\n\nIn the example video, the angle calculation is similar to that of the 3-Point Angle node, but the middle point is always the current point.",
                    videos="https://imgur.com/xp7Vbaj",
                ),
                MenuItem(
                    name="Point Distance",
                    description="Calculate the distance and the vector between the evaluating point and the point selected via the Index.\n\nIn the example video, each point is calculating a vector and a distance between itself and the indexed point. When the Point Mask node is used, this index is then on a per-group basis, so each point in the group points to just the group's corresponding point.",
                    videos="https://imgur.com/AykNvDz",
                ),
                # MenuItem(
                #     label="Cartoon Utilities",
                #     name=".MN_utils_style_cartoon",
                #     description="The underlying node group which powers the cartoon style",
                # ),
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
                    videos="https://imgur.com/jGgMSd4",
                ),
                MenuItem(
                    label="Style Wire",
                    name="Style Density Wire",
                    description="A wire surface made from the electron density given a certain threshold value.",
                    videos="https://imgur.com/jGgMSd4",
                ),
                MenuItem(
                    label="Sample Nearest Attribute",
                    name="Sample Nearest Atoms",
                    description="Sample the nearest atoms from another object, to get the colors or other attributes and apply them to a volume mesh.",
                    videos="https://imgur.com/UzNwLv2",
                ),
            ],
        ),
    ]
)
