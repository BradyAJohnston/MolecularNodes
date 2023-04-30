import bpy
import numpy as np
from . import data
from . import coll
from .load import create_object, add_attribute
import warnings

class TrajectorySelectionList(bpy.types.PropertyGroup):
    """Group of properties for custom selections for MDAnalysis import."""
    
    name: bpy.props.StringProperty(
        name="Attribute Name", 
        description="Attribute", 
        default="custom_selection"
    )
    
    selection: bpy.props.StringProperty(
        name="Selection String", 
        description="String that provides a selection through MDAnalysis", 
        default = "name CA"
    )

class MOL_UL_TrajectorySelectionListUI(bpy.types.UIList):
    """UI List"""
    
    def draw_item(self, context, layout, data, item, 
                  icon, active_data, active_propname, index):
        custom_icon = "VIS_SEL_11"
        
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            layout.label(text = item.name, icon = custom_icon)
        
        elif self.layout_type in {'GRID'}:
            layout.alignment = 'CENTER'
            layout.label(text = "", icon = custom_icon)
            

class TrajectorySelection_OT_NewItem(bpy.types.Operator):
    """Add a new custom selection to the list."""
    
    bl_idname = "trajectory_selection_list.new_item"
    bl_label = "+"
    
    def execute(self, context):
        context.scene.trajectory_selection_list.add()
        return {'FINISHED'}

class TrajectorySelection_OT_DeleteIem(bpy.types.Operator):
    
    bl_idname = "trajectory_selection_list.delete_item"
    bl_label = "-"
    
    @classmethod
    def poll(cls, context):
        return context.scene.trajectory_selection_list
    def execute(self, context):
        my_list = context.scene.trajectory_selection_list
        index = context.scene.list_index
        
        my_list.remove(index)
        context.scene.list_index = min(max(0, index - 1), len(my_list) - 1)
        
        return {'FINISHED'}

def load_trajectory(file_top, 
                    file_traj,
                    md_start = 1, 
                    md_end = 50, 
                    md_step = 1,
                    world_scale = 0.01, 
                    include_bonds = False, 
                    del_solvent = False,
                    selection = "not (name H* or name OW)",
                    name = "default",
                    custom_selections = None,
                    ):
    
    import MDAnalysis as mda
    import MDAnalysis.transformations as trans
    
    # initially load in the trajectory
    if file_traj == "":
        univ = mda.Universe(file_top)
    else:
        univ = mda.Universe(file_top, file_traj)
        
    # separate the trajectory, separate to the topology or the subsequence selections
    traj = univ.trajectory[md_start:md_end:md_step]
    
    # if there is a non-blank selection, apply the selection text to the universe for 
    # later use. This also affects the trajectory, even though it has been separated earlier
    if selection != "":
        try:
            univ = univ.select_atoms(selection)
        except:
            warnings.warn(f"Unable to apply selection: '{selection}'. Loading entire topology.")
    
    # Try and extract the elements from the topology. If the universe doesn't contain
    # the element information, then guess based on the atom names in the toplogy
    try:
        elements = univ.atoms.elements.tolist()
    except:
        try:
            elements = [mda.topology.guessers.guess_atom_element(x) for x in univ.atoms.names]
        except:
            pass
        
    
    
    if hasattr(univ, 'bonds') and include_bonds:

            # If there is a selection, we need to recalculate the bond indices
            if selection != "":
                index_map = { index:i for i, index in enumerate(univ.atoms.indices) }

                new_bonds = []
                for bond in univ.bonds.indices:
                    try:
                        new_index = [index_map[y] for y in bond]
                        new_bonds.append(new_index)
                    except KeyError:
                        # fragment - one of the atoms in the bonds was 
                        # deleted by the selection, so we shouldn't 
                        # pass this as a bond.  
                        pass
                    
                bonds = np.array(new_bonds)
            else:
                bonds = univ.bonds.indices

    else:
        bonds = []

    
    # create the initial model
    mol_object = create_object(
        name = name,
        collection = coll.mn(),
        locations = univ.atoms.positions * world_scale, 
        bonds = bonds
    )
    
    ## add the attributes for the model
    
    # The attributes for the model are initially defined as single-use functions. This allows
    # for a loop that attempts to add each attibute by calling the function. Only during this
    # loop will the call fail if the attribute isn't accessible, and the warning is reported
    # there rather than setting up a try: except: for each individual attribute which makes
    # some really messy code.
    
    def att_atomic_number():
        atomic_number = np.array(list(map(
            # if getting the element fails for some reason, return an atomic number of -1
            lambda x: data.elements.get(x, {"atomic_number": -1}).get("atomic_number"), 
            np.char.title(elements) 
        )))
        return atomic_number

    def att_vdw_radii():
        try:
            vdw_radii = np.array(list(map(
                lambda x: mda.topology.tables.vdwradii.get(x, 1), 
                np.char.upper(elements)
            )))
        except:
            # if fail to get radii, just return radii of 1 for everything as a backup
            vdw_radii = np.ones(len(univ.atoms.names))
            warnings.warn("Unable to extract VDW Radii. Defaulting to 1 for all points.")
        
        return vdw_radii * world_scale
    
    def att_res_id():
        return univ.atoms.resnums
    
    def att_res_name():
        res_names =  np.array(list(map(lambda x: x[0: 3], univ.atoms.resnames)))
        res_numbers = np.array(list(map(
            lambda x: data.residues.get(x, {'res_name_num': 0}).get('res_name_num'), 
            res_names
            )))
        return res_numbers
    
    def att_b_factor():
        return univ.atoms.tempfactors
    
    def att_chain_id():
        chain_id = univ.atoms.chainIDs
        chain_id_unique = np.unique(chain_id)
        chain_id_num = np.array(list(map(lambda x: np.where(x == chain_id_unique)[0][0], chain_id)))
        mol_object['chain_id_unique'] = chain_id_unique
        return chain_id_num
    
    # returns a numpy array of booleans for each atom, whether or not they are in that selection
    def bool_selection(selection):
        return np.isin(univ.atoms.ix, univ.select_atoms(selection).ix).astype(bool)
    
    def att_is_backbone():
        return bool_selection("backbone or nucleicbackbone or name BB")
    
    def att_is_alpha_carbon():
        return bool_selection('name CA or name BB')
    
    def att_is_lipid():
        return bool_selection('resname 23SM ABLIPA ABLIPB ADR ADRP ALIN ALINP APC APPC ARA ARAN ARANP ARAP ASM BCLIPA BCLIPB BCLIPC BEH BEHP BNSM BSM C6DHPC C7DHPC CER160 CER180 CER181 CER2 CER200 CER220 CER240 CER241 CER3E CHAPS CHAPSO CHL1 CHM1 CHNS CHOA CHSD CHSP CJLIPA CPC CTLIPA CYFOS3 CYFOS4 CYFOS5 CYFOS6 CYFOS7 CYSF CYSG CYSL CYSP DAPA DAPA DAPC DAPC DAPE DAPE DAPG DAPG DAPS DAPS DBPA DBPC DBPE DBPG DBPS DBSM DCPC DDA DDAO DDAOP DDAP DDMG DDOPC DDOPE DDOPS DDPC DEPA DEPC DEPE DEPG DEPS DFPA DFPC DFPE DFPG DFPS DGLA DGLAP DGPA DGPA DGPC DGPC DGPE DGPE DGPG DGPG DGPS DGPS DHA DHAP DHPC DHPCE DIPA DIPA DIPC DIPE DIPG DIPS DLIPC DLIPE DLIPI DLPA DLPA DLPC DLPC DLPE DLPE DLPG DLPG DLPS DLPS DMPA DMPC DMPCE DMPE DMPEE DMPG DMPI DMPI13 DMPI14 DMPI15 DMPI24 DMPI25 DMPI2A DMPI2B DMPI2C DMPI2D DMPI33 DMPI34 DMPI35 DMPS DNPA DNPA DNPC DNPC DNPE DNPE DNPG DNPG DNPS DNPS DOMG DOPA DOPA DOPC DOPC DOPCE DOPE DOPE DOPEE DOPG DOPG DOPP1 DOPP2 DOPP3 DOPS DOPS DPA DPAP DPC DPCE DPP1 DPP2 DPPA DPPA DPPC DPPC DPPE DPPE DPPEE DPPG DPPG DPPGK DPPI DPPS DPPS DPSM DPT DPTP DRPA DRPC DRPE DRPG DRPS DSPA DSPC DSPE DSPG DSPS DTPA DTPA DTPC DTPE DTPG DTPS DUPC DUPE DUPS DVPA DVPC DVPE DVPG DVPS DXCE DXPA DXPA DXPC DXPC DXPE DXPE DXPG DXPG DXPS DXPS DXSM DYPA DYPA DYPC DYPC DYPE DYPE DYPG DYPG DYPS DYPS ECLIPA ECLIPB ECLIPC EDA EDAP EICO EICOP EPA EPAP ERG ERU ERUP ETA ETAP ETE ETEP FOIS11 FOIS9 FOS10 FOS12 FOS13 FOS14 FOS15 FOS16 GLA GLAP GLYM HPA HPAP HPLIPA HPLIPB HTA HTAP IPC IPPC KPLIPA KPLIPB KPLIPC LAPAO LAPAOP LAU LAUP LDAO LDAOP LIGN LIGNP LILIPA LIN LINP LLPA LLPC LLPE LLPS LMPG LNACL1 LNACL2 LNBCL1 LNBCL2 LNCCL1 LNCCL2 LNDCL1 LNDCL2 LOACL1 LOACL2 LOCCL1 LOCCL2 LPC LPC12 LPC14 LPPA LPPC LPPC LPPE LPPG LPPG LPPS LSM LYSM MCLIPA MEA MEAP MYR MYRO MYROP MYRP NER NERP NGLIPA NGLIPB NGLIPC NSM OLE OLEP OPC OSM OSPE OYPE PADG PAL PALIPA PALIPB PALIPC PALIPD PALIPE PALO PALOP PALP PAPA PAPC PAPE PAPG PAPI PAPS PDOPC PDOPE PEPC PGPA PGPC PGPE PGPG PGPS PGSM PIDG PIM1 PIM2 PIPA PIPC PIPE PIPG PIPI PIPS PLPA PLPC PLPE PLPG PLPI PLPI13 PLPI14 PLPI15 PLPI24 PLPI25 PLPI2A PLPI2B PLPI2C PLPI2D PLPI33 PLPI34 PLPI35 PLPS PMCL1 PMCL2 PMPE PMPG PNCE PNPI PNPI13 PNPI14 PNPI15 PNPI24 PNPI25 PNPI2A PNPI2B PNPI2C PNPI2D PNPI33 PNPI34 PNPI35 PNSM PODG POP1 POP2 POP3 POPA POPA POPC POPC POPCE POPE POPE POPEE POPG POPG POPI POPI POPI13 POPI14 POPI15 POPI24 POPI25 POPI2A POPI2B POPI2C POPI2D POPI33 POPI34 POPI35 POPP1 POPP2 POPP3 POPS POPS POSM PPC PPPE PQPE PQPS PRPA PRPC PRPE PRPG PRPS PSM PSPG PUDG PUPA PUPC PUPE PUPI PUPS PVCL2 PVDG PVP1 PVP2 PVP3 PVPE PVPG PVPI PVSM PYPE PYPG PYPI PhPC QMPE SAPA SAPC SAPE SAPG SAPI SAPI13 SAPI14 SAPI15 SAPI24 SAPI25 SAPI2A SAPI2B SAPI2C SAPI2D SAPI33 SAPI34 SAPI35 SAPS SB3-10 SB3-12 SB3-14 SDA SDAP SDPA SDPC SDPE SDPG SDPS SDS SELIPA SELIPB SELIPC SFLIPA SITO SLPA SLPC SLPE SLPG SLPS SOPA SOPC SOPE SOPG SOPS SSM STE STEP STIG THA THAP THCHL THDPPC TIPA TLCL1 TLCL2 TMCL1 TMCL2 TOCL1 TOCL2 TPA TPAP TPC TPC TPT TPTP TRI TRIP TRIPAO TRPAOP TSPC TTA TTAP TXCL1 TXCL2 TYCL1 TYCL2 UDAO UDAOP UFOS10 UPC VCLIPA VCLIPB VCLIPC VCLIPD VCLIPE VPC XNCE XNSM YOPA YOPC YOPE YOPS YPLIPA YPLIPB')
    
    def att_is_solvent():
        return bool_selection('name OW or name HW1 or name HW2 or resname W')
    
    def att_atom_type():
        return np.array(univ.atoms.types, dtype = int)
    
    def att_is_nucleic():
        return bool_selection('nucleic')
    
    def att_is_peptide():
        return bool_selection('protein')

    attributes = (
        {'name': 'atomic_number',   'value': att_atomic_number,   'type': 'INT',     'domain': 'POINT'}, 
        {'name': 'vdw_radii',       'value': att_vdw_radii,       'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'res_id',          'value': att_res_id,          'type': 'INT',     'domain': 'POINT'}, 
        {'name': 'res_name',        'value': att_res_name,        'type': 'INT',     'domain': 'POINT'}, 
        {'name': 'b_factor',        'value': att_b_factor,        'type': 'float',   'domain': 'POINT'}, 
        {'name': 'chain_id',        'value': att_chain_id,        'type': 'INT',     'domain': 'POINT'}, 
        {'name': 'atom_types',      'value': att_atom_type,       'type': 'INT',     'domain': 'POINT'}, 
        {'name': 'is_backbone',     'value': att_is_backbone,     'type': 'BOOLEAN', 'domain': 'POINT'}, 
        {'name': 'is_alpha_carbon', 'value': att_is_alpha_carbon, 'type': 'BOOLEAN', 'domain': 'POINT'}, 
        {'name': 'is_solvent',      'value': att_is_solvent,      'type': 'BOOLEAN', 'domain': 'POINT'}, 
        {'name': 'is_nucleic',      'value': att_is_nucleic,      'type': 'BOOLEAN', 'domain': 'POINT'}, 
        {'name': 'is_peptide',      'value': att_is_peptide,      'type': 'BOOLEAN', 'domain': 'POINT'}, 
        {'name': 'is_lipid',        'value': att_is_lipid,        'type': 'BOOLEAN', 'domain': 'POINT'}, 
    )
    
    for att in attributes:
        # tries to add the attribute to the mesh by calling the 'value' function which returns
        # the required values do be added to the domain.
        try:
            add_attribute(mol_object, att['name'], att['value'](), att['type'], att['domain'])
        except:
            warnings.warn(f"Unable to add attribute: {att['name']}.")

    # add the custom selections if they exist
    if custom_selections:
        for sel in custom_selections:
            try:
                add_attribute(
                    object=mol_object, 
                    name=sel.name, 
                    data=bool_selection(sel.selection), 
                    type = "BOOLEAN", 
                    domain = "POINT"
                    )
            except:
                warnings.warn("Unable to add custom selection: {}".format(sel.name))

    coll_frames = coll.frames(name)
    
    add_occupancy = True
    for ts in traj:
        frame = create_object(
            name = name + "_frame_" + str(ts.frame),
            collection = coll_frames, 
            locations = univ.atoms.positions * world_scale
        )
        # adds occupancy data to each frame if it exists
        # This is mostly for people who want to store frame-specific information in the 
        # b_factor but currently neither biotite nor MDAnalysis give access to frame-specific
        # b_factor information. MDAnalysis gives frame-specific access to the `occupancy` 
        # so currently this is the only method to get frame-specific data into MN
        # for more details: https://github.com/BradyAJohnston/MolecularNodes/issues/128
        if add_occupancy:
            try:
                add_attribute(frame, 'occupancy', ts.data['occupancy'])
            except:
                add_occupancy = False
    
    # disable the frames collection from the viewer
    bpy.context.view_layer.layer_collection.children[coll.mn().name].children[coll_frames.name].exclude = True
    
    return mol_object, coll_frames
    
