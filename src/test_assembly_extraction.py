import bpy

assembly_id = 1
## From the outpudt of atomium.fetch("PDB").assemblies, use the assembly information to 
# generate a properties model that encodes the transformations matrices into the 
# vertex positions of the model

def assembly_to_vec_list(assembly):
    """
    Combines all of the transformation matrices into a list of triplet vectors, that
    can be used to create a mesh of vertices for including the transformation properties.
    """
    vec_list = []

    for mat in assembly.get("transformations"):
        # append the 3 rows of the rotation matrix
        for row in mat.get("matrix"):
            vec_list.append(row)
        # append the translation vector of thr translation matrix
        vec_list.append(mat.get("vector"))

    # return the list of vectors that will be used to create 
    return(vec_list)





coll_properties = bpy.data.collections.get(output_name + "_properties")
coll_mol = bpy.data.collections.get(output_name)

coll_assemblies = bpy.data.collections.get(output_name + "_assemblies")

if not coll_assemblies:
    coll_assemblies = bpy.data.collections.new(output_name + "_assemblies")
    coll_mol.children.link(coll_assemblies)
    
    # hide the assemblies collection

# bpy.context.view_layer.active_layer_collection.exclude = True
try:
    bpy.context.layer_collection.children[coll_assemblies.name].exclude = True
except:
    # warning("")
    pass


# col = bpy.data.collections.new(output_name)
# parent_coll.children.link(col)

def create_model(name, collection, locations, bonds=[], faces=[]):
    """
    Creates a mesh with the given name in the given collection, from the supplied
    values for the locationso of vertices, and if supplied, bonds and faces.
    """
    # create a new mesh
    atom_mesh = bpy.data.meshes.new(name)
    atom_mesh.from_pydata(locations, bonds, faces)
    new_object = bpy.data.objects.new(name, atom_mesh)
    collection.objects.link(new_object)
    return new_object

coll_assemblies_name = output_name + "_assemblies"


counter = 0
for assembly in assemblies:
    counter += 1

    new_prop_name = output_name + "_properties_assembly_" + str(counter)
    properties_vertices = assembly_to_vec_list(assembly)

    new_prop = coll_assemblies.all_objects.get(new_prop_name)

    if not new_prop:
        new_prop = create_model(
            name = new_prop_name, 
            collection = coll_assemblies, 
            locations = properties_vertices
            )
        
        if counter == 1:
            first_assembly = new_prop
