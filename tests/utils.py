import bpy
import molecularnodes as mn
import numpy as np
import random

def apply_mods(obj):
    """
    Applies the modifiers on the modifier stack
    
    This will realise the computations inside of any Geometry Nodes modifiers, ensuring
    that the result of the node trees can be compared by looking at the resulting 
    vertices of the object.
    """
    bpy.context.view_layer.objects.active = obj
    for modifier in obj.modifiers:
        bpy.ops.object.modifier_apply(modifier = modifier.name)


def sample_attribute(object,
                    attribute,
                    n = 100,
                    seed = 6):
    random.seed(seed)
    attribute = mn.obj.get_attribute(object, attribute)
    length = len(attribute)
    
    if n > length:
        idx = range(length)
    else:
        idx = random.sample(range(length), n)
    
    dimensions = len(np.shape(attribute))
    
    if dimensions == 1:
        array = attribute[idx]
    elif dimensions == 2:
        array = attribute[idx, :]
    else:
        Warning("Unable to sample higher dimensional attribute")
    
    return array


def sample_attribute_to_string(object,
                               attribute,
                               n = 100,
                               precision=3,
                               seed = 6):
    array = sample_attribute(object=object, attribute=attribute, n=n, seed=seed)
    length = len(array)
    threshold = 4 * length
    
    if n > length:
        idx = range(length)
    else:
        idx = random.sample(range(length), n)
    
    dimensions = len(np.shape(attribute))
    
    if dimensions == 1:
        array = attribute[idx]
    elif dimensions == 2:
        array = attribute[idx, :]
    else:
        Warning("Unable to sample higher dimensional attribute")
    
    return np.array2string(array, precision=precision, threshold=threshold)

def get_verts(obj, float_decimals=4, n_verts=100, apply_modifiers=True, seed=42):
    """
    Randomly samples a specified number of vertices from an object.

    Parameters
    ----------
    obj : object
        Object from which to sample vertices.
    float_decimals : int, optional
        Number of decimal places to round the vertex coordinates, defaults to 4.
    n_verts : int, optional
        Number of vertices to sample, defaults to 100.
    apply_modifiers : bool, optional
        Whether to apply all modifiers on the object before sampling vertices, defaults to True.
    seed : int, optional
        Seed for the random number generator, defaults to 42.

    Returns
    -------
    str
        String representation of the randomly selected vertices.

    Notes
    -----
    This function randomly samples a specified number of vertices from the given object.
    By default, it applies all modifiers on the object before sampling vertices. The
    random seed can be set externally for reproducibility.

    If the number of vertices to sample (`n_verts`) exceeds the number of vertices
    available in the object, all available vertices will be sampled.

    The vertex coordinates are rounded to the specified number of decimal places
    (`float_decimals`) before being included in the output string.

    Examples
    --------
    >>> obj = mn.load.molecule_rcsb('6n2y', starting_style=2)
    >>> get_verts(obj, float_decimals=3, n_verts=50, apply_modifiers=True, seed=42)
    '1.234,2.345,3.456\n4.567,5.678,6.789\n...'
    """

    import random

    random.seed(seed)

    
    if apply_modifiers:
        try:
            apply_mods(obj)
        except RuntimeError as ex:
            return str(ex)

    vert_list = [(v.co.x, v.co.y, v.co.z) for v in obj.data.vertices]

    if n_verts > len(vert_list):
        n_verts = len(vert_list)

    random_verts = random.sample(vert_list, n_verts)

    verts_string = ""
    for i, vert in enumerate(random_verts):
        if i < n_verts:
            rounded = [round(x, float_decimals) for x in vert]
            verts_string += "{},{},{}\n".format(rounded[0], rounded[1], rounded[2])

    return verts_string


def remove_all_molecule_objects(mda_session):
    for object in bpy.data.objects:
        try:
            obj_type = object["type"]
            if obj_type == "molecule":
                bpy.data.objects.remove(object)
        except KeyError:
            pass
    # remove frame change
    bpy.context.scene.frame_set(0)

    mda_session.universe_reps = {}
    mda_session.atom_reps = {}
    mda_session.rep_names = []

