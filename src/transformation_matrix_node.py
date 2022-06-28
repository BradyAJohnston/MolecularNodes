

vertex_list = list()

for mat in transformation:
    mat_rot = mat.get("matrix")
    mat_rot_trans = [mat_rot[0], mat_rot[1], mat_rot[2], mat.get("vector")]
    for vec in mat_rot_trans:
        vertex_list.append(vec)