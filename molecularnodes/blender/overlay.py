# intially taken from here: https://gist.github.com/PardhavMaradani/13d77d8c107e3d3c17228dc8598b11bd

import bpy
import blf  # type: ignore
import gpu  # type: ignore
from gpu_extras.batch import batch_for_shader  # type: ignore
from bpy_extras import view3d_utils  # type: ignore
from math import fabs, radians, sqrt, cos, sin

# External


def show_annotations():
    global _draw_handler
    _draw_handler = bpy.types.SpaceView3D.draw_handler_add(
        _draw_annotations_handler, (None, bpy.context), "WINDOW", "POST_PIXEL"
    )
    _redraw(bpy.context)


def hide_annotations():
    global _draw_handler
    if _draw_handler is not None:
        bpy.types.SpaceView3D.draw_handler_remove(_draw_handler, "WINDOW")
        _redraw(bpy.context)
        _draw_handler = None


def clear_annotations():
    global _annotations
    _annotations = []
    _redraw(bpy.context)


def print_annotations():
    print(_annotations)


def show_atom_info(obj_name, selection, options=None):
    if not _valid_object(obj_name):
        print("Invalid Object", obj_name)
        return
    selection_array = _make_selection_array(selection)
    ag = _get_selection(obj_name, selection_array)
    if ag is None:
        print("Invalid selection", obj_name, selection)
        return
    print(repr(ag))
    annotation = {
        "type": "atom_info",
        "obj_name": obj_name,
        "selection": selection_array,
        "options": options,
    }
    _annotations.append(annotation)
    _redraw(bpy.context)


def show_bond_angle(obj_name, selection, options=None):
    if not _valid_object(obj_name):
        print("Invalid Object", obj_name)
        return
    selection_array = _make_selection_array(selection)
    ag = _get_selection(obj_name, selection_array)
    if ag is None:
        print("Invalid selection", obj_name, selection)
        return
    if ag.n_atoms != 3:
        print("Need 3 atoms, selection returned", ag.n_atoms)
        return
    print(repr(ag))
    annotation = {
        "type": "bond_angle",
        "obj_name": obj_name,
        "selection": selection_array,
        "options": options,
    }
    _annotations.append(annotation)
    _redraw(bpy.context)


def show_com(obj_name, selection, options=None):
    if not _valid_object(obj_name):
        print("Invalid Object", obj_name)
        return
    selection_array = _make_selection_array(selection)
    ag = _get_selection(obj_name, selection_array)
    if ag is None:
        print("Invalid selection", obj_name, selection)
        return
    print(repr(ag))
    annotation = {
        "type": "com",
        "obj_name": obj_name,
        "selection": selection_array,
        "options": options,
    }
    _annotations.append(annotation)
    _redraw(bpy.context)


def show_com_distance(obj_name, selection1, selection2, options=None):
    if not _valid_object(obj_name):
        print("Invalid Object", obj_name)
        return
    selection1_array = _make_selection_array(selection1)
    ag1 = _get_selection(obj_name, selection1_array)
    if ag1 is None:
        print("Invalid selection", obj_name, selection1)
        return
    selection2_array = _make_selection_array(selection2)
    ag2 = _get_selection(obj_name, selection2_array)
    if ag1 is None:
        print("Invalid selection", obj_name, selection2)
        return
    print(repr(ag1), repr(ag2))
    annotation = {
        "type": "com_distance",
        "obj_name": obj_name,
        "selection1": selection1_array,
        "selection2": selection2_array,
        "options": options,
    }
    _annotations.append(annotation)
    _redraw(bpy.context)


def show_canonical_dihedrals(obj_name, resid, options=None):
    if not _valid_object(obj_name):
        print("Invalid Object", obj_name)
        return
    u = _get_universe(obj_name)
    try:
        r = u.residues[resid - 1]
        print(repr(r))
        phi = r.phi_selection()
        print("phi", phi.names, phi.dihedral.value())
        psi = r.psi_selection()
        print("psi", psi.names, psi.dihedral.value())
        omega = r.omega_selection()
        print("omega", omega.names, omega.dihedral.value())
        chi1 = r.chi1_selection()
        print("chi1", chi1.names, chi1.dihedral.value())
    except Exception as e:
        print("Exception :", e)
        return
    annotation = {
        "type": "dihedrals",
        "obj_name": obj_name,
        "resid": resid,
        "options": options,
    }
    _annotations.append(annotation)
    _redraw(bpy.context)


def get_universe(obj_name):
    return _get_universe(obj_name)


# Internal

_fmt = "%1.2f"
_annotations = []
_viewport = (0, 0)
_draw_handler = None
_rad45 = radians(45)
_rad315 = radians(315)
_bsf = 0.01
_default_color = (1, 1, 1, 1)
_default_font_size = 16
_default_line_width = 1.0
_default_arrow_size = 16
_default_text_align = "C"
_default_text_rotation = 0.0
_shader_line = (
    gpu.shader.from_builtin("POLYLINE_UNIFORM_COLOR")
    if not bpy.app.background
    else None
)


def _set_viewport():
    global _viewport
    _viewport = (bpy.context.region.width, bpy.context.region.height)


def _valid_object(obj_name):
    if obj_name in bpy.context.scene.objects and hasattr(
        bpy.context.scene.objects[obj_name], "uuid"
    ):
        uuid = bpy.context.scene.objects[obj_name].uuid
        if (
            hasattr(bpy.context.scene, "MNSession")
            and uuid in bpy.context.scene.MNSession.trajectories
        ):
            return True
    return False


def _make_selection_array(selection):
    selection_array = selection
    if type(selection) is str:
        selection_array = [selection]
    return selection_array


def _get_selection(obj_name, selection):
    u = _get_universe(obj_name)
    try:
        ag = u.select_atoms(*selection)
        return ag
    except Exception as e:
        print("Exception :", e)
    return None


def _redraw(context):
    if _draw_handler is not None:
        area = next(area for area in context.screen.areas if area.type == "VIEW_3D")
        area.tag_redraw()


def _draw_annotations_handler(self, context):
    _set_viewport()
    region = bpy.context.region
    if not context.space_data.region_quadviews:
        rv3d = bpy.context.space_data.region_3d
    else:
        if context.area.type != "VIEW_3D" or context.space_data.type != "VIEW_3D":
            return
        i = -1
        for region in context.area.regions:
            if region.type == "WINDOW":
                i += 1
                if context.region.id == region.id:
                    break
        else:
            return
        rv3d = context.space_data.region_quadviews[i]

    for annotation in _annotations:
        if annotation["type"] == "atom_info":
            _draw_atom_info(region, rv3d, annotation)
        elif annotation["type"] == "bond_angle":
            _draw_bond_angle(region, rv3d, annotation)
        elif annotation["type"] == "com":
            _draw_com(region, rv3d, annotation)
        elif annotation["type"] == "com_distance":
            _draw_com_distance(region, rv3d, annotation)
        elif annotation["type"] == "dihedrals":
            _draw_dihedrals(region, rv3d, annotation)


def _get_universe(obj_name):
    uuid = bpy.context.scene.objects[obj_name].uuid
    trajectory = bpy.context.scene.MNSession.trajectories[uuid]
    return trajectory.universe


def _draw_atom_info(region, rv3d, annotation):
    ag = _get_selection(annotation["obj_name"], annotation["selection"])
    options = annotation["options"]
    for a in ag:
        text = a.name
        if options is not None:
            if options.get("showRes"):
                text += "|res " + str(a.resid) + ", " + a.resname
            if options.get("showSeg"):
                text += "|seg " + a.segid
        _draw_text_3d(region, rv3d, a.position * _bsf, text, options=options)


def _draw_bond_angle(region, rv3d, annotation):
    ag = _get_selection(annotation["obj_name"], annotation["selection"])
    options = annotation["options"]
    for i in range(3):
        a = ag.atoms[i]
        text = a.name
        _draw_text_3d(region, rv3d, a.position * _bsf, text, options=options)
    v1 = ag.atoms[0].position * _bsf
    v2 = ag.atoms[1].position * _bsf
    v3 = ag.atoms[2].position * _bsf
    dist = _distance(v2, v1)
    fdist = (_fmt % (dist / _bsf)) + " Å"
    v21o = _interpolate3d(v2, v1, dist * 0.2)
    _draw_arrow_line_3d(
        region,
        rv3d,
        v2,
        v1,
        text=fdist,
        leftArrow=False,
        rightArrow=True,
        options=options,
    )
    dist = _distance(v2, v3)
    fdist = (_fmt % (dist / _bsf)) + " Å"
    _draw_arrow_line_3d(
        region,
        rv3d,
        v2,
        v3,
        text=fdist,
        leftArrow=False,
        rightArrow=True,
        options=options,
    )
    v23o = _interpolate3d(v2, v3, dist * 0.2)
    angle = _fmt % ag.angle.value() + " °"
    _draw_arrow_line_3d(region, rv3d, v21o, v23o, text=angle, options=options)


def _draw_com(region, rv3d, annotation):
    ag = _get_selection(annotation["obj_name"], annotation["selection"])
    com = ag.center_of_mass() * _bsf
    text = ",".join(annotation["selection"]) + "|" + "COM"
    options = annotation["options"]
    if options is not None:
        text = options.get("text") or text
    _draw_text_3d(region, rv3d, com, text, options=options)


def _draw_com_distance(region, rv3d, annotation):
    ag1 = _get_selection(annotation["obj_name"], annotation["selection1"])
    com1 = ag1.center_of_mass() * _bsf
    text = ",".join(annotation["selection1"]) + "|" + "COM"
    options = annotation["options"]
    if options is not None:
        text = options.get("text1") or text
    _draw_text_3d(region, rv3d, com1, text, options=options)
    ag2 = _get_selection(annotation["obj_name"], annotation["selection2"])
    com2 = ag2.center_of_mass() * _bsf
    text = ",".join(annotation["selection2"]) + "|" + "COM"
    if options is not None:
        text = options.get("text2") or text
    _draw_text_3d(region, rv3d, com2, text, options=options)
    dist = _distance(com1, com2) / _bsf
    fdist = ("%1.2f" % dist) + " Å"
    _draw_arrow_line_3d(
        region,
        rv3d,
        com1,
        com2,
        text=fdist,
        leftArrow=True,
        rightArrow=True,
        options=options,
    )


def _draw_dihedrals(region, rv3d, annotation):
    u = _get_universe(annotation["obj_name"])
    r = u.residues[annotation["resid"] - 1]
    selections = [
        r.phi_selection(),
        r.psi_selection(),
        r.omega_selection(),
        r.chi1_selection(),
    ]
    symbols = ["ϕ", "ψ", "ω", "χ1"]
    options = annotation["options"]
    for i in range(4):
        for j in range(4):
            a = selections[i].atoms[j]
            _draw_text_3d(region, rv3d, a.position * _bsf, a.name, options=options)
            if j != 3:
                b = selections[i].atoms[j + 1]
                text = None
                if j == 1:
                    text = (
                        symbols[i]
                        + " = "
                        + (_fmt % selections[i].dihedral.value())
                        + " °"
                    )
                _draw_arrow_line_3d(
                    region,
                    rv3d,
                    a.position * _bsf,
                    b.position * _bsf,
                    text=text,
                    rightArrow=(j == 1),
                    options=options,
                )


def _in_viewport(pos2d):
    x, y = pos2d
    vw, vh = _viewport
    if x < 0 or x > vw or y < 0 or y > vh:
        return False
    return True


def _get_2d_point(region, rv3d, point3d):
    return view3d_utils.location_3d_to_region_2d(region, rv3d, point3d)


def _draw_text_2d(pos2d, display_text, options=None):
    if pos2d is None or not _in_viewport(pos2d):
        return
    rgba = _default_color
    fsize = _default_font_size
    align = _default_text_align
    text_rot = _default_text_rotation
    if options is not None:
        rgba = options.get("fontColor") or _default_color
        fsize = options.get("fontSize") or _default_font_size
        align = options.get("align") or _default_text_align
        text_rot = options.get("textRotation") or _default_text_rotation
    gap = 12
    x_pos, y_pos = pos2d
    font_id = 0
    ui_scale = bpy.context.preferences.system.ui_scale
    blf.size(font_id, round(fsize * ui_scale))
    # height of one line
    mwidth, mheight = blf.dimensions(font_id, "Tp")  # uses high/low letters
    # split lines
    mylines = display_text.split("|")
    idx = len(mylines) - 1
    # -------------------
    # Draw all text lines
    # -------------------
    for line in mylines:
        text_width, text_height = blf.dimensions(font_id, line)
        if align == "C":
            new_x = x_pos - text_width / 2
        elif align == "R":
            new_x = x_pos - text_width - gap
        else:
            new_x = x_pos
            blf.enable(font_id, blf.ROTATION)
            blf.rotation(font_id, text_rot)
        # calculate new Y position
        new_y = y_pos + (mheight * idx)
        # Draw
        blf.position(font_id, new_x, new_y, 0)
        blf.color(font_id, rgba[0], rgba[1], rgba[2], rgba[3])
        blf.draw(font_id, " " + line)
        # sub line
        idx -= 1
    if align == "L":
        blf.disable(font_id, blf.ROTATION)


def _draw_text_3d(region, rv3d, v, display_text, options=None):
    pos2dv = _get_2d_point(region, rv3d, v)
    _draw_text_2d(pos2dv, display_text, options)


def _distance(v1, v2):
    d = sqrt((v2[0] - v1[0]) ** 2 + (v2[1] - v1[1]) ** 2 + (v2[2] - v1[2]) ** 2)
    return d


def _interpolate3d(v1, v2, d1):
    # calculate vector
    v = (v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2])
    # calculate distance between points
    d0 = _distance(v1, v2)
    # calculate interpolate factor (distance from origin / distance total)
    # if d1 > d0, the point is projected in 3D space
    if d0 > 0:
        x = d1 / d0
    else:
        x = d1
    final = (v1[0] + (v[0] * x), v1[1] + (v[1] * x), v1[2] + (v[2] * x))
    return final


def _draw_line_2d(v1, v2, options=None):
    if v1 is None or v2 is None:
        return
    rgba = _default_color
    lineWidth = _default_line_width
    if options is not None:
        rgba = options.get("lineColor") or _default_color
        lineWidth = options.get("lineWidth") or _default_line_width
    coords = [(v1[0], v1[1], 0), (v2[0], v2[1], 0)]
    gpu.state.blend_set("ALPHA")
    batch = batch_for_shader(_shader_line, "LINES", {"pos": coords})
    try:
        _shader_line.bind()
        _shader_line.uniform_float("color", rgba)
        _shader_line.uniform_float("lineWidth", lineWidth)
        _shader_line.uniform_float("viewportSize", _viewport)
        batch.draw(_shader_line)
    except Exception as e:
        print(e)


def _draw_arrow_line_2d(v1, v2, leftArrow=False, rightArrow=False, options=None):
    if v1 is None or v2 is None:
        return
    arrowSize = _default_arrow_size
    if options is not None:
        arrowSize = options.get("arrowSize") or _default_arrow_size
    _draw_line_2d(v1, v2, options=options)
    if leftArrow:
        v = _interpolate3d((v1[0], v1[1], 0.0), (v2[0], v2[1], 0.0), arrowSize)
        v1i = (v[0] - v1[0], v[1] - v1[1])
        v1a = (
            int(v1i[0] * cos(_rad45) - v1i[1] * sin(_rad45) + v1[0]),
            int(v1i[1] * cos(_rad45) + v1i[0] * sin(_rad45)) + v1[1],
        )
        v1b = (
            int(v1i[0] * cos(_rad315) - v1i[1] * sin(_rad315) + v1[0]),
            int(v1i[1] * cos(_rad315) + v1i[0] * sin(_rad315) + v1[1]),
        )
        _draw_line_2d(v1, v1a, options=options)
        _draw_line_2d(v1, v1b, options=options)
    if rightArrow:
        v = _interpolate3d((v2[0], v2[1], 0.0), (v1[0], v1[1], 0.0), arrowSize)
        v2i = (v[0] - v2[0], v[1] - v2[1])
        v2a = (
            int(v2i[0] * cos(_rad45) - v2i[1] * sin(_rad45) + v2[0]),
            int(v2i[1] * cos(_rad45) + v2i[0] * sin(_rad45)) + v2[1],
        )
        v2b = (
            int(v2i[0] * cos(_rad315) - v2i[1] * sin(_rad315) + v2[0]),
            int(v2i[1] * cos(_rad315) + v2i[0] * sin(_rad315) + v2[1]),
        )
        _draw_line_2d(v2, v2a, options=options)
        _draw_line_2d(v2, v2b, options=options)


def _draw_arrow_line_3d(
    region, rv3d, v1, v2, text=None, leftArrow=False, rightArrow=False, options=None
):
    if v1 is None or v2 is None:
        return
    pos2dv1 = _get_2d_point(region, rv3d, v1)
    pos2dv2 = _get_2d_point(region, rv3d, v2)
    _draw_arrow_line_2d(pos2dv1, pos2dv2, leftArrow, rightArrow, options=options)
    if text is not None:
        dist = _distance(v1, v2)
        midpoint3d = _interpolate3d(v1, v2, fabs(dist / 2))
        _draw_text_3d(region, rv3d, midpoint3d, text, options=options)
