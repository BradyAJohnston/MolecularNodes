from contextlib import contextmanager
from pathlib import Path
import bmesh
import bpy
import databpy as db
import numpy as np


def path_resolve(path: str | Path) -> Path:
    if isinstance(path, str):
        return Path(bpy.path.abspath(path))
    elif isinstance(path, Path):
        return Path(bpy.path.abspath(str(path)))
    else:
        raise ValueError(f"Unable to resolve path: {path}")


def set_obj_active(
    obj: bpy.types.Object, context: bpy.types.Context | None = None
) -> None:
    if not context:
        context = bpy.context

    context.view_layer.objects.active = obj  # type: ignore


def set_object_visibility(object: bpy.types.Object, visible: bool) -> None:
    """Set visibility of Blender object"""
    if object.name not in bpy.context.view_layer.objects:
        return
    try:
        object.hide_set(not visible)
        # obj.hide_viewport = hide
        # icon incorrect, but causes blender bug for volumes with above
        object.hide_render = not visible
    except RuntimeError:
        # Keyframing object visibility can at times lead to:
        # RuntimeError: Object can't be hidden because it is not in View Layer 'ViewLayer'!
        pass


def viewport_tag_redraw() -> None:
    """Tag all the viewport areas for a redraw"""
    for window in bpy.context.window_manager.windows:
        for area in window.screen.areas:
            if area.type == "VIEW_3D":
                area.tag_redraw()


def get_viewport_region_from_context(context) -> tuple:
    """Get the 3D viewport region and region data from context"""
    region = context.region
    rv3d = None
    if context.space_data is None:
        return (region, rv3d)
    if not context.space_data.region_quadviews:
        rv3d = context.space_data.region_3d
    else:
        # handle quadview case
        if context.area.type != "VIEW_3D" or context.space_data.type != "VIEW_3D":
            return (region, rv3d)
        i = -1
        for region in context.area.regions:
            if region.type == "WINDOW":
                i += 1
                if context.region == region:
                    break
        else:
            return (region, rv3d)
        rv3d = context.space_data.region_quadviews[i]
    return (region, rv3d)


# the corners of a unit cube centred on the origin, used to grow a point out
# into the geometry drawn around it
_CUBE_CORNERS = np.array(
    [[x, y, z] for x in (-1.0, 1.0) for y in (-1.0, 1.0) for z in (-1.0, 1.0)]
)


def _component_positions(data) -> np.ndarray | None:
    """
    The positions of a geometry component, if it has any.

    A point cloud is drawn as a sphere at every point, so its points are grown
    out by their radius - framing the centres alone cuts the outermost spheres
    in half at the edge of the frame.
    """
    attributes = getattr(data, "attributes", None)
    if attributes is None or "position" not in attributes:
        return None
    positions = db.Attribute(attributes["position"]).as_array()
    if len(positions) == 0:
        return None
    positions = np.asarray(positions, dtype=np.float64).reshape(-1, 3)

    if isinstance(data, bpy.types.PointCloud) and "radius" in attributes:
        radii = np.asarray(
            db.Attribute(attributes["radius"]).as_array(), dtype=np.float64
        ).reshape(-1, 1)
        if radii.any():
            offsets = _CUBE_CORNERS[None, :, :] * radii[:, :, None]
            positions = (positions[:, None, :] + offsets).reshape(-1, 3)

    return positions


def _reference_corners(reference) -> np.ndarray:
    """
    The 8 corners of a box around a geometry that is being instanced.

    Framing on the corners rather than every vertex of every copy keeps the
    point count to 8 per instance while still covering the whole of each one -
    a sphere instanced onto every atom reaches a radius past its centre, and
    framing the centres alone would clip the outermost spheres.
    """
    chunks = []
    for data in (
        getattr(reference, "mesh", None),
        getattr(reference, "pointcloud", None),
        getattr(reference, "curves", None),
    ):
        if data is None:
            continue
        positions = _component_positions(data)
        if positions is not None:
            chunks.append(positions)

    if chunks:
        points = np.concatenate(chunks)
        low, high = points.min(axis=0), points.max(axis=0)
    elif isinstance(reference, bpy.types.Object):
        corners = np.array([co[:] for co in reference.bound_box[:]], dtype=np.float64)
        low, high = corners.min(axis=0), corners.max(axis=0)
    else:
        # a reference we can't read, so treat it as a point at its own origin
        return np.zeros((1, 3))

    return np.array(
        [
            [x, y, z]
            for x in (low[0], high[0])
            for y in (low[1], high[1])
            for z in (low[2], high[2])
        ]
    )


def _instance_points(instances, references) -> np.ndarray:
    """Place a box around each instanced geometry at each of its instances."""
    transforms = np.asarray(
        db.Attribute(instances.attributes["instance_transform"]).as_array(),
        dtype=np.float64,
    )
    indices = np.asarray(
        db.Attribute(instances.attributes[".reference_index"]).as_array(), dtype=int
    )

    chunks = []
    for index, reference in enumerate(references):
        placements = transforms[indices == index]
        if len(placements) == 0:
            continue
        corners = _reference_corners(reference)
        corners = np.hstack([corners, np.ones((len(corners), 1))])
        # the transforms are stored row-vector style, with the translation in
        # the last row rather than the last column
        placed = np.einsum("kj,nji->nki", corners, placements)[..., :3]
        chunks.append(placed.reshape(-1, 3))

    if not chunks:
        return np.zeros((0, 3))
    return np.concatenate(chunks)


def evaluated_points(obj: bpy.types.Object) -> np.ndarray:
    """
    The world-space positions of the geometry an object actually renders.

    Reads back the evaluated geometry, so what comes out is what a style built -
    the vertices of a cartoon ribbon, say - rather than the atoms it was built
    from. This is what makes framing match what ends up on screen.

    Parameters
    ----------
    obj : bpy.types.Object
        The object to read.

    Returns
    -------
    ndarray
        An ``(N, 3)`` array of positions. Falls back to the corners of the
        object's bounds if it renders nothing readable.

    Notes
    -----
    Geometry nodes can hand back several components at once, and an object only
    ever exposes one of them through its ``data`` - a molecule styled with
    spheres evaluates to an empty mesh with the actual point cloud alongside it,
    invisible from ``obj.data``. Every component is collected here, instanced
    geometry included.
    """
    # evaluated geometry is stale until the depsgraph catches up
    bpy.context.view_layer.update()
    geometry = db.GeometrySet(obj)

    chunks = []
    for name, data in geometry.components().items():
        if name == "INSTANCES":
            chunks.append(_instance_points(data, geometry.instance_references))
        else:
            positions = _component_positions(data)
            if positions is not None:
                chunks.append(positions)

    chunks = [chunk for chunk in chunks if len(chunk) > 0]
    if chunks:
        points = np.concatenate(chunks)
    else:
        # nothing readable, so fall back to the bounds Blender reports, which
        # do account for every component
        points = np.array([co[:] for co in obj.bound_box[:]], dtype=np.float64)

    matrix = np.array(geometry.evaluated_object.matrix_world)
    return points @ matrix[:3, :3].T + matrix[:3, 3]


@contextmanager
def new_bmesh() -> bmesh.types.BMesh:  # type: ignore
    # create a new bmesh
    bm = bmesh.new()
    try:
        # return the new bmesh
        yield bm
    finally:
        # free bmesh
        bm.free()
