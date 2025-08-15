import inspect
import types
import typing
import bpy
from PIL import Image


def get_all_class_annotations(cls) -> dict:
    """Get all the annotations from the class including base classes"""
    all_annotations = {}
    for base_cls in cls.mro():
        if hasattr(base_cls, "__annotations__"):
            current_annotations = inspect.get_annotations(base_cls, eval_str=True)
            all_annotations.update(current_annotations)
    if "interface" in all_annotations:
        del all_annotations["interface"]
    return all_annotations


def get_blender_supported_type(atype):
    supported_types = (
        str,
        bool,
        int,
        float,
        tuple[float, float],
        tuple[float, float, float],
        tuple[float, float, float, float],
    )
    if atype in supported_types:
        return atype
    if typing.get_origin(atype) is typing.Union or type(atype) is types.UnionType:
        for utype in atype.__args__:
            if utype in supported_types:
                return utype
    return None


def get_view_matrix(obj):
    if obj._render_mode:
        return obj._scene.camera.matrix_world.inverted()
    else:
        return obj._rv3d.view_matrix


def is_perspective_projection(obj):
    if obj._render_mode:
        return obj._scene.camera.data.type != "ORTHO"
    else:
        return obj._rv3d.is_perspective


def render_annotations(
    scene: bpy.types.Scene, render_scale: float, image: Image, image_scale: float
) -> None:
    """Render annotations of all entities to an image"""
    session = scene.MNSession
    session.prune()  # remove any invalid session entities
    for entity in session.entities.values():
        if not hasattr(entity, "annotations"):
            continue
        manager = entity.annotations
        manager._enable_render_mode(scene, render_scale, image, image_scale)
        manager._draw_annotations()
        manager._disable_render_mode()
