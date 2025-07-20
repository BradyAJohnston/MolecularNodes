import inspect
import types
import typing


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
    supported_types = (str, bool, int, float)
    if atype in supported_types:
        return atype
    if typing.get_origin(atype) is typing.Union or type(atype) is types.UnionType:
        for utype in atype.__args__:
            if utype in supported_types:
                return utype
    return None
