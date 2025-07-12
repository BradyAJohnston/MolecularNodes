import inspect


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
