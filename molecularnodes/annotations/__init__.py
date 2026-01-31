import importlib
import sys
from . import base, interface, manager
from .utils import get_all_subclasses


def _reregister_annotations():
    subclasses = get_all_subclasses(manager.BaseAnnotationManager)
    # reload required annotation modules
    importlib.reload(sys.modules[manager.__name__])
    for cls in subclasses:
        importlib.reload(sys.modules[cls.__module__])


__all__ = [
    "base",
    "interface",
    "manager",
]
