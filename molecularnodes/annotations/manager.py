import functools
import inspect
from abc import ABCMeta
from uuid import uuid1
import bpy
from ..blender.utils import viewport_tag_redraw
from .base import BaseAnnotation
from .interface import AnnotationInterface
from .props import (
    BaseAnnotationProperties,
    create_annotation_type_inputs,
    create_property_interface,
)
from .utils import get_all_class_annotations


def _validate_annotation_update(self, context, prop_name):
    """Update callback when annotation inputs change (both API and GUI)"""
    session = context.scene.MNSession
    entity = session.get(self.id_data.uuid)
    interface = entity.annotations._interfaces.get(self.uuid)
    try:
        if not interface._instance.validate():
            raise ValueError(f"Invalid input {prop_name}")
    except Exception as exception:
        self.valid_inputs = False
        raise exception
    self.valid_inputs = True


class BaseAnnotationManager(metaclass=ABCMeta):
    """
    Base class for Annotation Manager

    This is the base class for the Annotation Manager, which manages all the
    annotation classes, instances and interfaces

    Entities that need annotation support have to derive from this base class and
    set the class attribute '_classes' and instance attribute '_interfaces' to
    empty dictionaries. Derived classes will have to pass the entity instance as
    part of its '__init__' as well.

    Attributes
    ----------
    visible: bool
        Whether to show or hide all annotations

    """

    def __init__(self, entity):
        # Access to the entity to which this manager is attached
        self._entity = entity

    def __iter__(self) -> iter:
        """To support iteration"""
        return iter(self._interfaces.values())

    def __len__(self) -> int:
        """Number of annotations"""
        return len(self._interfaces)

    def __getitem__(self, name: str) -> AnnotationInterface:
        """Return annotation interface based on name or index"""
        if isinstance(name, int):
            # numeric index based array access
            index = int(name)
            max = len(self._interfaces)
            if 0 <= index < max:
                key = list(self._interfaces.keys())[index]
                return self._interfaces[key]
            else:
                raise ValueError(f"Invalid index {index} of {max} items")
        return self.get(name)

    def _ipython_key_completions_(self) -> list[str]:
        """Return annotation names"""
        return [i.name for i in self._interfaces.values()]

    @classmethod
    def register(cls, annotation_class) -> None:
        """
        Register an annotation class

        This method adds the annotation class to the entity specific class
        registry (_classes) and adds a new method (add_<annotation_type>) to
        the manager with a signature that matches the annotation inputs.

        """
        cls._validate_annotation_class(annotation_class)
        # Add to Entity class specific registry
        cls._classes[annotation_class.annotation_type] = annotation_class
        # Add method to Entity specific manager
        method_name = f"add_{annotation_class.annotation_type}"
        method = functools.partialmethod(
            cls._create_annotation_instance, annotation_class
        )
        parameters = [
            inspect.Parameter("self", inspect.Parameter.POSITIONAL_ONLY),
            inspect.Parameter("annotation_class", inspect.Parameter.POSITIONAL_ONLY),
        ]
        # Add the annotation inputs as keyword params
        py_annotations = get_all_class_annotations(annotation_class)
        for key, value in py_annotations.items():
            default = inspect.Parameter.empty
            # set the default if the param in annotation class is a class
            # attribute and not just an annotation
            if hasattr(annotation_class, key):
                default = getattr(annotation_class, key)
            param = inspect.Parameter(
                key, inspect.Parameter.KEYWORD_ONLY, annotation=value, default=default
            )
            parameters.append(param)
        # Create method signature to match annotation class inputs
        method.func.__signature__ = inspect.Signature(parameters)
        # Update annotation properties
        cls._update_annotation_props(annotation_class)
        setattr(cls, method_name, method)

    @classmethod
    def _update_annotation_props(cls, annotation_class: BaseAnnotation):
        """Update annotation properties attached to Object"""
        AnnotationProperties = type(
            "AnnotationProperties", (BaseAnnotationProperties,), {}
        )
        # Add each annotation type inputs as a pointer to a separate property group
        for annotation_class in cls._classes.values():
            annotation_type = annotation_class.annotation_type
            update_callback = None
            # Add an update callback for annotation input properties to validate
            if hasattr(annotation_class, "validate") and callable(
                getattr(annotation_class, "validate")
            ):
                update_callback = _validate_annotation_update
            AnnotationInputs = create_annotation_type_inputs(
                annotation_class, update_callback=update_callback
            )
            bpy.utils.register_class(AnnotationInputs)
            AnnotationProperties.__annotations__[annotation_type] = (
                bpy.props.PointerProperty(type=AnnotationInputs)
            )
        # Re-register the new AnnotationProperties class
        bpy.utils.register_class(AnnotationProperties)
        # Re-assign the annotation properties to Object - old data is retained
        bpy.types.Object.mn_annotations = bpy.props.CollectionProperty(
            type=AnnotationProperties
        )

    @classmethod
    def unregister(cls, annotation_class) -> None:
        """
        Unregister a registered annotation class

        This method removes the annotation class from the entity speicific
        class registry and removes the 'add_<>' method from the manager

        """
        if annotation_class.annotation_type not in cls._classes:
            raise ValueError(f"{annotation_class} is not registered")
        cls._validate_annotation_class(annotation_class)
        del cls._classes[annotation_class.annotation_type]
        # Remove method from Entity specific manager
        method_name = f"add_{annotation_class.annotation_type}"
        delattr(cls, method_name)

    def get(self, name: str) -> AnnotationInterface:
        """Get an annotation by name"""
        for instance in self._interfaces.values():
            if instance.name == name:
                return instance
        raise ValueError(f"Annotation with name '{name}' not found")

    def remove(self, annotation: str | AnnotationInterface) -> None:
        """
        Remove an annotation by name or instance

        When a name is used, all annotations that match the name will be removed

        """
        instances = []
        if isinstance(annotation, str):
            for instance in self._interfaces.values():
                if instance.name == annotation:
                    instances.append(instance)
        else:
            if not isinstance(annotation, AnnotationInterface):
                raise ValueError("Invalid annotation instance")
            if annotation._uuid not in self._interfaces:
                raise ValueError(f"Instance for {annotation._uuid} not found")
            instances.append(annotation)
        if not instances:
            raise ValueError(f"No annotations found for {annotation}")
        for instance in instances:
            self._remove_annotation_instance(instance)

    def clear(self) -> None:
        """Remove all annotations"""
        instances = []
        for instance in self._interfaces.values():
            instances.append(instance)
        for instance in instances:
            self._remove_annotation_instance(instance)

    @property
    def visible(self) -> bool:
        """Visibility of all annotations - getter"""
        return self._entity.object.mn.annotations_visible

    @visible.setter
    def visible(self, value: bool) -> None:
        """Visibility of all annotations - setter"""
        self._entity.object.mn.annotations_visible = value
        viewport_tag_redraw()

    @classmethod
    def _validate_annotation_class(cls, annotation_class: BaseAnnotation) -> None:
        """Internal method to validate an annotation class"""
        if not issubclass(annotation_class, BaseAnnotation):
            raise ValueError(f"{annotation_class} is not an annotation class")
        if not hasattr(annotation_class, "annotation_type"):
            raise ValueError("annotation_type not specified in class")
        if annotation_class.draw.__qualname__ == "BaseAnnotation.draw":
            raise ValueError("draw method not found in class")

    def _create_annotation_instance(
        self, annotation_class: BaseAnnotation, **kwargs
    ) -> AnnotationInterface:
        """
        Create annotation instance and return a dynamic Annotation Interface

        This is the method that is bound to the 'add_<annotation_type>' methods
        of the manager. This method is responsible for creating the actual
        annotation instance, creating a Geometry Node with inputs that match the
        annotation type and creating an interface that binds the properties to
        the corresponding Annotation Geometry Node inputs.

        """
        # validations for required and invalid inputs
        py_annotations = get_all_class_annotations(annotation_class)
        for key, value in py_annotations.items():
            if not hasattr(annotation_class, key) and key not in kwargs:
                raise ValueError(f"{key} is a required parameter")
        for key, value in kwargs.items():
            if key not in py_annotations:
                raise ValueError(
                    f"Unknown input {key}. Valid values are {py_annotations}"
                )
        # create an annotations instance
        annotation_instance = annotation_class(self._entity)
        # create a new dynamic interface class
        inteface_class_name = f"{annotation_class.__name__}_interface"
        DynamicInterface = type(inteface_class_name, (AnnotationInterface,), {})
        interface = DynamicInterface(annotation_instance)
        # set the interface attribute of the annotation instance so that
        # users can access the inputs and common params using this interface
        setattr(annotation_instance, "interface", interface)
        # call the validate method in the annotation class if specified for
        # any annotation specific custom validation
        if hasattr(annotation_instance, "validate") and callable(
            getattr(annotation_instance, "validate")
        ):
            # set the annotation inputs to interface for validation
            for key, value in py_annotations.items():
                # set the input value based on what is passed
                # if input value is not passed, use the value from the annotation class
                value = kwargs.get(key, None)
                if value is None:
                    value = getattr(annotation_class, key, None)
                if value is not None:
                    setattr(interface.__class__, key, value)
            # validate
            if not annotation_instance.validate():
                raise ValueError("Invalid annotation inputs")
        # only after all validations pass start doing real stuff like creating
        # properties and adding to the interface list
        uuid = str(uuid1())
        # add the new interface to the entity specific interfaces registry
        self._interfaces[uuid] = interface
        # create a new annotation property
        object = self._entity.object
        prop = object.mn_annotations.add()
        # use the instance uuid as name for later lookups using find()
        prop.name = uuid
        prop.type = annotation_class.annotation_type
        # set the annotations active index for GUI
        object.mn.annotations_active_index = len(object.mn_annotations) - 1
        # add the uuid as internal attribute
        # _uuid will be used to lookup this interface in the interfaces registry
        setattr(interface, "_uuid", uuid)
        # set the annotation name and create property interface
        value = kwargs.get("name", None)
        if value is not None:
            setattr(prop, "label", value)
        else:
            if object.mn.annotations_next_index:
                prop.label = f"Annotation.{object.mn.annotations_next_index:03d}"
            else:
                prop.label = "Annotation"
            object.mn.annotations_next_index += 1
        prop_interface = create_property_interface(prop, "label")
        setattr(interface.__class__, "name", prop_interface)
        # iterate though all the annotation inputs, set passed values and
        # create property interfaces
        inputs = getattr(prop, prop.type, None)
        if inputs is not None:
            inputs.uuid = uuid  # add annotation uuid for lookup in update callback
            for key, value in py_annotations.items():
                # name is a special case that is already added above
                if key == "name":
                    continue
                # set the input value based on what is passed
                # if input value is not passed, use the value from the annotation class
                value = kwargs.get(key, None)
                if value is None:
                    value = getattr(annotation_class, key, None)
                if value is not None:
                    setattr(inputs, key, value)
                prop_interface = create_property_interface(inputs, key)
                setattr(interface.__class__, key, prop_interface)
        # create property interfaces for the common annotation params
        for item in prop.bl_rna.properties:
            if not item.is_runtime:
                continue
            prop_path = item.identifier
            if prop_path in ("name", "label", "type"):
                continue
            if item.type == "POINTER":
                continue
            prop_interface = create_property_interface(prop, prop_path)
            setattr(interface.__class__, prop_path, prop_interface)
        # call the defaults method in the annotation class if specified
        if hasattr(annotation_instance, "defaults") and callable(
            getattr(annotation_instance, "defaults")
        ):
            annotation_instance.defaults()
        # tag redraw of viewport to refresh values in GUI
        viewport_tag_redraw()
        # return the newly created dynamic annotation interface
        return interface

    def _remove_annotation_instance(self, instance) -> None:
        """Actual method to remove annotation instance"""
        uuid = instance._uuid
        object = self._entity.object
        index = object.mn_annotations.find(uuid)
        if index != -1:
            object.mn_annotations.remove(index)
        object.mn.annotations_active_index = len(object.mn_annotations) - 1
        # tag redraw of viewport to refresh values in GUI
        viewport_tag_redraw()
        # remove from instances registry
        del self._interfaces[uuid]

    def _remove_annotation_by_uuid(self, uuid: str) -> None:
        """Remove annotation instance by uuid - used by GUI"""
        if uuid not in self._interfaces:
            raise ValueError(f"Instance for {uuid} not found")
        self._remove_annotation_instance(self._interfaces[uuid])
