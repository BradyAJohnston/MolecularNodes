import functools
import inspect
from abc import ABCMeta
from uuid import uuid1
import bpy
from databpy.object import LinkedObjectError
from mathutils import Vector
from ..blender.utils import get_viewport_region_from_context, viewport_tag_redraw
from .base import BaseAnnotation
from .interface import AnnotationInterface
from .props import (
    BaseAnnotationProperties,
    create_annotation_type_inputs,
    create_property_interface,
)
from .utils import (
    get_all_class_annotations,
    get_blender_supported_type,
    get_view_matrix,
    is_perspective_projection,
)


def _validate_annotation_update(self, context, attr):
    """Update callback when annotation inputs change (both API and GUI)"""
    session = context.scene.MNSession
    entity = session.get(self.id_data.uuid)
    interface = entity.annotations._interfaces.get(self.uuid)
    instance = interface._instance
    # delete non blender attribute as blender attribute updated
    nbattr = f"_{attr}"
    if hasattr(instance, nbattr):
        delattr(instance, nbattr)
    try:
        if not instance.validate():
            raise ValueError(f"Invalid input {attr}")
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
    set the class attribute _entity_type to the entity tpye, the class attribute
    '_classes' and instance attribute '_interfaces' to empty dictionaries. Derived
    classes will have to pass the entity instance as part of its '__init__' as well.

    Attributes
    ----------
    visible: bool
        Whether to show or hide all annotations

    """

    _entity_type = None  # Derived classes need to specify
    _classes = {}  # All annotation classes across all entities

    def __init__(self, entity):
        # Access to the entity to which this manager is attached
        self._entity = entity
        self._draw_handler = None
        self._scene = None
        self._render_mode = False
        # viewport params
        self._region = None
        self._rv3d = None
        # render params
        self._render_scale = 1.0
        self._image = None
        self._image_scale = 1
        # add draw handler by default
        self._draw_handler_add()

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
        annotation_type = annotation_class.annotation_type
        if annotation_type in cls._classes:
            raise ValueError(f"Annotation type {annotation_type} already registered")
        cls._classes[annotation_type] = annotation_class
        # Add to all annotation classes
        if cls._entity_type is None:
            raise ValueError("Entity type is needed in derived class")
        entity_annotation_type = f"{cls._entity_type.value}_{annotation_type}"
        if entity_annotation_type in BaseAnnotationManager._classes:
            raise ValueError(
                f"Annotation type {entity_annotation_type} already registered"
            )
        BaseAnnotationManager._classes[entity_annotation_type] = annotation_class
        # Add method to Entity specific manager
        method_name = f"add_{annotation_type}"

        # Add a dynamic wrapper to ensure custom signature is retained
        def dynamic_wrapper(cls, annotation_class, **kwargs):
            return cls._create_annotation_instance(annotation_class, **kwargs)

        method = functools.partialmethod(dynamic_wrapper, annotation_class)
        parameters = [
            inspect.Parameter("self", inspect.Parameter.POSITIONAL_ONLY),
            inspect.Parameter("annotation_class", inspect.Parameter.POSITIONAL_ONLY),
        ]
        # Add the annotation inputs as keyword params
        py_annotations = get_all_class_annotations(annotation_class)
        for attr, atype in py_annotations.items():
            default = inspect.Parameter.empty
            # set the default if the param in annotation class is a class
            # attribute and not just an annotation
            if hasattr(annotation_class, attr):
                default = getattr(annotation_class, attr)
            param = inspect.Parameter(
                attr, inspect.Parameter.KEYWORD_ONLY, annotation=atype, default=default
            )
            parameters.append(param)
        # Create method signature to match annotation class inputs
        method.func.__signature__ = inspect.Signature(
            parameters, return_annotation=AnnotationInterface
        )
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
        for annotation_type, annotation_class in BaseAnnotationManager._classes.items():
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
        annotation_type = annotation_class.annotation_type
        if annotation_type not in cls._classes:
            raise ValueError(f"{annotation_class} is not registered")
        cls._validate_annotation_class(annotation_class)
        # Delete from Entity class specific registry
        del cls._classes[annotation_type]
        # Delete from all annotation classes
        entity_annotation_type = f"{cls._entity_type.value}_{annotation_type}"
        del BaseAnnotationManager._classes[entity_annotation_type]
        # Remove method from Entity specific manager
        method_name = f"add_{annotation_type}"
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
        if value:
            self._draw_handler_add()
        else:
            self._draw_handler_remove()
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
        for attr in py_annotations.keys():
            if not hasattr(annotation_class, attr) and attr not in kwargs:
                raise ValueError(f"{attr} is a required parameter")
        for attr in kwargs.keys():
            if attr not in py_annotations:
                raise ValueError(
                    f"Unknown input {attr}. Valid values are {py_annotations}"
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
            for attr in py_annotations.keys():
                # set the input value based on what is passed
                # if input value is not passed, use the value from the annotation class
                value = kwargs.get(attr, None)
                if value is None:
                    value = getattr(annotation_class, attr, None)
                if value is not None:
                    setattr(interface.__class__, attr, value)
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
        prop_interface = create_property_interface(self._entity, uuid, "label")
        setattr(interface.__class__, "name", prop_interface)
        # iterate though all the annotation inputs, set passed values and
        # create property interfaces
        entity_annotation_type = f"{self._entity_type.value}_{prop.type}"
        inputs = getattr(prop, entity_annotation_type, None)
        if inputs is not None:
            inputs.uuid = uuid  # add annotation uuid for lookup in update callback
            # link to the valid inputs property for use in draw handler
            prop_interface = create_property_interface(
                self._entity,
                uuid,
                "valid_inputs",
                annotation_type=entity_annotation_type,
            )
            setattr(interface.__class__, "_valid_inputs", prop_interface)
            # all annotation inputs as defined in class
            for attr, atype in py_annotations.items():
                # name is a special case that is already added above
                if attr == "name":
                    continue
                stype = get_blender_supported_type(atype)
                prop_interface = create_property_interface(
                    self._entity,
                    uuid,
                    attr,
                    stype,
                    annotation_instance,
                    entity_annotation_type,
                )
                setattr(interface.__class__, attr, prop_interface)
                # set the input value based on what is passed
                # if input value is not passed, use the value from the annotation class
                value = kwargs.get(attr, None)
                if value is None:
                    value = getattr(annotation_class, attr, None)
                if value is not None:
                    setattr(interface, attr, value)
        # create property interfaces for the common annotation params
        for item in prop.bl_rna.properties:
            if not item.is_runtime:
                continue
            prop_path = item.identifier
            if prop_path in ("name", "label", "type"):
                continue
            if item.type == "POINTER":
                continue
            prop_interface = create_property_interface(self._entity, uuid, prop_path)
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

    def _draw_handler_add(self):
        if bpy.app.background:
            return
        if self._draw_handler is not None:
            return
        self._draw_handler = bpy.types.SpaceView3D.draw_handler_add(
            self._draw_annotations_handler,
            (bpy.context,),
            "WINDOW",
            "POST_PIXEL",
        )
        self._redraw()

    def _draw_handler_remove(self):
        if self._draw_handler is not None:
            bpy.types.SpaceView3D.draw_handler_remove(self._draw_handler, "WINDOW")
            self._redraw()
            self._draw_handler = None

    def _redraw(self):
        if self._draw_handler is not None:
            viewport_tag_redraw()

    def _is_valid_entity(self) -> bool:
        try:
            _name = self._entity.name
            return True
        except LinkedObjectError:
            # remove any registered draw handler
            self._draw_handler_remove()
            return False

    def _draw_annotations_handler(self, context):
        if self._draw_handler is None:
            return
        if not self._is_valid_entity():
            return
        region, rv3d = get_viewport_region_from_context(context)
        # for viewport drawing, region and rv3d are required
        if region is None or rv3d is None:
            return
        self._region = region
        self._rv3d = rv3d
        self._scene = context.scene
        self._draw_annotations()

    def _enable_render_mode(self, scene, render_scale, image, image_scale):
        self._scene = scene
        self._render_scale = render_scale
        self._image = image
        self._image_scale = image_scale
        self._render_mode = True

    def _disable_render_mode(self):
        self._render_mode = False

    def _draw_annotations(self):
        # all annotations visibilty
        if not self.visible:
            return
        # object visibility
        object = self._entity.object
        if object.hide_get():
            return
        # check if object is in current scene
        if object.name not in self._scene.objects:
            return
        # calculate the min and max viewport distance of the bounding box verts
        # find bounding box verts
        bb_verts = [co for co in bpy.data.objects[object.name].bound_box]
        # viewport distance of the bounding box verts
        bb_verts_distances = []
        for vert in bb_verts:
            view_matrix = get_view_matrix(self)
            if is_perspective_projection(self):
                # perspective
                dist = (view_matrix @ Vector(vert)).length
            else:
                # orthographic
                dist = -(view_matrix @ Vector(vert)).z
            bb_verts_distances.append(dist)
        # min and max
        min_dist = min(bb_verts_distances)
        max_dist = max(bb_verts_distances)
        # iterate over all annotations
        for interface in self._interfaces.values():
            # annotation specific visibility
            if not interface.visible:
                continue
            # annotations input validity
            if not getattr(interface, "_valid_inputs", True):
                continue
            # set distance params
            interface._instance._set_distance_params(min_dist, max_dist)
            if self._render_mode:
                # set the render params
                interface._instance._set_render_params(
                    self._scene, self._render_scale, self._image, self._image_scale
                )
            else:
                # set the viewport region params
                interface._instance._set_viewport_params(
                    self._scene, self._region, self._rv3d
                )
            # handle exceptions to allow other annotations to be drawn
            try:
                interface._instance.draw()
            except Exception:
                pass
