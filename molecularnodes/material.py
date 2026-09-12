from types import ModuleType
from typing import Generic, TypeVar, overload
import bpy
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy import shader as sh
from .nodes.materials import (
    ambient_occlusion,
    default,
    flat,
    flat_outline,
    squishy,
    transparent_outline,
)
from .nodes.shader import ColorAO
from .nodes.shader import FlatInternal as FlatShader
from .nodes.shader import TransparentOutlineInternal as TransparentOutlineShader

# The dumped recipe modules under nodes/materials/ are the source of truth for
# the pre-built materials: the asset library's material datablocks are built
# from the same recipes, so a material created here matches what a user drags
# in from the asset browser.
RECIPES: dict[str, ModuleType] = {
    module.MATERIAL_NAME: module
    for module in (
        default,
        flat,
        flat_outline,
        squishy,
        transparent_outline,
        ambient_occlusion,
    )
}
MATERIAL_NAMES = list(RECIPES)


def _build_recipe(module: ModuleType, name: str | None = None) -> bpy.types.Material:
    "Create a material datablock by running a dumped recipe module's build."
    material = bpy.data.materials.new(
        name if name is not None else module.MATERIAL_NAME
    )
    if material.node_tree is None:  # pragma: no cover - pre-5.x Blender
        material.use_nodes = True
    tree = material.node_tree
    assert tree is not None
    tree.nodes.clear()
    builder = module.MATERIAL.__new__(module.MATERIAL)
    with TreeBuilder(tree) as wrapped:
        builder._build_group(wrapped)
    for key, value in getattr(module, "MATERIAL_PROPERTIES", {}).items():
        setattr(material, key, value)
    return material


def append_material(name: str) -> bpy.types.Material:
    "Get the named pre-built material, creating it from its recipe if needed."
    existing = bpy.data.materials.get(name)
    if existing is not None:
        return existing
    try:
        module = RECIPES[name]
    except KeyError:
        raise KeyError(
            f"No pre-built material {name!r}; available: {MATERIAL_NAMES}"
        ) from None
    return _build_recipe(module)


def add_all_materials() -> dict[str, bpy.types.Material]:
    "Ensure all pre-built materials exist in the file."
    materials = {name: append_material(name) for name in MATERIAL_NAMES}
    # a preset that no style uses yet has zero users, so without a fake user it
    # would be dropped on save/reload or swept up by an orphan purge before the
    # user has had a chance to select it
    for mat in materials.values():
        mat.use_fake_user = True
    return materials


T = TypeVar("T")


class SocketValue(Generic[T]):
    """
    Expose a node input's ``default_value`` as a property on a `PresetMaterial`.

    The descriptor reads and writes the input socket of a node handle stored on
    the preset instance, so parameters remain tweakable after the material has
    been created.
    """

    def __init__(self, node_attr: str, input_name: str, doc: str = ""):
        self._node_attr = node_attr
        self._input_name = input_name
        self.__doc__ = doc or None

    def _socket(self, obj: "PresetMaterial"):
        return getattr(getattr(obj, self._node_attr).i, self._input_name)

    @overload
    def __get__(self, obj: None, objtype: type) -> "SocketValue[T]": ...
    @overload
    def __get__(self, obj: "PresetMaterial", objtype: type | None = ...) -> T: ...

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self
        return self._socket(obj).default_value

    def __set__(self, obj: "PresetMaterial", value: T) -> None:
        self._socket(obj).default_value = value


class PresetMaterial:
    """
    Base class for the pre-built MolecularNodes materials.

    Instantiating a preset creates a new, independent material datablock built
    from the preset's dumped recipe module (``nodes/materials/``) — the same
    recipe that builds the asset library — so tweaking one instance never
    affects the materials of other styles. The key parameter nodes are kept as
    attributes on the instance and exposed as properties that read and write
    the underlying node inputs directly; parameters left as ``None`` keep the
    recipe's values.
    """

    name = "Material"
    recipe: ModuleType

    def _build(self, name: str | None = None) -> None:
        "Create the material datablock and run the preset's recipe into it."
        self._material = _build_recipe(
            type(self).recipe, name if name is not None else type(self).name
        )

    def _handle(self, node_cls):
        """A typed nodebpy handle for the recipe-built node of ``node_cls``,
        so SocketValue properties can address its inputs by name."""
        for node in self._material.node_tree.nodes:
            if node.bl_idname != node_cls._bl_idname:
                continue
            if node.bl_idname == "ShaderNodeGroup":
                group = node.node_tree
                wanted = node_cls._name
                if group is None or not (
                    group.name == wanted or group.name.startswith(wanted + ".")
                ):
                    continue
            return node_cls._from_node(node)
        raise RuntimeError(
            f"{type(self).name}: recipe built no {node_cls.__name__} node"
        )

    @property
    def material(self) -> bpy.types.Material:
        "The underlying Blender material datablock."
        return self._material

    @property
    def tree(self) -> TreeBuilder:
        "The nodebpy tree builder for the material's node tree."
        return TreeBuilder(self._material.node_tree)

    def node(self) -> g.Material:
        "Add a `Material` node to the active GeometryNodeTree and set it to this material."
        return g.Material(material=self.material)

    def __repr__(self) -> str:
        return f"{type(self).__name__}(material={self.material.name!r})"


class Default(PresetMaterial):
    """
    The default MolecularNodes material.

    A principled BSDF with the atom colors darkened by a small amount of
    ambient occlusion, so crevices between atoms read clearly.

    Parameters
    ----------
    roughness : float, optional
        Roughness of the surface. Defaults to the recipe's value.
    ao_distance : float, optional
        Distance (in world units) that the ambient occlusion samples.
        Defaults to the recipe's value.
    ao_exponent : float, optional
        Exponent applied to the occlusion factor; higher values darken
        crevices more aggressively. Defaults to the recipe's value.
    name : str, optional
        Name for the created material datablock. Defaults to ``"Default"``.
    """

    name = "Default"
    recipe = default

    def __init__(
        self,
        roughness: float | None = None,
        ao_distance: float | None = None,
        ao_exponent: float | None = None,
        *,
        name: str | None = None,
    ):
        self._build(name)
        self.ao = self._handle(ColorAO)
        self.bsdf = self._handle(sh.PrincipledBSDF)
        if roughness is not None:
            self.roughness = roughness
        if ao_distance is not None:
            self.ao_distance = ao_distance
        if ao_exponent is not None:
            self.ao_exponent = ao_exponent

    roughness = SocketValue[float]("bsdf", "roughness", "Roughness of the surface.")
    ao_distance = SocketValue[float](
        "ao", "distance", "Distance that the ambient occlusion samples."
    )
    ao_exponent = SocketValue[float](
        "ao", "exponent", "Exponent applied to the occlusion factor."
    )


class AmbientOcclusion(PresetMaterial):
    """
    Colors shaded only by ambient occlusion.

    An emission shader that ignores scene lighting entirely, so it is cheap to
    render and looks similar in Cycles and EEVEE.

    Parameters
    ----------
    distance : float, optional
        Distance in metres (nm) that the ambient occlusion samples.
        Defaults to the recipe's value.
    exponent : float, optional
        Exponent applied to the occlusion factor; higher values darken
        crevices more aggressively. Defaults to the recipe's value.
    name : str, optional
        Name for the created material datablock. Defaults to
        ``"Ambient Occlusion"``.
    """

    name = "Ambient Occlusion"
    recipe = ambient_occlusion

    def __init__(
        self,
        distance: float | None = None,
        exponent: float | None = None,
        *,
        name: str | None = None,
    ):
        self._build(name)
        self.ao = self._handle(ColorAO)
        if distance is not None:
            self.distance = distance
        if exponent is not None:
            self.exponent = exponent

    distance = SocketValue[float](
        "ao", "distance", "Distance that the ambient occlusion samples."
    )
    exponent = SocketValue[float](
        "ao", "exponent", "Exponent applied to the occlusion factor."
    )


class Flat(PresetMaterial):
    """
    Flat cartoon-like shading with an optional dark outline.

    Parameters
    ----------
    outline : bool, optional
        Whether to render the dark outline around the object. Defaults to the
        recipe's value.
    threshold : float, optional
        Threshold for the edge detection that forms the outline. Defaults to
        the recipe's value.
    thickness : float, optional
        Thickness of the outline. Defaults to the recipe's value.
    name : str, optional
        Name for the created material datablock. Defaults to ``"Flat"``.
    """

    name = "Flat"
    recipe = flat

    def __init__(
        self,
        outline: bool | None = None,
        threshold: float | None = None,
        thickness: float | None = None,
        *,
        name: str | None = None,
    ):
        self._build(name)
        self.node = self._handle(FlatShader)
        if outline is not None:
            self.outline = outline
        if threshold is not None:
            self.threshold = threshold
        if thickness is not None:
            self.thickness = thickness

    threshold = SocketValue[float](
        "node", "threshold", "Threshold for the outline edge detection."
    )
    thickness = SocketValue[float]("node", "thickness", "Thickness of the outline.")

    @property
    def outline(self) -> bool:
        "Whether the dark outline is rendered."
        return self.node.i.outline.default_value == "Outline"

    @outline.setter
    def outline(self, value: bool) -> None:
        self.node.i.outline.default_value = "Outline" if value else "None"


class Squishy(PresetMaterial):
    """
    A soft, subsurface-scattering material that makes molecules look jelly-like.

    Parameters
    ----------
    subsurface_scale : float, optional
        Scale of the subsurface scattering radius; larger values look softer.
        Defaults to the recipe's value.
    roughness : float, optional
        Roughness of the surface. Defaults to the recipe's value.
    name : str, optional
        Name for the created material datablock. Defaults to ``"Squishy"``.
    """

    name = "Squishy"
    recipe = squishy

    def __init__(
        self,
        subsurface_scale: float | None = None,
        roughness: float | None = None,
        *,
        name: str | None = None,
    ):
        self._build(name)
        self.bsdf = self._handle(sh.PrincipledBSDF)
        if subsurface_scale is not None:
            self.subsurface_scale = subsurface_scale
        if roughness is not None:
            self.roughness = roughness

    subsurface_scale = SocketValue[float](
        "bsdf", "subsurface_scale", "Scale of the subsurface scattering radius."
    )
    roughness = SocketValue[float]("bsdf", "roughness", "Roughness of the surface.")


class TransparentOutline(PresetMaterial):
    """
    A partially transparent material with an optional solid outline.

    Parameters
    ----------
    alpha : float, optional
        Opacity of the surface; ``0`` is fully transparent, ``1`` fully
        opaque. Defaults to the recipe's value.
    outline : bool, optional
        Whether to render the solid outline around the object. Defaults to
        the recipe's value.
    outline_color : tuple[float, float, float, float], optional
        Color of the outline. Defaults to the recipe's value.
    threshold : float, optional
        Threshold for the edge detection that forms the outline. Defaults to
        the recipe's value.
    thickness : float, optional
        Thickness of the outline. Defaults to the recipe's value.
    name : str, optional
        Name for the created material datablock. Defaults to
        ``"Transparent Outline"``.
    """

    name = "Transparent Outline"
    recipe = transparent_outline

    def __init__(
        self,
        alpha: float | None = None,
        outline: bool | None = None,
        outline_color: tuple[float, float, float, float] | None = None,
        threshold: float | None = None,
        thickness: float | None = None,
        *,
        name: str | None = None,
    ):
        self._build(name)
        self.node = self._handle(TransparentOutlineShader)
        if alpha is not None:
            self.alpha = alpha
        if outline is not None:
            self.outline = outline
        if outline_color is not None:
            self.outline_color = outline_color
        if threshold is not None:
            self.threshold = threshold
        if thickness is not None:
            self.thickness = thickness

    alpha = SocketValue[float]("node", "alpha", "Opacity of the surface.")
    outline_color = SocketValue[tuple[float, float, float, float]](
        "node", "outline_color", "Color of the outline."
    )
    threshold = SocketValue[float](
        "node", "threshold", "Threshold for the outline edge detection."
    )
    thickness = SocketValue[float]("node", "thickness", "Thickness of the outline.")

    @property
    def outline(self) -> bool:
        "Whether the solid outline is rendered."
        return self.node.i.menu.default_value == "Outline"

    @outline.setter
    def outline(self, value: bool) -> None:
        self.node.i.menu.default_value = "Outline" if value else "Transparent"
