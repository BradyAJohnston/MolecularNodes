from typing import Generic, TypeVar, overload
import bpy
from databpy.material import append_from_blend
from nodebpy import geometry as g
from nodebpy import shader as sh
from nodebpy.builder import MaterialBuilder
from ..assets import MN_DATA_FILE
from .interface import remove_linked
from .shader import ColorAO, MNColor
from .shader import Flat as FlatShader
from .shader import TransparentOutline as TransparentOutlineShader

MATERIAL_NAMES = [
    "MN Default",
    "MN Flat Outline",
    "MN Squishy",
    "MN Transparent Outline",
    "MN Ambient Occlusion",
]


def append_material(name: str) -> bpy.types.Material:
    "Append a material from the MN_DATA_FILE."
    return append_from_blend(name, str(MN_DATA_FILE))


def add_all_materials() -> dict[str, bpy.types.Material]:
    "Append all pre-defined materials from the MN_DATA_FILE."
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
    with `nodebpy`, so tweaking one instance never affects the materials of
    other styles. The nodes used to build the material are kept as attributes
    on the instance and the key parameters are exposed as properties that read
    and write the underlying node inputs directly.
    """

    name = "Material"

    def _create(self, name: str | None = None) -> MaterialBuilder:
        "Create the material datablock and return its tree builder for a `with` block."
        self._tree = sh.material(name if name is not None else type(self).name)
        self._tree.nodes.clear()
        return self._tree

    @property
    def material(self) -> bpy.types.Material:
        "The underlying Blender material datablock."
        return self._tree.material

    @property
    def tree(self) -> MaterialBuilder:
        "The nodebpy tree builder for the material's node tree."
        return self._tree

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
    roughness : float
        Roughness of the surface.
    ao_distance : float
        Distance (in world units) that the ambient occlusion samples.
    ao_exponent : float
        Exponent applied to the occlusion factor; higher values darken
        crevices more aggressively.
    name : str, optional
        Name for the created material datablock. Defaults to ``"Default"``.
    """

    name = "Default"

    def __init__(
        self,
        roughness: float = 0.264,
        ao_distance: float = 0.5,
        ao_exponent: float = 1.0,
        *,
        name: str | None = None,
    ):
        with self._create(name):
            color = MNColor()
            self.ao = ColorAO(color.o.color, distance=ao_distance, exponent=ao_exponent)
            self.bsdf = sh.PrincipledBSDF(
                base_color=self.ao.o.result,
                roughness=roughness,
                ior=1.45,
                alpha=color.o.alpha,
            )
            self.bsdf.o.bsdf >> sh.MaterialOutput().i.surface

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
    distance : float
        Distance in metres (nm) that the ambient occlusion samples.
    exponent : float
        Exponent applied to the occlusion factor; higher values darken
        crevices more aggressively.
    name : str, optional
        Name for the created material datablock. Defaults to
        ``"Ambient Occlusion"``.
    """

    name = "Ambient Occlusion"

    def __init__(
        self,
        distance: float = 1.0,
        exponent: float = 2.0,
        *,
        name: str | None = None,
    ):
        with self._create(name):
            color = MNColor()
            self.ao = ColorAO(color.o.color, distance=distance, exponent=exponent)
            emission = sh.Emission(color=self.ao.o.result)
            emission.o.emission >> sh.MaterialOutput().i.surface

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
    outline : bool
        Whether to render the dark outline around the object.
    threshold : float
        Threshold for the edge detection that forms the outline.
    thickness : float
        Thickness of the outline.
    name : str, optional
        Name for the created material datablock. Defaults to ``"Flat"``.
    """

    name = "Flat"

    def __init__(
        self,
        outline: bool = True,
        threshold: float = 0.8,
        thickness: float = 0.15,
        *,
        name: str | None = None,
    ):
        with self._create(name):
            self.node = FlatShader(
                outline="Outline" if outline else "None",
                threshold=threshold,
                thickness=thickness,
            )
            self.node.o.emission >> sh.MaterialOutput().i.surface

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
    subsurface_scale : float
        Scale of the subsurface scattering radius; larger values look softer.
    roughness : float
        Roughness of the surface.
    name : str, optional
        Name for the created material datablock. Defaults to ``"Squishy"``.
    """

    name = "Squishy"

    def __init__(
        self,
        subsurface_scale: float = 0.2,
        roughness: float = 1.0,
        *,
        name: str | None = None,
    ):
        with self._create(name):
            color = MNColor()
            self.bsdf = sh.PrincipledBSDF(
                base_color=color.o.color,
                alpha=color.o.alpha,
                roughness=roughness,
                ior=1.05,
                diffuse_roughness=1.0,
                subsurface_weight=1.0,
                subsurface_radius=(1.0, 0.2, 0.1),
                subsurface_scale=subsurface_scale,
                coat_weight=1.0,
                coat_roughness=0.246,
            )
            self.bsdf.o.bsdf >> sh.MaterialOutput().i.surface

    subsurface_scale = SocketValue[float](
        "bsdf", "subsurface_scale", "Scale of the subsurface scattering radius."
    )
    roughness = SocketValue[float]("bsdf", "roughness", "Roughness of the surface.")


class TransparentOutline(PresetMaterial):
    """
    A partially transparent material with an optional solid outline.

    Parameters
    ----------
    alpha : float
        Opacity of the surface; ``0`` is fully transparent, ``1`` fully opaque.
    outline : bool
        Whether to render the solid outline around the object.
    outline_color : tuple[float, float, float, float]
        Color of the outline.
    threshold : float
        Threshold for the edge detection that forms the outline.
    thickness : float
        Thickness of the outline.
    name : str, optional
        Name for the created material datablock. Defaults to
        ``"Transparent Outline"``.
    """

    name = "Transparent Outline"

    def __init__(
        self,
        alpha: float = 0.7,
        outline: bool = True,
        outline_color: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 1.0),
        threshold: float = 0.2,
        thickness: float = 0.15,
        *,
        name: str | None = None,
    ):
        with self._create(name):
            self.node = TransparentOutlineShader(
                alpha=alpha,
                menu="Outline" if outline else "Transparent",
                outline_color=outline_color,
                threshold=threshold,
                thickness=thickness,
            )
            self.node.o.shader >> sh.MaterialOutput().i.surface
        self.material.surface_render_method = "BLENDED"

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


def set_socket_material(
    socket: bpy.types.NodeSocketMaterial,
    mat: bpy.types.Material
    | bpy.types.NodeSocketMaterial
    | PresetMaterial
    | MaterialBuilder
    | str
    | None,
) -> None:
    remove_linked(socket)
    if mat is None:
        socket.default_value = None
    elif isinstance(mat, (PresetMaterial, MaterialBuilder)):
        socket.default_value = mat.material
    elif isinstance(mat, bpy.types.Material):
        socket.default_value = mat
    elif isinstance(mat, bpy.types.NodeSocketMaterial):
        socket.node.id_data.links.new(mat, socket)  # type: ignore
    elif isinstance(mat, str):
        mat = bpy.data.materials[mat]
        socket.default_value = mat
    else:
        raise TypeError("Invalid type for setting of a material: " + str(type(mat)))


def assign_material(
    node: bpy.types.GeometryNodeGroup,
    new_material: bpy.types.Material
    | bpy.types.NodeSocketMaterial
    | PresetMaterial
    | MaterialBuilder
    | str
    | None = "default",
) -> None:
    add_all_materials()

    if isinstance(new_material, str):
        if new_material not in bpy.data.materials:
            try:
                append_material(new_material)
            except Exception:
                try:
                    new_material = "MN " + new_material.title().strip()
                    append_material(new_material)
                except Exception:
                    raise ValueError(
                        f"Material {new_material} not found in this file of the included MN preset file."
                    )
    try:
        set_socket_material(
            socket=node.inputs["Material"],
            mat=new_material,
        )
    except KeyError:
        return "Material input not found on node."
