# Node-group asset "Style Spheres" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MaterialSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputMaterial,
    InputMenu,
)
from ._shared.mn_utils_style_spheres_icosphere import MN_utils_style_spheres_icosphere
from ._shared.mn_utils_style_spheres_points import MN_utils_style_spheres_points
from .evaluate_on_atoms import EvaluateOnAtoms


class StyleSpheres(AssetGeometryGroup):
    """
    Style Spheres

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this style to, discarding unselected points
    sphere : InputMenu | Literal["Point", "Instance", "Mesh"]
        Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
    quality : InputInteger
        Number of subdicisions when using _Instances_ or _Mesh_ to represent atoms
    scale : InputFloat
        Scale the `vdw_radii` of the atom when setting the radius of the spheres
    shade_smooth : InputBoolean
        Apply smooth shading when using _Instances_ or _Mesh_
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this style to, discarding unselected points
    i.sphere : MenuSocket
        Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
    i.quality : IntegerSocket
        Number of subdicisions when using _Instances_ or _Mesh_ to represent atoms
    i.scale : FloatSocket
        Scale the `vdw_radii` of the atom when setting the radius of the spheres
    i.shade_smooth : BooleanSocket
        Apply smooth shading when using _Instances_ or _Mesh_
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        The generated geometry for the style node group
    """

    _name = "Style Spheres"
    _asset_name = "Style Spheres"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "node_tool_idname": "geometry.style_spheres",
        "show_modifier_manage_panel": False,
        "is_modifier": True,
    }

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this style to, discarding unselected points"""
        sphere: MenuSocket
        """Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles."""
        quality: IntegerSocket
        """Number of subdicisions when using _Instances_ or _Mesh_ to represent atoms"""
        scale: FloatSocket
        """Scale the `vdw_radii` of the atom when setting the radius of the spheres"""
        shade_smooth: BooleanSocket
        """Apply smooth shading when using _Instances_ or _Mesh_"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """The generated geometry for the style node group"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        sphere: InputMenu | Literal["Point", "Instance", "Mesh"] = "Point",
        quality: InputInteger = 2,
        scale: InputFloat = 0.8,
        shade_smooth: InputBoolean = True,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Sphere": sphere,
                "Quality": quality,
                "Scale": scale,
                "Shade Smooth": shade_smooth,
                "Material": material,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this style to, discarding unselected points",
            hide_value=True,
        )
        sphere = tree.inputs.menu(
            "Sphere",
            description="Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.",
            expanded=True,
            optional_label=True,
        )
        quality = tree.inputs.integer(
            "Quality",
            2,
            description="Number of subdicisions when using _Instances_ or _Mesh_ to represent atoms",
            min_value=0,
            max_value=5,
        )
        scale = tree.inputs.float(
            "Scale",
            0.8,
            description="Scale the `vdw_radii` of the atom when setting the radius of the spheres",
            min_value=0.0,
            max_value=2.0,
        )
        with tree.inputs.panel("Material", default_closed=True):
            shade_smooth = tree.inputs.boolean(
                "Shade Smooth",
                True,
                description="Apply smooth shading when using _Instances_ or _Mesh_",
            )
            material = tree.inputs.material(
                "Material",
                description="Material to apply to the resulting geometry",
                optional_label=True,
            )
        geometry = tree.outputs.geometry(
            "Geometry", description="The generated geometry for the style node group"
        )

        closure_zone = g.ClosureZone()
        atoms_1 = closure_zone.inputs.geometry("Atoms")
        geometry_1 = closure_zone.outputs.geometry("Geometry")
        group = MN_utils_style_spheres_icosphere(
            atoms=atoms_1,
            selection=selection,
            scale=scale,
            subdivisions=quality,
            shade_smooth=shade_smooth,
            material=material,
        )
        group_1 = MN_utils_style_spheres_points(
            atoms=atoms_1, selection=selection, scale=scale, material=material
        )
        menu_switch = g.MenuSwitch.geometry(
            sphere,
            {
                "Point": group_1,
                "Instance": group,
                "Mesh": g.RealizeInstances(
                    geometry=group, realize_to_point_domain=True
                ),
            },
        )
        menu_switch >> geometry_1
        EvaluateOnAtoms(geometry=atoms, closure=closure_zone.closure) >> geometry

        sphere.default_value = "Point"


ASSET = StyleSpheres

ASSET_METADATA = {
    "catalog_id": "541e6649-2ea6-4225-b1ee-5c0da6f5f1f6",
}
