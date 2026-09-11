# Node-group asset "Style Sticks" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from ._shared.mn_units import MNUnits
from ._shared.mn_utils_style_sticks import MN_utils_style_sticks
from .evaluate_on_atoms import EvaluateOnAtoms
from .style_spheres import StyleSpheres


class StyleSticks(AssetGeometryGroup):
    """
    Style Sticks

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this style to, discarding unselected points
    sphere : InputMenu | Literal["Point", "Instance", "Mesh"]
        Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
    quality : InputInteger
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    scale : InputFloat
        Radius of the sticks in Angstroms
    shade_smooth : InputBoolean
        Apply smooth shading to the created geometry
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
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    i.scale : FloatSocket
        Radius of the sticks in Angstroms
    i.shade_smooth : BooleanSocket
        Apply smooth shading to the created geometry
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        The generated geometry for the style node group
    """

    _name = "Style Sticks"
    _asset_name = "Style Sticks"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "node_tool_idname": "geometry.style_sticks",
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
        """A lower value results in less geometry, with a higher value meaning better looking but more dense geometry"""
        scale: FloatSocket
        """Radius of the sticks in Angstroms"""
        shade_smooth: BooleanSocket
        """Apply smooth shading to the created geometry"""
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
        sphere: InputMenu | Literal["Point", "Instance", "Mesh"] = "Instance",
        quality: InputInteger = 3,
        scale: InputFloat = 0.2,
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
            3,
            description="A lower value results in less geometry, with a higher value meaning better looking but more dense geometry",
            min_value=0,
            max_value=5,
        )
        scale = tree.inputs.float(
            "Scale",
            0.2,
            description="Radius of the sticks in Angstroms",
            min_value=0.0,
            max_value=1.0,
        )
        with tree.inputs.panel("Material", default_closed=True):
            shade_smooth = tree.inputs.boolean(
                "Shade Smooth",
                True,
                description="Apply smooth shading to the created geometry",
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
        separate_geometry = atoms_1 >> g.SeparateGeometry.point(selection=selection)
        group = MN_utils_style_sticks(
            atoms=separate_geometry.o.selection,
            radius=scale,
            resolution=g.Math.multiply(quality, g.Integer(integer=8)),
            scale_extra_bond_radius=0.37,
            extra_bond_offset=scale / 2.0,
            shade_smooth=shade_smooth,
            material=material,
        )
        store_named_attribute = g.StoreNamedAttribute.point.float(
            separate_geometry.o.selection,
            name="vdw_radii",
            value=MNUnits(value=1.0).o.angstrom,
        )
        group_1 = StyleSpheres(
            atoms=store_named_attribute,
            sphere=sphere,
            quality=quality,
            scale=scale,
            shade_smooth=shade_smooth,
            material=material,
        )
        g.JoinGeometry(geometry=(group_1, group)) >> geometry_1
        EvaluateOnAtoms(geometry=atoms, closure=closure_zone.closure) >> geometry

        sphere.default_value = "Instance"


ASSET = StyleSticks

ASSET_METADATA = {
    "catalog_id": "541e6649-2ea6-4225-b1ee-5c0da6f5f1f6",
}
