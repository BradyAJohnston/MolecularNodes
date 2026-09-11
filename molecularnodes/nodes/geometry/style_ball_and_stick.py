# Node-group asset "Style Ball and Stick" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from ._shared.mn_utils_style_sticks import MN_utils_style_sticks
from .evaluate_on_atoms import EvaluateOnAtoms
from .find_bonds import FindBonds
from .style_spheres import StyleSpheres


class StyleBallAndStick(AssetGeometryGroup):
    """
    Style Ball and Stick

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this style to, discarding unselected points
    quality : InputInteger
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    sphere : InputMenu | Literal["Point", "Instance", "Mesh"]
        Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
    scale : InputFloat
        Scale the `vdw_radii` attribute before setting the radius for the spheres
    bond_split : InputMenu | Literal["Single", "Double"]
        Bond Split
    bond_scale : InputFloat
        Set the radius for the generated bonds in Angstroms
    bond_find : InputBoolean
        Find possible bonds for the selected atoms based on a distance search. Unselected atoms maintain any bonds they already have. Bonds that are found are all treated as single bonds
    bond_find_scale : InputFloat
        Scale the VDW radii of the atoms when searching for bonds
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
    i.quality : IntegerSocket
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    i.sphere : MenuSocket
        Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
    i.scale : FloatSocket
        Scale the `vdw_radii` attribute before setting the radius for the spheres
    i.bond_split : MenuSocket
        Bond Split
    i.bond_scale : FloatSocket
        Set the radius for the generated bonds in Angstroms
    i.bond_find : BooleanSocket
        Find possible bonds for the selected atoms based on a distance search. Unselected atoms maintain any bonds they already have. Bonds that are found are all treated as single bonds
    i.bond_find_scale : FloatSocket
        Scale the VDW radii of the atoms when searching for bonds
    i.shade_smooth : BooleanSocket
        Apply smooth shading to the created geometry
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        The generated geometry for the style node group
    """

    _name = "Style Ball and Stick"
    _asset_name = "Style Ball and Stick"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "node_tool_idname": "geometry.style_ball_and_stick",
        "show_modifier_manage_panel": False,
        "is_modifier": True,
    }

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this style to, discarding unselected points"""
        quality: IntegerSocket
        """A lower value results in less geometry, with a higher value meaning better looking but more dense geometry"""
        sphere: MenuSocket
        """Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles."""
        scale: FloatSocket
        """Scale the `vdw_radii` attribute before setting the radius for the spheres"""
        bond_split: MenuSocket
        """Bond Split"""
        bond_scale: FloatSocket
        """Set the radius for the generated bonds in Angstroms"""
        bond_find: BooleanSocket
        """Find possible bonds for the selected atoms based on a distance search. Unselected atoms maintain any bonds they already have. Bonds that are found are all treated as single bonds"""
        bond_find_scale: FloatSocket
        """Scale the VDW radii of the atoms when searching for bonds"""
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
        quality: InputInteger = 2,
        sphere: InputMenu | Literal["Point", "Instance", "Mesh"] = "Instance",
        scale: InputFloat = 0.3,
        bond_split: InputMenu | Literal["Single", "Double"] = "Double",
        bond_scale: InputFloat = 0.3,
        bond_find: InputBoolean = False,
        bond_find_scale: InputFloat = 1.0,
        shade_smooth: InputBoolean = True,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Quality": quality,
                "Sphere": sphere,
                "Scale": scale,
                "Bond Split": bond_split,
                "Bond Scale": bond_scale,
                "Bond Find": bond_find,
                "Bond Find Scale": bond_find_scale,
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
        quality = tree.inputs.integer(
            "Quality",
            2,
            description="A lower value results in less geometry, with a higher value meaning better looking but more dense geometry",
            min_value=0,
            max_value=8,
            structure_type="SINGLE",
            force_non_field=True,
        )
        with tree.inputs.panel("Sphere", default_closed=True):
            sphere = tree.inputs.menu(
                "Sphere",
                description="Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.",
                expanded=True,
                optional_label=True,
            )
            scale = tree.inputs.float(
                "Scale",
                0.3,
                description="Scale the `vdw_radii` attribute before setting the radius for the spheres",
                min_value=0.0,
                max_value=2.0,
            )
        with tree.inputs.panel("Bond", default_closed=True):
            bond_split = tree.inputs.menu(
                "Bond Split", expanded=True, optional_label=True
            )
            bond_scale = tree.inputs.float(
                "Bond Scale",
                0.3,
                description="Set the radius for the generated bonds in Angstroms",
                min_value=0.0,
                max_value=1.0,
            )
            with tree.inputs.panel(
                "Bond Find",
                description="Search for nearby atoms to create bonds based on VDW radii. Removes existing bonds.",
                default_closed=True,
            ):
                bond_find = tree.inputs.boolean(
                    "Bond Find",
                    False,
                    description="Find possible bonds for the selected atoms based on a distance search. Unselected atoms maintain any bonds they already have. Bonds that are found are all treated as single bonds",
                    is_panel_toggle=True,
                )
                bond_find_scale = tree.inputs.float(
                    "Bond Find Scale",
                    1.0,
                    description="Scale the VDW radii of the atoms when searching for bonds",
                    min_value=0.0,
                    max_value=2.0,
                    subtype="FACTOR",
                )
        with tree.inputs.panel("Material"):
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
        group = StyleSpheres(
            atoms=separate_geometry.o.selection,
            sphere=sphere,
            quality=quality,
            scale=scale,
            shade_smooth=shade_smooth,
            material=material,
        )
        switch = bond_find.switch.geometry(
            separate_geometry.o.selection,
            FindBonds(atoms=separate_geometry.o.selection, scale=bond_find_scale),
        )
        group_1 = MN_utils_style_sticks(
            atoms=switch,
            radius=bond_scale,
            resolution=quality * 6,
            scale_extra_bond_radius=0.45,
            extra_bond_offset=scale * 0.8,
            menu=bond_split,
            shade_smooth=shade_smooth,
            material=material,
        )
        (
            g.SetMaterial(
                geometry=g.JoinGeometry(geometry=(group, group_1)), material=material
            )
            >> geometry_1
        )
        EvaluateOnAtoms(geometry=atoms, closure=closure_zone.closure) >> geometry

        sphere.default_value = "Instance"
        bond_split.default_value = "Double"


ASSET = StyleBallAndStick

ASSET_METADATA = {
    "catalog_id": "541e6649-2ea6-4225-b1ee-5c0da6f5f1f6",
}

DATABLOCK_DEPENDENCIES = {
    "materials": ("MN Default.old",),
}
