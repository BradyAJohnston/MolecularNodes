# Node-group asset "Style Surface" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    CustomGeometryGroup,
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
from ._shared.mn_constants_atom_name_nucleic import MN_constants_atom_name_nucleic
from ._shared.mn_world_scale import MN_world_scale
from .angstrom_to_world import AngstromToWorld
from .atom_name import AtomName
from .color import Color
from .edge_length import EdgeLength
from .evaluate_on_atoms import EvaluateOnAtoms
from .evaluate_ordered_bundles import EvaluateOrderedBundles
from .evaluate_per_group import EvaluatePerGroup
from .evluate_while_planar import EvluateWhilePlanar
from .is_alpha_carbon import IsAlphaCarbon
from .mn_typed_bundles import MNTypedBundles
from .set_color import SetColor
from .vdw_radii import VDWRadii


class MN_utils_style_surface_new(CustomGeometryGroup):
    _name = ".MN_utils_style_surface_new"
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry._mn_utils_style_surface_new"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        shade_smooth = tree.inputs.boolean(
            "Shade Smooth",
            True,
            description="Apply smooth shading to the created geometry",
        )
        material = tree.inputs.material(
            "Material", description="Material to apply to the resulting geometry"
        )
        _relaxation_steps = tree.inputs.integer("Relaxation Steps", 30, min_value=0)
        grid_create = tree.inputs.closure("Grid Create")
        mesh_processing = tree.inputs.bundle("Mesh Processing")
        _threshold = tree.inputs.float(
            "Threshold",
            0.0,
            description="Values larger than the threshold are inside the generated mesh",
        )
        geometry = tree.outputs.geometry("Geometry")

        evaluate_closure = g.EvaluateClosure(grid_create, define_signature=True)
        evaluate_closure.inputs.geometry("Atoms", atoms)
        grid = evaluate_closure.outputs.geometry("Grid", structure_type="SINGLE")
        _field_to_list = g.FieldToList(count=g.ListLength.string(""), items={})
        (
            EvaluateOrderedBundles(geometry=grid, bundles=mesh_processing)
            >> g.SetShadeSmooth.face(shade_smooth=shade_smooth)
            >> g.SetMaterial(material=material)
            >> geometry
        )
        _sort_list = g.SortList.float()


class Surface_compute_density_from_points(CustomGeometryGroup):
    _name = ".surface_compute_density_from_points"
    _tree_properties = {
        "node_tool_idname": "geometry._surface_compute_density_from_points"
    }

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        scale_radius = tree.inputs.float(
            "Scale Radius", 1.0, min_value=-10_000.0, max_value=10_000.0
        )
        probe_size = tree.inputs.float(
            "Probe Size", 0.0, min_value=0.0, max_value=10_000.0
        )
        result = tree.outputs.boolean("Result")
        distance = tree.outputs.float("Distance")

        position = g.Position()
        sample_nearest = g.SampleNearest.point(atoms)
        sample_index = g.SampleIndex(
            geometry=atoms,
            value=position,
            index=sample_nearest,
            data_type="FLOAT_VECTOR",
        )
        sample_index_1 = g.SampleIndex(
            geometry=atoms,
            value=VDWRadii().o.vdw_radii * scale_radius,
            index=sample_nearest,
        )
        math_1 = probe_size + sample_index_1 - sample_index.o.value.distance(position)
        (math_1 > 0.0) >> result

        math_1 >> distance


class Utils_bounding_box(CustomGeometryGroup):
    _name = ".utils_bounding_box"
    _tree_properties = {"node_tool_idname": "geometry._utils_bounding_box"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        subdivisions = tree.inputs.float(
            "Subdivisions", 16.7, min_value=-10_000.0, max_value=10_000.0
        )
        min = tree.outputs.vector("Min")
        max = tree.outputs.vector("Max")
        x = tree.outputs.integer("X")
        y = tree.outputs.integer("Y")
        z = tree.outputs.integer("Z")

        group = MN_world_scale()
        bounding_box = g.BoundingBox(geometry=geometry)
        math_1 = group.o.world_scale * 2.0
        vector_math = g.VectorMath.snap(bounding_box.o.min, group).o.vector - math_1
        vector_math_1 = g.VectorMath.snap(bounding_box.o.max, group).o.vector + math_1
        vector = (vector_math_1 - vector_math) * subdivisions
        vector.x.max(2.0) >> x
        vector.y.max(2.0) >> y
        vector.z.max(2.0) >> z

        vector_math >> min
        vector_math_1 >> max


class MN_surface_smooth_bumps(CustomGeometryGroup):
    _name = ".MN_surface_smooth_bumps"
    _tree_properties = {"node_tool_idname": "geometry._mn_surface_smooth_bumps"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        geometry_1 = tree.outputs.geometry("Geometry")

        with g.Frame("Smoothen out weird bumps from meshing"):
            edge_vertices = g.EdgeVertices()
            compare = g.Compare.integer.equal(g.VertexNeighbors().o.vertex_count, 3)
            boolean_math = compare.o.result.point.at(
                edge_vertices.o.vertex_index_1
            ) & compare.o.result.point.at(edge_vertices.o.vertex_index_2)
            face_group_boundaries = g.FaceGroupBoundaries(
                face_set=g.EdgesToFaceGroups(boundary_edges=boolean_math)
            )
        (
            geometry
            >> g.SetPosition(
                selection=face_group_boundaries,
                position=g.BlurAttribute.vector(g.Position(), 4),
            )
            >> geometry_1
        )


class MN_utils_style_surface_sdf(CustomGeometryGroup):
    _name = ".MN_utils_style_surface_sdf"
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry._mn_utils_style_surface_new"}

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        quality = tree.inputs.integer("Quality", 12, min_value=1, max_value=15)
        scale_radii = tree.inputs.float(
            "Scale Radii", 1.0, min_value=0.0, max_value=10.0
        )
        probe_size = tree.inputs.float(
            "Probe Size", 0.6, min_value=0.0, max_value=10_000.0
        )
        color_source = tree.inputs.menu("Color Source", optional_label=True)
        color_blur = tree.inputs.integer(
            "Color Blur",
            1,
            description="Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating",
            min_value=0,
            max_value=20,
        )
        shade_smooth = tree.inputs.boolean(
            "Shade Smooth",
            True,
            description="Apply smooth shading to the created geometry",
        )
        material = tree.inputs.material(
            "Material", description="Material to apply to the resulting geometry"
        )
        relaxation_steps = tree.inputs.integer("Relaxation Steps", 30, min_value=0)
        surface_geometry = tree.outputs.geometry("Surface Geometry")

        compare = g.Compare.integer.equal(
            MN_constants_atom_name_nucleic().o.side_chain_joint_carbon, AtomName()
        )
        separate_geometry = g.SeparateGeometry.point(atoms, selection)
        group = Surface_compute_density_from_points(
            Atoms=separate_geometry.o.selection,
            **{"Scale Radius": scale_radii, "Probe Size": probe_size},
        )
        group_1 = Utils_bounding_box(
            Geometry=separate_geometry.o.selection, Subdivisions=quality * 5.0
        )
        separate_geometry_1 = g.SeparateGeometry.point(
            separate_geometry.o.selection, IsAlphaCarbon().o.selection | compare
        )
        compare_1 = g.Compare.integer.equal(
            g.DomainSize(geometry=separate_geometry_1.o.selection).o.point_count, 0
        )
        switch = compare_1.o.result.switch.geometry(
            separate_geometry_1.o.selection, separate_geometry_1.o.inverted
        )
        menu_switch = g.MenuSwitch.geometry(
            color_source,
            {"Alpha Carbon": switch, "Nearest": separate_geometry.o.selection},
        )
        with g.Frame("Generate Surface from Measurements"):
            volume_cube = g.VolumeCube(
                density=group.o.result,
                min=group_1.o.min,
                max=group_1.o.max,
                resolution_x=group_1.o.x,
                resolution_y=group_1.o.y,
                resolution_z=group_1.o.z,
            )
            volume_to_mesh = g.VolumeToMesh(
                volume=volume_cube, voxel_size=0.01, threshold=0.1
            )
        with g.Frame("Pull in surface towards atoms"):
            position = g.Position()
            capture = g.CaptureAttribute.point(geometry=volume_to_mesh)
            value = capture.items.integer(
                "Value", g.SampleNearest.point(separate_geometry.o.selection)
            )
            sample_index = g.SampleIndex(
                geometry=separate_geometry.o.selection,
                value=position,
                index=value.output,
                data_type="FLOAT_VECTOR",
            )
            sample_index_1 = g.SampleIndex(
                geometry=separate_geometry.o.selection,
                value=scale_radii * VDWRadii(),
                index=value.output,
            )
            set_position = capture.o.geometry >> g.SetPosition(
                position=sample_index.o.value
                + (position.o.position - sample_index).normalize() * sample_index_1
            )
        with g.Frame("smoothing of tightened surface"):
            blur_attribute = g.BlurAttribute.float(
                g.EdgeAngle().o.signed_angle.map_range(-0.2, -1.0).face.evaluate(), 2
            )
            position_1 = g.Position()
            capture_1 = g.CaptureAttribute.point(geometry=set_position)
            value_1 = capture_1.items.vector("Value", position_1)
            vector_math = g.Normal(legacy_corner_normals=True).o.normal * (
                value_1.output.distance(position_1) * -1.8
            )
            set_position_1 = (
                capture_1.o.geometry
                >> g.SetPosition(position=g.BlurAttribute.vector(position_1, 2))
                >> g.SetPosition(offset=vector_math)
                >> g.SetPosition(position=g.BlurAttribute.vector(g.Position(), 2))
                >> g.SetPosition(
                    position=g.BlurAttribute.vector(
                        g.Position(), relaxation_steps, blur_attribute
                    )
                )
            )
            group_2 = MN_surface_smooth_bumps(Geometry=set_position_1)
        triangulate = group_2 >> g.Triangulate(quad_method="Beauty")
        with g.Frame("Sample colors of nearest atom"):
            sample_index_2 = g.SampleIndex(
                geometry=menu_switch,
                value=Color(),
                index=g.SampleNearest.point(menu_switch),
                data_type="FLOAT_COLOR",
            )
            group_3 = SetColor(
                atoms=triangulate,
                color=g.BlurAttribute.color(sample_index_2, color_blur),
            )
        (
            group_3
            >> g.SetShadeSmooth.face(shade_smooth=shade_smooth)
            >> g.SetMaterial(material=material)
            >> surface_geometry
        )

        color_source.default_value = "Alpha Carbon"


class TriangulateMesh(CustomGeometryGroup):
    _name = ".Triangulate Mesh"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        step = tree.inputs.integer("step", 4)
        bundle = tree.outputs.bundle("Bundle")

        closure_zone = g.ClosureZone()
        geometry = closure_zone.inputs.geometry("Geometry")
        geometry_1 = closure_zone.outputs.geometry("Geometry")
        group = MN_surface_smooth_bumps(Geometry=geometry)
        group.node.mute = True
        triangulate = group >> g.Triangulate(quad_method="Beauty")
        triangulate.node.mute = True
        triangulate >> geometry_1
        (
            MNTypedBundles(
                closure=closure_zone.closure, step=step, path="Triangulate Mesh"
            )
            >> bundle
        )


class RelaxSurface(CustomGeometryGroup):
    _name = ".Relax Surface"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        relaxation_steps = tree.inputs.integer("Relaxation Steps", 30, min_value=0)
        step = tree.inputs.integer("step", 3)
        bundle = tree.outputs.bundle("Bundle")

        closure_zone = g.ClosureZone()
        geometry = closure_zone.inputs.geometry("Geometry")
        geometry_1 = closure_zone.outputs.geometry("Geometry")
        group = EdgeLength()
        field_min_max = g.FieldMinAndMax.point.float(group)
        _map_range = group.o.length.map_range(field_min_max.o.min, field_min_max.o.max)
        blur_attribute = g.BlurAttribute.vector(
            g.Position(),
            relaxation_steps,
            g.BlurAttribute.float(
                g.EdgeAngle().o.signed_angle.map_range(-0.2, 0.3, 1.0, 0.0)
            ),
        )
        (
            MN_surface_smooth_bumps(
                Geometry=geometry >> g.SetPosition(position=blur_attribute)
            )
            >> geometry_1
        )
        (
            MNTypedBundles(
                closure=closure_zone.closure, step=step, path="Relax Surface"
            )
            >> bundle
        )


class SampleColors(CustomGeometryGroup):
    _name = ".Sample Colors"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        color_source = tree.inputs.menu("Color Source", optional_label=True)
        atoms = tree.inputs.geometry("Atoms")
        step = tree.inputs.integer("step", 1)
        blur = tree.inputs.integer(
            "Blur",
            1,
            description="How many times to blur the values for all elements",
            min_value=0,
        )
        bundle = tree.outputs.bundle("Bundle")

        closure_zone = g.ClosureZone()
        geometry = closure_zone.inputs.geometry("Geometry")
        geometry_1 = closure_zone.outputs.geometry("Geometry")
        compare = g.Compare.integer.equal(
            MN_constants_atom_name_nucleic().o.side_chain_joint_carbon, AtomName()
        )
        with g.Frame("Don't transfer these attributes"):
            field_to_list = g.FieldToList(
                count=2,
                items={
                    "Attribute Names": g.IndexSwitch.string(
                        g.Index(), ("bond_type", "position")
                    )
                },
            )
        separate_geometry = g.SeparateGeometry.point(
            atoms, IsAlphaCarbon().o.selection | compare
        )
        compare_1 = g.Compare.integer.equal(
            g.DomainSize(geometry=separate_geometry.o.selection).o.point_count, 0
        )
        switch = compare_1.o.result.switch.geometry(
            separate_geometry.o.selection, separate_geometry.o.inverted
        )
        menu_switch = g.MenuSwitch.geometry(
            color_source, {"Alpha Carbon": switch, "Nearest": atoms}
        )
        capture = g.CaptureAttribute.point(geometry=geometry)
        index = capture.items.integer("Index", g.SampleNearest.point(menu_switch))
        transfer_attributes = capture.o.geometry >> g.TransferAttributes(
            target_point_id=index.output,
            source=menu_switch,
            attribute_names=field_to_list,
            pattern_mode="Exact",
            exclude_names=True,
        )
        (
            SetColor(
                atoms=transfer_attributes, color=g.BlurAttribute.color(Color(), blur)
            )
            >> geometry_1
        )
        (
            MNTypedBundles(
                closure=closure_zone.closure, step=step, path="Sample Colors"
            )
            >> bundle
        )

        color_source.default_value = "Alpha Carbon"


class SurfaceToRadius(CustomGeometryGroup):
    _name = ".Surface to Radius"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        scale_radii = tree.inputs.float(
            "Scale Radii", 1.47, min_value=0.0, max_value=10.0
        )
        atoms = tree.inputs.geometry("Atoms")
        step = tree.inputs.integer("step", 0)
        bundle = tree.outputs.bundle("Bundle")

        closure_zone = g.ClosureZone()
        geometry = closure_zone.inputs.geometry("Geometry")
        geometry_1 = closure_zone.outputs.geometry("Geometry")
        position = g.Position()
        capture = g.CaptureAttribute.point(geometry=geometry)
        value = capture.items.integer("Value", g.SampleNearest.point(atoms))
        sample_index = g.SampleIndex(
            geometry=atoms, value=position, index=value.output, data_type="FLOAT_VECTOR"
        )
        capture_1 = g.CaptureAttribute.point(geometry=capture.o.geometry)
        value_001 = capture_1.items.vector("Value.001", sample_index)
        vector_math = (
            position.o.position - value_001.output
        ).normalize() * g.SampleIndex(
            geometry=atoms, value=scale_radii * VDWRadii(), index=value.output
        )
        mix = g.Mix(
            a_vector=value_001.output + vector_math,
            b_vector=g.Position(),
            factor_float=0.0,
            data_type="VECTOR",
            clamp_factor=True,
        )
        (
            capture_1.o.geometry
            >> g.SetPosition(position=mix.o.result_vector)
            >> geometry_1
        )
        (
            MNTypedBundles(
                closure=closure_zone.closure, step=step, path="Surface to Radius"
            )
            >> bundle
        )


class StyleSurface(AssetGeometryGroup):
    """
    Style Surface

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this style to, discarding unselected points
    quality : InputInteger
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    scale : InputFloat
        Scale the VDW radii of the atoms when creating the surface
    relax : InputInteger
        Relax
    offset : InputFloat
        Object-space distance to offset the SDF surface
    fillet : InputInteger
        Number of iterations to apply the filter
    mean_width : InputInteger
        Filter kernel radius in voxels
    mean_iterations : InputInteger
        Number of iterations to apply the filter
    separate_by : InputMenu | Literal["chain_id", "Group ID"]
        Separate By
    group_id : InputInteger
        Group ID
    color_source : InputMenu | Literal["Alpha Carbon", "Nearest"]
        Color Source
    color_blur : InputInteger
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
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
    i.scale : FloatSocket
        Scale the VDW radii of the atoms when creating the surface
    i.relax : IntegerSocket
        Relax
    i.offset : FloatSocket
        Object-space distance to offset the SDF surface
    i.fillet : IntegerSocket
        Number of iterations to apply the filter
    i.mean_width : IntegerSocket
        Filter kernel radius in voxels
    i.mean_iterations : IntegerSocket
        Number of iterations to apply the filter
    i.separate_by : MenuSocket
        Separate By
    i.group_id : IntegerSocket
        Group ID
    i.color_source : MenuSocket
        Color Source
    i.color_blur : IntegerSocket
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    i.shade_smooth : BooleanSocket
        Apply smooth shading to the created geometry
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        The generated geometry for the style node group
    """

    _name = "Style Surface"
    _asset_name = "Style Surface"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "node_tool_idname": "geometry.style_surface",
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
        scale: FloatSocket
        """Scale the VDW radii of the atoms when creating the surface"""
        relax: IntegerSocket
        """Relax"""
        offset: FloatSocket
        """Object-space distance to offset the SDF surface"""
        fillet: IntegerSocket
        """Number of iterations to apply the filter"""
        mean_width: IntegerSocket
        """Filter kernel radius in voxels"""
        mean_iterations: IntegerSocket
        """Number of iterations to apply the filter"""
        separate_by: MenuSocket
        """Separate By"""
        group_id: IntegerSocket
        """Group ID"""
        color_source: MenuSocket
        """Color Source"""
        color_blur: IntegerSocket
        """Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating"""
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
        quality: InputInteger = 3,
        scale: InputFloat = 1.25,
        relax: InputInteger = 5,
        offset: InputFloat = 0.15,
        fillet: InputInteger = 0,
        mean_width: InputInteger = 1,
        mean_iterations: InputInteger = 1,
        separate_by: InputMenu | Literal["chain_id", "Group ID"] = "chain_id",
        group_id: InputInteger = 0,
        color_source: InputMenu | Literal["Alpha Carbon", "Nearest"] = "Alpha Carbon",
        color_blur: InputInteger = 2,
        shade_smooth: InputBoolean = True,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Quality": quality,
                "Scale": scale,
                "Relax": relax,
                "Offset": offset,
                "Fillet": fillet,
                "Mean Width": mean_width,
                "Mean Iterations": mean_iterations,
                "Separate By": separate_by,
                "Group ID": group_id,
                "Color Source": color_source,
                "Color Blur": color_blur,
                "Shade Smooth": shade_smooth,
                "Material": material,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms",
            description="Atomic geometry that contains vertices and edges",
            hide_value=True,
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this style to, discarding unselected points",
            hide_value=True,
        )
        quality = tree.inputs.integer(
            "Quality",
            3,
            description="A lower value results in less geometry, with a higher value meaning better looking but more dense geometry",
            min_value=0,
            max_value=6,
        )
        with tree.inputs.panel("Surface", default_closed=True):
            scale = tree.inputs.float(
                "Scale",
                1.25,
                description="Scale the VDW radii of the atoms when creating the surface",
                min_value=0.0,
                max_value=10.0,
            )
            relax = tree.inputs.integer("Relax", 5, min_value=0)
            with tree.inputs.panel("SDF", default_closed=True):
                offset = tree.inputs.float(
                    "Offset",
                    0.15,
                    description="Object-space distance to offset the SDF surface",
                    min_value=0.0,
                    subtype="DISTANCE",
                )
                fillet = tree.inputs.integer(
                    "Fillet",
                    0,
                    description="Number of iterations to apply the filter",
                    min_value=0,
                )
                with tree.inputs.panel("Mean", default_closed=True):
                    mean_width = tree.inputs.integer(
                        "Mean Width",
                        1,
                        description="Filter kernel radius in voxels",
                        min_value=0,
                    )
                    mean_iterations = tree.inputs.integer(
                        "Mean Iterations",
                        1,
                        description="Number of iterations to apply the filter",
                        min_value=0,
                    )
            with tree.inputs.panel("Separate", default_closed=True):
                separate_by = tree.inputs.menu(
                    "Separate By", expanded=True, optional_label=True
                )
                group_id = tree.inputs.integer("Group ID", 0, hide_value=True)
        with tree.inputs.panel("Material", default_closed=True):
            color_source = tree.inputs.menu(
                "Color Source", expanded=True, optional_label=True
            )
            color_blur = tree.inputs.integer(
                "Color Blur",
                2,
                description="Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating",
                min_value=0,
                max_value=20,
            )
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

        _group = MN_utils_style_surface_new()
        with g.Frame("Create grid and turn into mesh"):
            with g.Frame("Old-style 'manual SDF creation"):
                closure_zone = g.ClosureZone()
                atoms_1 = closure_zone.inputs.geometry("Atoms")
                geometry_1 = closure_zone.outputs.geometry(
                    "Geometry", structure_type="SINGLE"
                )
                group_1 = Utils_bounding_box(
                    Geometry=atoms_1, Subdivisions=quality * 5.0
                )
                cube_grid_topology = g.CubeGridTopology(
                    bounds_min=group_1.o.min,
                    bounds_max=group_1.o.max,
                    resolution_x=group_1.o.x,
                    resolution_y=group_1.o.y,
                    resolution_z=group_1.o.z,
                )
                field_to_grid = g.FieldToGrid.boolean(topology=cube_grid_topology)
                density = field_to_grid.items.float(
                    "Density",
                    Surface_compute_density_from_points(
                        Atoms=atoms_1, **{"Scale Radius": scale}
                    ).o.distance,
                )
                grid_mean = g.SetGridBackground.float(density.grid).o.grid.mean()
                grid_mean.node.mute = True
                subdivision_surface = grid_mean.median().to_mesh(
                    AngstromToWorld(angstrom=4.3)
                ) >> g.SubdivisionSurface(limit_surface=False, quality=1)
                subdivision_surface >> geometry_1
            with g.Frame("New SDF Node"):
                closure_zone_1 = g.ClosureZone()
                atoms_2 = closure_zone_1.inputs.geometry("Atoms")
                geometry_2 = closure_zone_1.outputs.geometry("Geometry")
                points_to_sdf_grid = atoms_2 >> g.PointsToSDFGrid(
                    radius=scale * VDWRadii(),
                    voxel_size=AngstromToWorld(
                        angstrom=2.0 / g.Switch.float(quality, 0.5, quality)
                    ),
                )
                sdf_grid_fillet = (
                    points_to_sdf_grid.o.sdf_grid.sdf_offset(offset)
                    .sdf_mean(mean_width, mean_iterations)
                    .sdf_fillet(fillet)
                )
                sdf_grid_fillet.to_mesh(0.0) >> geometry_2
        closure_zone_2 = g.ClosureZone()
        atoms_3 = closure_zone_2.inputs.geometry("Atoms")
        geometry_3 = closure_zone_2.outputs.geometry("Geometry")
        closure_zone_3 = g.ClosureZone()
        geometry_4 = closure_zone_3.inputs.geometry("Geometry")
        geometry_5 = closure_zone_3.outputs.geometry("Geometry")
        group_2 = MN_utils_style_surface_sdf(
            Atoms=geometry_4,
            Selection=selection,
            Quality=quality,
            **{
                "Scale Radii": scale,
                "Color Blur": color_blur,
                "Shade Smooth": shade_smooth,
            },
            Material=material,
        )
        group_2 >> geometry_5
        menu_switch = g.MenuSwitch.closure(
            "SDF", {"Current": closure_zone.closure, "SDF": closure_zone_1.closure}
        )
        with g.Frame("Mesh processing"):
            group_3 = TriangulateMesh()
            group_3.node.mute = True
            closure_zone_4 = g.ClosureZone()
            geometry_6 = closure_zone_4.inputs.geometry("Geometry")
            geometry_7 = closure_zone_4.outputs.geometry("Geometry")
            join_bundle = g.JoinBundle(
                bundle=(
                    SurfaceToRadius(**{"Scale Radii": scale}, Atoms=geometry_6),
                    SampleColors(
                        **{"Color Source": color_source},
                        Atoms=geometry_6,
                        Blur=color_blur,
                    ),
                    RelaxSurface(**{"Relaxation Steps": relax}),
                    group_3,
                )
            )
            evaluate_closure = g.EvaluateClosure(
                menu_switch.o.output, define_signature=True
            )
            evaluate_closure.inputs.geometry("Atoms", geometry_6)
            geometry_8 = evaluate_closure.outputs.geometry(
                "Geometry", structure_type="SINGLE"
            )
            set_material = (
                EvaluateOrderedBundles(geometry=geometry_8, bundles=join_bundle)
                >> g.SetShadeSmooth.face(shade_smooth=shade_smooth)
                >> g.SetMaterial(material=material)
            )
            set_material >> geometry_7
        with g.Frame('Create surface while points are as "flat" as possible'):
            closure_zone_5 = g.ClosureZone()
            geometry_9 = closure_zone_5.inputs.geometry("Geometry")
            group_id_1 = closure_zone_5.inputs.integer("group_id")
            geometry_10 = closure_zone_5.outputs.geometry("Geometry")
            _string = g.String(
                string="We get better performance if we first orient the structure better inside of a bounding box for more efficient use of grid space & voxels!"
            )
            store_named_attribute = EvluateWhilePlanar(
                geometry=geometry_9, closure=closure_zone_4.closure
            ) >> g.StoreNamedAttribute.point.integer(name="chain_id", value=group_id_1)
            store_named_attribute >> geometry_10
        group_4 = EvaluatePerGroup(
            geometry=atoms_3,
            closure=closure_zone_5.closure,
            group=separate_by,
            group_id=group_id,
        )
        group_4 >> geometry_3
        (
            EvaluateOnAtoms(
                geometry=atoms, selection=selection, closure=closure_zone_2.closure
            )
            >> geometry
        )

        separate_by.default_value = "chain_id"
        color_source.default_value = "Alpha Carbon"


ASSET = StyleSurface

ASSET_METADATA = {
    "catalog_id": "541e6649-2ea6-4225-b1ee-5c0da6f5f1f6",
}
