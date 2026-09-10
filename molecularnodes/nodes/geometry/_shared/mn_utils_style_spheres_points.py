# Node group '.MN_utils_style_spheres_points' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    MaterialSocket,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry, InputMaterial
from ..vdw_radii import VDWRadii


class MN_utils_style_spheres_points(CustomGeometryGroup):
    """
    .MN_utils_style_spheres_points

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this node to
    scale : InputFloat
        Scale
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.scale : FloatSocket
        Scale
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.point_cloud : GeometrySocket
        Point Cloud
    """

    _name = ".MN_utils_style_spheres_points"
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry._mn_utils_style_spheres_points"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        scale: FloatSocket
        """Scale"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""

    class _Outputs(SocketAccessor):
        point_cloud: GeometrySocket
        """Point Cloud"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        scale: InputFloat = 0.8,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Scale": scale,
                "Material": material,
            }
        )

    def _build_group(self, tree):
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        scale = tree.inputs.float("Scale", 0.8, min_value=0.0, max_value=10_000.0)
        material = tree.inputs.material(
            "Material", description="Material to apply to the resulting geometry"
        )
        point_cloud = tree.outputs.geometry("Point Cloud")

        (
            atoms
            >> g.MeshToPoints(selection=selection, radius=scale * VDWRadii())
            >> g.SetMaterial(material=material)
            >> point_cloud
        )
