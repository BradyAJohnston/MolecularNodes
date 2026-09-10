# Node group '.MN_utils_style_spheres_icosphere' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MaterialSocket,
    SocketAccessor,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputMaterial,
)
from ..vdw_radii import VDWRadii
from .set_instancer import SetInstancer


class MN_utils_style_spheres_icosphere(CustomGeometryGroup):
    """
    .MN_utils_style_spheres_icosphere

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this node to
    scale : InputFloat
        Scale the VDW radii of the atoms.
    subdivisions : InputInteger
        Subdivisions
    shade_smooth : InputBoolean
        Apply smooth shading to the created geometry
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.scale : FloatSocket
        Scale the VDW radii of the atoms.
    i.subdivisions : IntegerSocket
        Subdivisions
    i.shade_smooth : BooleanSocket
        Apply smooth shading to the created geometry
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.instances : GeometrySocket
        Instances
    """

    _name = ".MN_utils_style_spheres_icosphere"
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "node_tool_idname": "geometry._mn_utils_style_spheres_icosphere"
    }

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        scale: FloatSocket
        """Scale the VDW radii of the atoms."""
        subdivisions: IntegerSocket
        """Subdivisions"""
        shade_smooth: BooleanSocket
        """Apply smooth shading to the created geometry"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""

    class _Outputs(SocketAccessor):
        instances: GeometrySocket
        """Instances"""

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
        subdivisions: InputInteger = 2,
        shade_smooth: InputBoolean = True,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Scale": scale,
                "Subdivisions": subdivisions,
                "Shade Smooth": shade_smooth,
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
        scale = tree.inputs.float(
            "Scale",
            0.8,
            description="Scale the VDW radii of the atoms.",
            min_value=0.0,
            max_value=10_000.0,
        )
        subdivisions = tree.inputs.integer("Subdivisions", 2, min_value=0, max_value=5)
        shade_smooth = tree.inputs.boolean(
            "Shade Smooth",
            True,
            description="Apply smooth shading to the created geometry",
        )
        material = tree.inputs.material(
            "Material", description="Material to apply to the resulting geometry"
        )
        instances = tree.outputs.geometry("Instances")

        switch = g.Switch.geometry(
            subdivisions,
            g.TransformGeometry(
                geometry=g.Cube(), rotation=(math.pi / 4, math.pi / 4, math.pi / 4)
            ),
            g.IcoSphere(subdivisions=subdivisions),
        )
        (
            SetInstancer(geometry=atoms)
            >> g.InstanceOnPoints(
                selection=selection,
                instance=switch,
                scale=VDWRadii().o.vdw_radii * scale,
            )
            >> g.SetShadeSmooth.face(shade_smooth=shade_smooth)
            >> g.SetMaterial(material=material)
            >> instances
        )
