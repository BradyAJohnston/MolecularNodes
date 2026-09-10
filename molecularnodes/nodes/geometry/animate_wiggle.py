# Node-group asset 'Animate Wiggle' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry
from ._shared.mn_select_res_name_peptide import MN_select_res_name_peptide
from ._shared.mn_utils_aa_atom_pos import MN_utils_aa_atom_pos
from .is_peptide import IsPeptide


class MN_animate_wiggle_mask_res(CustomGeometryGroup):
    _name = ".MN_animate_wiggle_mask_res"
    _tree_properties = {"node_tool_idname": "geometry._mn_animate_wiggle_mask_res"}

    def _build_group(self, tree):
        a = tree.inputs.integer("A", 0)
        result = tree.outputs.boolean("Result")

        group = MN_select_res_name_peptide(
            ala=True,
            arg=True,
            asn=True,
            asp=True,
            glu=True,
            gln=True,
            his=True,
            leu=True,
            lys=True,
            met=True,
            phe=True,
            trp=True,
        )
        group_1 = MN_select_res_name_peptide(
            ala=True,
            arg=True,
            asn=True,
            asp=True,
            cys=True,
            glu=True,
            gln=True,
            his=True,
            ile=True,
            leu=True,
            lys=True,
            met=True,
            phe=True,
            ser=True,
            thr=True,
            trp=True,
            tyr=True,
            val=True,
        )
        group_2 = MN_select_res_name_peptide(
            ala=True,
            arg=True,
            asn=True,
            asp=True,
            cys=True,
            glu=True,
            gln=True,
            his=True,
            ile=True,
            leu=True,
            lys=True,
            met=True,
            phe=True,
            ser=True,
            thr=True,
            trp=True,
            tyr=True,
            val=True,
        )
        (
            g.IndexSwitch.boolean(
                a,
                (
                    group_2.o.selection,
                    group_1.o.selection,
                    group.o.selection,
                    MN_select_res_name_peptide(
                        arg=True, gln=True, lys=True
                    ).o.selection,
                    MN_select_res_name_peptide(ile=True, lys=True).o.selection,
                ),
            )
            >> result
        )


class MN_animate_wiggle_mask_length(CustomGeometryGroup):
    _name = ".MN_animate_wiggle_mask_length"
    _tree_properties = {"node_tool_idname": "geometry._mn_animate_wiggle_mask_length"}

    def _build_group(self, tree):
        a = tree.inputs.integer("A", 0)
        result = tree.outputs.integer("Result")

        g.IndexSwitch.integer(a, (2, 5, 6, 12, 20)) >> result


class MN_animate_noise_repeat(CustomGeometryGroup):
    _name = "MN_animate_noise_repeat"
    _color_tag = "TEXTURE"
    _tree_properties = {"node_tool_idname": "geometry.mn_animate_noise_repeat"}

    def _build_group(self, tree):
        amplitude = tree.inputs.float(
            "Amplitude", 1.0, min_value=-10_000.0, max_value=10_000.0
        )
        detail = tree.inputs.float("Detail", 0.5, min_value=0.0, max_value=15.0)
        roughness = tree.inputs.float(
            "Roughness", 0.5, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        distortion = tree.inputs.float(
            "Distortion", 0.0, min_value=-1000.0, max_value=1000.0
        )
        vector = tree.inputs.vector(
            "Vector",
            (0.0, 0.0, 0.0),
            min_value=-10_000.0,
            max_value=10_000.0,
            hide_value=True,
            default_attribute="position",
        )
        speed = tree.inputs.float("Speed", 0.5, min_value=-10_000.0, max_value=10_000.0)
        animate_0_1 = tree.inputs.float(
            "Animate 0..1", 0.5, min_value=-10_000.0, max_value=10_000.0
        )
        noise_float = tree.outputs.float("Noise Float")
        noise_vector = tree.outputs.vector("Noise Vector")

        _value = g.Value(4.0)
        vector_math = (
            vector
            + g.CombineXYZ(x=animate_0_1 * speed)
            + g.RandomValue.vector((-10.0, -10.0, -10.0), (-1.0, 10.0, 10.0), vector)
        )
        vector_1 = vector_math * (math.tau / speed)
        noise_texture = g.NoiseTexture(
            vector=g.CombineXYZ(
                x=vector_1.x.sin(), y=vector_1.x.cos(), z=vector_1.y.sin()
            ),
            w=vector_1.y.cos(),
            scale=speed * g.Value(0.2),
            detail=detail,
            roughness=roughness.clamp(),
            distortion=distortion,
            noise_dimensions="4D",
            normalize=True,
        )
        map_range = g.MapRange(
            vector=noise_texture.o.color,
            to_min_float3=(-1.0, -1.0, -1.0),
            clamp=True,
            data_type="FLOAT_VECTOR",
        )
        noise_texture.o.factor.map_range(to_min=-1.0) * amplitude >> noise_float
        map_range.o.vector * amplitude >> noise_vector


class MN_utils_rotate_res(CustomGeometryGroup):
    _name = ".MN_utils_rotate_res"
    _tree_properties = {"node_tool_idname": "geometry._mn_utils_rotate_res"}

    def _build_group(self, tree):
        selection = tree.inputs.boolean(
            "Selection", False, description="Selection of atoms to apply this node to"
        )
        atom_name_rotation = tree.inputs.integer("atom_name rotation", 0)
        atom_name_axis = tree.inputs.integer("atom_name axis", 2)
        scale_b_factor = tree.inputs.float("Scale b_factor", 0.0, subtype="FACTOR")
        amplitude = tree.inputs.float("Amplitude", 1.0, min_value=0.0, max_value=10.0)
        amp_axis = tree.inputs.float(
            "Amp. Axis", 1.0, min_value=-10_000.0, max_value=10_000.0
        )
        amp_euler = tree.inputs.float(
            "Amp. Euler", 1.0, min_value=-10_000.0, max_value=10_000.0
        )
        speed = tree.inputs.float(
            "Speed", 10.0, min_value=-10_000.0, max_value=10_000.0
        )
        animate_0_1 = tree.inputs.float(
            "Animate 0..1", 0.5, min_value=-10_000.0, max_value=10_000.0
        )
        selection_1 = tree.outputs.boolean(
            "Selection", description="The calculated selection"
        )
        position = tree.outputs.vector("Position")

        group = MN_utils_aa_atom_pos(atom_name=atom_name_rotation)
        mix = g.Mix(
            factor_float=scale_b_factor,
            b_float=group.o.b_factor.map_range(
                1.0, 100.0, interpolation_type="SMOOTHERSTEP"
            ),
            a_float=1.0,
            clamp_factor=True,
        )
        math_1 = mix.o.result_float * amplitude
        named_attribute = g.NamedAttribute.integer("atom_name")
        boolean_math = (named_attribute.o.attribute > 4) & (
            (named_attribute.o.attribute > atom_name_rotation) & IsPeptide().o.selection
        )
        (
            (
                g.Compare.integer.not_equal(named_attribute.o.attribute, 38).o.result
                & boolean_math
                & selection
            )
            >> selection_1
        )
        random_value = g.RandomValue.vector(
            (-13.0, -13.0, -13.0),
            (13.9, 13.9, 13.9),
            group.o.group_index,
            atom_name_rotation,
        )
        group_1 = MN_animate_noise_repeat(
            Detail=2.0,
            Roughness=1.0,
            Distortion=1.98,
            Vector=random_value,
            Speed=speed,
            **{"Animate 0..1": animate_0_1},
        )
        group_2 = MN_animate_noise_repeat(
            Amplitude=amp_euler,
            Detail=1.0,
            Roughness=1.0,
            Distortion=1.98,
            Vector=random_value,
            Speed=speed,
            **{"Animate 0..1": animate_0_1},
        )
        vector_math = g.VectorMath.scale(
            group_2.o.noise_vector * math_1, g.Compare.integer.equal(group.o.integer, 3)
        )
        vector_math.node.mute = True
        vector_rotate = g.VectorRotate(
            vector=g.Position(),
            center=group.o.position,
            axis=group.o.position
            - MN_utils_aa_atom_pos(atom_name=atom_name_axis).o.position,
            angle=math_1 * amp_axis * group_1.o.noise_float,
        )
        vector_rotate_1 = g.VectorRotate.euler(
            vector_rotate, group.o.position, vector_math.o.vector
        )

        vector_rotate_1 >> position


class AnimateWiggle(AssetGeometryGroup):
    """
    Animate Wiggle

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this node to
    b_factor : InputFloat
        Amount that `b_factor` changeds the amplitude of wiggling
    amplitude : InputFloat
        Overall amplitude of the wiggling
    amp_axis : InputFloat
        Aplitude for the rotation around the bond axes
    amp_euler : InputFloat
        Amplitude for applying euler rotations separate to the axis
    speed : InputFloat
        Speed at which the wiggle is applied, 3 will repeat 3 times
    animate : InputFloat
        Controls the animation of the wiggle, repeating every `1.00`

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.b_factor : FloatSocket
        Amount that `b_factor` changeds the amplitude of wiggling
    i.amplitude : FloatSocket
        Overall amplitude of the wiggling
    i.amp_axis : FloatSocket
        Aplitude for the rotation around the bond axes
    i.amp_euler : FloatSocket
        Amplitude for applying euler rotations separate to the axis
    i.speed : FloatSocket
        Speed at which the wiggle is applied, 3 will repeat 3 times
    i.animate : FloatSocket
        Controls the animation of the wiggle, repeating every `1.00`

    Outputs
    -------
    o.atoms : GeometrySocket
        The animated atomic geometry
    """

    _name = "Animate Wiggle"
    _asset_name = "Animate Wiggle"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.animate_wiggle"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        b_factor: FloatSocket
        """Amount that `b_factor` changeds the amplitude of wiggling"""
        amplitude: FloatSocket
        """Overall amplitude of the wiggling"""
        amp_axis: FloatSocket
        """Aplitude for the rotation around the bond axes"""
        amp_euler: FloatSocket
        """Amplitude for applying euler rotations separate to the axis"""
        speed: FloatSocket
        """Speed at which the wiggle is applied, 3 will repeat 3 times"""
        animate: FloatSocket
        """Controls the animation of the wiggle, repeating every `1.00`"""

    class _Outputs(SocketAccessor):
        atoms: GeometrySocket
        """The animated atomic geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        b_factor: InputFloat = 1.0,
        amplitude: InputFloat = 1.0,
        amp_axis: InputFloat = 1.0,
        amp_euler: InputFloat = 0.4,
        speed: InputFloat = 3.0,
        animate: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "b_factor": b_factor,
                "Amplitude": amplitude,
                "Amp. Axis": amp_axis,
                "Amp. Euler": amp_euler,
                "Speed": speed,
                "Animate": animate,
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
        b_factor = tree.inputs.float(
            "b_factor",
            1.0,
            description="Amount that `b_factor` changeds the amplitude of wiggling",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        amplitude = tree.inputs.float(
            "Amplitude",
            1.0,
            description="Overall amplitude of the wiggling",
            min_value=0.0,
            max_value=10.0,
        )
        amp_axis = tree.inputs.float(
            "Amp. Axis",
            1.0,
            description="Aplitude for the rotation around the bond axes",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        amp_euler = tree.inputs.float(
            "Amp. Euler",
            0.4,
            description="Amplitude for applying euler rotations separate to the axis",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        speed = tree.inputs.float(
            "Speed",
            3.0,
            description="Speed at which the wiggle is applied, 3 will repeat 3 times",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        animate = tree.inputs.float(
            "Animate",
            0.0,
            description="Controls the animation of the wiggle, repeating every `1.00`",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        atoms_1 = tree.outputs.geometry(
            "Atoms", description="The animated atomic geometry"
        )

        repeat_zone = g.RepeatZone(5)
        geometry = repeat_zone.items.geometry("Geometry", atoms)
        integer = repeat_zone.items.integer("Integer")
        group = MN_utils_rotate_res(
            Selection=MN_animate_wiggle_mask_res(A=integer.current).o.result
            & selection,
            **{
                "atom_name rotation": MN_animate_wiggle_mask_length(A=integer.current),
                "atom_name axis": MN_animate_wiggle_mask_length(
                    A=integer.current - 1.0
                ),
                "Scale b_factor": b_factor,
            },
            Amplitude=amplitude,
            **{"Amp. Axis": amp_axis, "Amp. Euler": amp_euler},
            Speed=speed,
            **{"Animate 0..1": animate},
        )
        set_position = g.SetPosition(
            geometry=geometry.current,
            selection=group.o.selection,
            position=group.o.position,
        )
        set_position >> geometry.next
        integer.current + 1.0 >> integer.next

        geometry.result >> atoms_1


ASSET = AnimateWiggle

ASSET_METADATA = {
    "catalog_id": "85730213-4c2e-469f-b333-52ac53adf274",
}
