# Node-group asset "Animate Reveal" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry, InputMenu
from .vdw_radii import VDWRadii


class AnimateReveal(AssetGeometryGroup):
    """
    Animate Reveal

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this node to
    factor : InputFloat
        How revealed each atom is, from 0 (hidden) to 1 (shown). Link `Animate Value` or `Animate Stagger` to animate it
    mode : InputMenu | Literal["Alpha", "Scale", "Cull"]
        Alpha multiplies the alpha of the `Color` attribute (honoured by the Default and Squishy materials), Scale multiplies `vdw_radii` so sphere styles shrink, Cull deletes atoms whose factor is 0
    invert : InputBoolean
        Hide instead of reveal, using 1 - `Factor`

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.factor : FloatSocket
        How revealed each atom is, from 0 (hidden) to 1 (shown). Link `Animate Value` or `Animate Stagger` to animate it
    i.mode : MenuSocket
        Alpha multiplies the alpha of the `Color` attribute (honoured by the Default and Squishy materials), Scale multiplies `vdw_radii` so sphere styles shrink, Cull deletes atoms whose factor is 0
    i.invert : BooleanSocket
        Hide instead of reveal, using 1 - `Factor`

    Outputs
    -------
    o.atoms : GeometrySocket
        Atomic geometry revealed by the factor
    """

    _name = "Animate Reveal"
    _asset_name = "Animate Reveal"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        factor: FloatSocket
        """How revealed each atom is, from 0 (hidden) to 1 (shown). Link `Animate Value` or `Animate Stagger` to animate it"""
        mode: MenuSocket
        """Alpha multiplies the alpha of the `Color` attribute (honoured by the Default and Squishy materials), Scale multiplies `vdw_radii` so sphere styles shrink, Cull deletes atoms whose factor is 0"""
        invert: BooleanSocket
        """Hide instead of reveal, using 1 - `Factor`"""

    class _Outputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry revealed by the factor"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        factor: InputFloat = 1.0,
        mode: InputMenu | Literal["Alpha", "Scale", "Cull"] = "Alpha",
        invert: InputBoolean = False,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Factor": factor,
                "Mode": mode,
                "Invert": invert,
            }
        )

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
        factor = tree.inputs.float(
            "Factor",
            1.0,
            description="How revealed each atom is, from 0 (hidden) to 1 (shown). Link `Animate Value` or `Animate Stagger` to animate it",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        mode = tree.inputs.menu(
            "Mode",
            description="Alpha multiplies the alpha of the `Color` attribute (honoured by the Default and Squishy materials), Scale multiplies `vdw_radii` so sphere styles shrink, Cull deletes atoms whose factor is 0",
        )
        invert = tree.inputs.boolean(
            "Invert", False, description="Hide instead of reveal, using 1 - `Factor`"
        )
        atoms_1 = tree.outputs.geometry(
            "Atoms", description="Atomic geometry revealed by the factor"
        )

        clamp = factor.clamp()
        switch = invert.switch.float(clamp, 1.0 - clamp)
        with g.Frame("Alpha"):
            attribute = g.NamedAttribute.color("Color").o.attribute
            _string = g.String(
                string="Keeps the colour and multiplies the existing alpha of the Color attribute, so an earlier alpha (or a second Animate Reveal) still applies."
            )
            combine_color = g.CombineColor(
                red=attribute.r,
                green=attribute.g,
                blue=attribute.b,
                alpha=attribute.a * switch,
            )
            store_named_attribute = g.StoreNamedAttribute.point.color(
                atoms, selection, "Color", combine_color
            )
        with g.Frame("Scale"):
            _string_1 = g.String(
                string="Style Spheres and Style Ball and Stick size their spheres from vdw_radii, so scaling it grows the atoms in. Bonds, ribbons and cartoons are unaffected."
            )
            store_named_attribute_1 = g.StoreNamedAttribute.point.float(
                atoms, selection, "vdw_radii", VDWRadii().o.vdw_radii * switch
            )
        with g.Frame("Cull"):
            delete_geometry = g.DeleteGeometry.point(atoms, selection & (switch <= 0.0))
        (
            g.MenuSwitch.geometry(
                mode,
                {
                    "Alpha": store_named_attribute,
                    "Scale": store_named_attribute_1,
                    "Cull": delete_geometry,
                },
            )
            >> atoms_1
        )

        mode.default_value = "Alpha"


ASSET = AnimateReveal

ASSET_METADATA = {
    "description": "Reveal or hide atoms by a 0..1 factor through the Color alpha, the sphere radius or by deleting them",
    "catalog_id": "85730213-4c2e-469f-b333-52ac53adf274",
}
