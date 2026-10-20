from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import CustomGeometryGroup
from molecularnodes.nodes.geometry.angstrom_to_world import AngstromToWorld


class StyleIllustration(CustomGeometryGroup):
    """Flat vector illustration of atoms and bonds as Grease Pencil strokes."""

    _name = "Style Illustration"
    _color_tag = "GEOMETRY"

    def __init__(self, atoms=None, **kwargs):
        names = {k: k.replace("_", " ").title() for k in kwargs}
        super().__init__(**{"Atoms": atoms, **{names[k]: v for k, v in kwargs.items()}})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry("Atoms")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        circles = tree.inputs.boolean(
            "Circles",
            True,
            hide_value=True,
            description="Atoms drawn as circles. Bonds still connect atoms without one",
        )
        camera = tree.inputs.object(
            "Camera",
            description="Circles face this camera and strokes are depth sorted from it",
        )
        atom_scale = tree.inputs.float(
            "Atom Scale",
            0.3,
            min_value=0.0,
            description="Scale of `vdw_radii` for the circle radius",
        )
        outline = tree.inputs.float(
            "Outline Thickness",
            0.06,
            min_value=0.0,
            description="Circle stroke width in Angstroms",
        )
        bond_thickness = tree.inputs.float(
            "Bond Thickness",
            0.12,
            min_value=0.0,
            description="Bond stroke width in Angstroms",
        )
        bond_gap = tree.inputs.float(
            "Bond Gap",
            0.0,
            description="Extra space between a bond and its atom's circle, in Angstroms",
        )
        bond_orders = tree.inputs.boolean(
            "Bond Orders",
            True,
            description="Draw double and triple bonds as parallel lines",
        )
        bond_spacing = tree.inputs.float(
            "Bond Spacing",
            0.3,
            min_value=0.0,
            description="Distance between parallel lines of multiple bonds, in Angstroms",
        )
        fill = tree.inputs.boolean(
            "Fill",
            True,
            description="Fill circles. Each circle is its own fill group (`fill_id`) so overlaps do not cut holes",
        )
        with tree.inputs.panel("Color"):
            fill_color = tree.inputs.color("Fill Color", (1.0, 1.0, 1.0, 1.0))
            fill_element = tree.inputs.float(
                "Fill Element",
                1.0,
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
                description="Mix the fill towards the atom's `Color`",
            )
            stroke_color = tree.inputs.color("Stroke Color", (0.0, 0.0, 0.0, 1.0))
            stroke_element = tree.inputs.float(
                "Stroke Element", 0.0, min_value=0.0, max_value=1.0, subtype="FACTOR"
            )
            bond_color = tree.inputs.color("Bond Color", (0.0, 0.0, 0.0, 1.0))
            bond_element = tree.inputs.float(
                "Bond Element",
                0.0,
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
                description="Mix each half of a bond towards its atom's `Color`",
            )
        with tree.inputs.panel("Material", default_closed=True):
            atom_material = tree.inputs.material(
                "Atom Material", description="Grease Pencil material for circles"
            )
            bond_material = tree.inputs.material(
                "Bond Material", description="Grease Pencil material for bonds"
            )
        grease_pencil = tree.outputs.geometry("Grease Pencil")

        cam = g.ObjectInfo(camera, transform_space="RELATIVE")
        view = g.RotateVector((0.0, 0.0, -1.0), cam.o.rotation)
        atoms_sel = g.SeparateGeometry.point(atoms, selection=selection)
        is_circle = g.NamedAttribute.boolean("_circle").o.attribute
        color = g.NamedAttribute.color("Color").o.attribute
        radius = g.NamedAttribute.float("_radius").o.attribute
        depth = g.NamedAttribute.float("_depth").o.attribute
        position = g.Position()

        with g.Frame("Atoms"):
            g.String(
                string="Every atom gets a world-space radius (0 without a circle) and its "
                "depth along the camera view, used to trim and to sort bonds and circles."
            )
            atoms_r = (
                atoms_sel.o.selection
                >> g.StoreNamedAttribute.point.boolean(name="_circle", value=circles)
                >> g.StoreNamedAttribute.point.float(
                    name="_radius",
                    value=g.Switch.float(
                        is_circle,
                        0.0,
                        g.NamedAttribute.float("vdw_radii").o.attribute * atom_scale,
                    ),
                )
                >> g.StoreNamedAttribute.point.float(
                    name="_depth", value=(position - cam.o.location).dot(view)
                )
            )

        with g.Frame("Circles"):
            g.String(
                string="Unit circles instanced on each atom, rotated to face the camera "
                "and scaled by the atom's radius."
            )
            atom_curves = (
                g.InstanceOnPoints(
                    points=g.MeshToPoints(atoms_r, selection=is_circle),
                    instance=g.CurveCircle(resolution=32, radius=1.0),
                    rotation=cam.o.rotation,
                    scale=radius,
                )
                >> g.RealizeInstances()
                >> g.StoreNamedAttribute.spline.boolean(name="_circle", value=True)
                >> g.SetCurveRadius(radius=AngstromToWorld(angstrom=outline) * 0.5)
            )

        with g.Frame("Bonds"):
            g.String(
                string="Each bond becomes 2 halves times its order, each a two-point "
                "curve from its own atom (trimmed to the projected circle) to the bond "
                "midpoint. A half is sorted just behind its own atom so circles cover "
                "their bonds while nearer atoms' bonds still cross farther atoms."
            )
            edge_points = g.MeshToPoints(
                atoms_r, mode="EDGES"
            ) >> g.StoreNamedAttribute.point.integer(name="_edge", value=g.Index())
            bond_type = g.NamedAttribute.integer("bond_type").o.attribute
            is_double = g.Compare.integer.equal(bond_type, 2) | g.Compare.integer.equal(
                bond_type, 6
            )
            is_triple = g.Compare.integer.equal(bond_type, 3) | g.Compare.integer.equal(
                bond_type, 7
            )
            order = g.Switch.integer(
                bond_orders,
                1,
                g.Switch.integer(is_triple, g.Switch.integer(is_double, 1, 2), 3),
            )
            dupes = g.DuplicateElements(edge_points, amount=order * 4)
            d = dupes.o.duplicate_index
            edge = g.NamedAttribute.integer("_edge").o.attribute
            endpoint = d % 2
            half = (d // 2) % 2
            line = d // 4
            n_lines = g.SampleIndex.point.integer(edge_points, order, edge)
            ev = g.EdgeVertices()
            v_own = g.Switch.integer(
                half == 1,
                g.SampleIndex.edge.integer(atoms_r, ev.o.vertex_index_1, edge),
                g.SampleIndex.edge.integer(atoms_r, ev.o.vertex_index_2, edge),
            )
            v_other = g.Switch.integer(
                half == 1,
                g.SampleIndex.edge.integer(atoms_r, ev.o.vertex_index_2, edge),
                g.SampleIndex.edge.integer(atoms_r, ev.o.vertex_index_1, edge),
            )
            p_own = g.SampleIndex.point.vector(atoms_r, position, v_own)
            p_other = g.SampleIndex.point.vector(atoms_r, position, v_other)
            r_own = g.SampleIndex.point.float(atoms_r, radius, v_own)
            bond = p_other - p_own
            direction = bond.normalize()
            side = direction.cross(view)
            sin_theta = g.Math.maximum(side.length(), 0.2)
            half_length = bond.length() * 0.5
            trim = g.Math.minimum(
                (r_own + AngstromToWorld(angstrom=bond_gap)) / sin_theta, half_length
            )
            along = g.Switch.float(endpoint == 1, trim, half_length)
            offset = (g.Math.subtract(line, (n_lines - 1) * 0.5)) * AngstromToWorld(
                angstrom=bond_spacing
            )
            bond_points = (
                dupes.o.geometry
                >> g.SetPosition(
                    position=p_own + direction * along + side.normalize() * offset
                )
                >> g.StoreNamedAttribute.point.float(
                    name="_depth",
                    value=g.SampleIndex.point.float(atoms_r, depth, v_own) + 0.0001,
                )
                >> g.StoreNamedAttribute.point.color(
                    name="Color", value=g.SampleIndex.point.color(atoms_r, color, v_own)
                )
            )
            bond_curves = (
                g.PointsToCurves(
                    bond_points, curve_group_id=edge * 8 + d // 2, weight=endpoint
                )
                >> g.DeleteGeometry(
                    selection=g.SplineLength().o.length < 0.001, domain="CURVE"
                )
                >> g.SetCurveRadius(
                    radius=AngstromToWorld(angstrom=bond_thickness) * 0.5
                )
                >> g.StoreNamedAttribute.spline.boolean(name="_circle", value=False)
            )

        is_atom = g.NamedAttribute.boolean("_circle").o.attribute
        stroke = g.Switch.color(
            is_atom,
            g.Mix.color(bond_element, bond_color, color),
            g.Mix.color(stroke_element, stroke_color, color),
        )
        (
            g.JoinGeometry((bond_curves, atom_curves))
            >> g.SortElements(sort_weight=-depth, domain="CURVE")
            >> g.StoreNamedAttribute.spline.integer(
                name="fill_id", value=g.Switch.integer(is_atom & fill, 0, g.Index() + 1)
            )
            >> g.CurvesToGreasePencil(instances_as_layers=False)
            >> g.SetMaterial(selection=is_atom, material=atom_material)
            >> g.SetMaterial(selection=~is_atom, material=bond_material)
            >> g.SetGreasePencilColor(color=stroke, mode="STROKE")
            >> g.SetGreasePencilColor(
                selection=is_atom,
                color=g.Mix.color(fill_element, fill_color, color),
                mode="FILL",
            )
            >> g.SetGreasePencilDepth(depth_order="2D")
            >> grease_pencil
        )
