from pathlib import Path
import bpy
import MDAnalysis as mda
from bpy.props import (
    BoolProperty,
    CollectionProperty,
    EnumProperty,
    FloatProperty,
    IntProperty,
    StringProperty,
)
from bpy.types import Context, Operator
from .. import entities
from ..annotations.props import create_annotation_type_inputs
from ..blender.utils import path_resolve
from ..download import CACHE_DIR, FileDownloadPDBError
from ..entities import (
    Molecule,
    StreamingTrajectory,
    Trajectory,
    density,
    ensemble,
    trajectory,
)
from ..nodes import geometry as g
from ..nodes import nodes
from ..nodes.material import add_all_materials
from ..nodes.node_management import (
    remove_style_node,
)
from ..scene.compositor import setup_compositor
from ..session import get_session
from .pref import addon_preferences
from .style import STYLE_ITEMS


def _add_node(node_name, context, show_options=False, material="default"):
    """
    Add a node group to the node tree and set the values.

    intended to be called upon button press in the node tree, and not for use in general scripting
    """

    # actually invoke the operator to add a node to the current node tree
    # use_transform=True ensures it appears where the user's mouse is and is currently
    # being moved so the user can place it where they wish
    bpy.ops.node.add_node(
        "INVOKE_DEFAULT", type="GeometryNodeGroup", use_transform=True
    )
    node = context.active_node
    node.node_tree = bpy.data.node_groups[node_name]
    # set the label to the node tree name by default
    node.label = node.node_tree.name
    node.width = nodes.NODE_WIDTH
    node.show_options = show_options
    node.name = node_name

    # if added node has a 'Material' input, set it to the default MN material
    nodes.assign_material(node, new_material=material)


_STYLE_NODE = {
    "spheres": g.StyleSpheres,
    "ribbon": g.StyleRibbon,
    "cartoon": g.StyleCartoon,
    "sticks": g.StyleSticks,
    "ball_and_stick": g.StyleBallAndStick,
    "surface": g.StyleSurface,
}


class Import_Molecule(bpy.types.Operator):
    style: EnumProperty(  # type: ignore
        name="Style",
        default="spheres",
        description="Starting style for the structure on import",
        items=STYLE_ITEMS,
    )

    node_setup: BoolProperty(  # type: ignore
        name="Node Setup",
        default=True,
        description="Whether to setup the starting default node tree on import",
    )
    assembly: BoolProperty(  # type: ignore
        default=False,
        name="Build Biological Assembly",
        description="Build the biological assembly for the structure on import",
    )

    @property
    def style_node(
        self,
    ) -> type[
        g.StyleSpheres
        | g.StyleSurface
        | g.StyleRibbon
        | g.StyleSticks
        | g.StyleBallAndStick
        | g.StyleCartoon
    ]:
        "Helper to get the selected node class for adding to the tree"
        return _STYLE_NODE[self.style]

    def draw(self, context):
        layout = self.layout
        assert layout
        row = layout.row()
        row.prop(self, "node_setup", text="")
        col = row.column()
        col.prop(self, "style")
        col.enabled = self.node_setup
        # row = layout.row()
        layout.prop(self, "assembly")

        return layout


class MN_OT_Import_Molecule(Import_Molecule):
    bl_idname = "mn.import_molecule"
    bl_label = "Import a Molecule"
    bl_options = {"REGISTER", "UNDO"}

    directory: StringProperty(  # type: ignore
        subtype="FILE_PATH", options={"SKIP_SAVE", "HIDDEN"}
    )
    files: CollectionProperty(  # type: ignore
        type=bpy.types.OperatorFileListElement, options={"SKIP_SAVE", "HIDDEN"}
    )

    def draw(self, context):
        layout = self.layout
        assert layout
        layout.label(text=f"Importing {len(self.files)} molecules")
        layout = super().draw(context)

    def execute(self, context):
        if not self.directory:
            return {"CANCELLED"}

        for file in self.files:
            try:
                mol = Molecule.load(
                    file_path=Path(self.directory, file.name), name=file.name
                )

                with mol.tree.reset() as (atoms, join):
                    (
                        atoms
                        >> self.style_node(material=add_all_materials()["MN Default"])
                        >> (
                            g.AssemblyInstance(data_object=mol.create_data_object())
                            if self.assembly
                            else None
                        )
                        >> join
                    )
            except Exception as e:
                print(f"Failed importing {file}: {e}")

        return {"FINISHED"}

    def invoke(self, context, event):
        if context.area and context.area.type == "VIEW_3D":
            context.window_manager.invoke_props_dialog(self)
        else:
            context.window_manager.fileselect_add(self)
        return {"RUNNING_MODAL"}


class MN_FH_Import_Molecule(bpy.types.FileHandler):
    bl_idname = "MN_FH_import_molecule"
    bl_label = "File handler for import molecular data files."
    bl_import_operator = "mn.import_molecule"
    bl_file_extensions = ".pdb;.cif;.mmcif;.bcif;.pdbx"

    @classmethod
    def poll_drop(cls, context):
        return context.area and context.area.type == "VIEW_3D"


DOWNLOAD_FORMATS = (
    ("bcif", ".bcif", "Binary compressed .cif file, fastest for downloading"),
    ("cif", ".cif", "The new standard of .cif / .mmcif"),
    ("pdb", ".pdb", "The classic (and depcrecated) PDB format"),
)


# operator that is called by the 'button' press which calls the fetch function


class MN_OT_Import_Fetch(Import_Molecule, bpy.types.Operator):
    bl_idname = "mn.import_fetch"
    bl_label = "Import Molecule"
    bl_description = "Open a local structure file or download one from a database"
    bl_options = {"REGISTER", "UNDO"}

    database: EnumProperty(  # type: ignore
        name="Database",
        default="wwpdb",
        items=(
            (
                "local",
                "Local File",
                "Open a structure file already on disk",
            ),
            (
                "wwpdb",
                "wwPDB",
                "The world-wide Protein Data Bank (wwPDB)",
            ),
            (
                "alphafold",
                "AlphaFold",
                "The AlphaFold computational structure database",
            ),
        ),
    )
    code: StringProperty(  # type: ignore
        name="PDB",
        description="The PDB code to download (4-character e.g. '1abc' or 12-character e.g. 'pdb_00001abc')",
        options={"TEXTEDIT_UPDATE"},
    )
    filepath: StringProperty(  # type: ignore
        name="File",
        description="Path of the local structure file to open",
        subtype="FILE_PATH",
    )
    file_format: EnumProperty(  # type: ignore
        name="Format",
        description="Format to download as from the PDB",
        default="bcif",
        items=DOWNLOAD_FORMATS,
    )

    cache_dir: StringProperty(  # type: ignore
        name="Cache Directory",
        description="Where to store the structures downloaded from the Protein Data Bank",
        default=str(CACHE_DIR),
        subtype="DIR_PATH",
    )

    def draw(self, context):
        layout = self.layout
        assert layout
        layout.prop_tabs_enum(self, "database")
        if self.database == "local":
            layout.prop(self, "filepath")
        else:
            row = layout.row().split(factor=0.7)
            row.prop(self, "code")
            # file format only applies to wwPDB downloads; others pick their own
            if self.database == "wwpdb":
                row.prop(self, "file_format", text="")
        row = layout.row()
        row.prop(self, "node_setup", text="")
        col = row.column()
        col.prop(self, "style")
        col.enabled = self.node_setup
        layout.prop(self, "assembly")

    def invoke(self, context, event):
        prefs = addon_preferences()
        self.cache_dir = str(prefs.cache_dir) if prefs is not None else bpy.app.tempdir
        return context.window_manager.invoke_props_dialog(self)

    def execute(self, context):
        try:
            if self.database == "local":
                mol = Molecule.load(file_path=path_resolve(self.filepath))
                message = f"Imported '{self.filepath}' as {mol.name}"
            else:
                mol = entities.Molecule.fetch(
                    code=self.code,
                    cache=self.cache_dir,
                    format=self.file_format,
                    database=self.database,
                )
                message = f"Downloaded {self.code} as {mol.name}"
        except FileDownloadPDBError as e:
            self.report({"ERROR"}, str(e))
            if self.file_format == "pdb":
                self.report(
                    {"ERROR"},
                    "There may not be a `.pdb` formatted file available - try a different download format.",
                )
            return {"CANCELLED"}

        if self.assembly:
            nodes.assembly_data_object_from_obj(mol.object)

        with mol.tree.reset() as (atoms, join):
            (
                atoms
                >> g.SetColor(color=g.ColorElement(c=g.RandomColor(g.ChainID(), 3)))
                >> self.style_node(material=add_all_materials()["MN Default"])
                >> (
                    g.AssemblyInstance(data_object=mol.create_data_object())
                    if self.assembly
                    else None
                )
                >> join
            )

        try:
            bpy.context.view_layer.objects.active = mol.object  # type: ignore
        except RuntimeError:
            message += " - Molecular Nodes collection is disabled"

        self.report({"INFO"}, message=message)

        return {"FINISHED"}


ENSEMBLE_TYPES = (
    ("starfile", "Starfile", "Import a .star mapback file"),
    ("cellpack", "CellPack", "Import a CellPack .cif / .bcif model"),
)


class MN_OT_Import_Ensemble(bpy.types.Operator):
    bl_idname = "mn.import_ensemble"
    bl_label = "Import Ensemble"
    bl_description = "Import an ensemble as a Starfile or CellPack model"
    bl_options = {"REGISTER", "UNDO"}

    ensemble_type: EnumProperty(  # type: ignore
        name="Type",
        description="The kind of ensemble to import",
        default="starfile",
        items=ENSEMBLE_TYPES,
    )
    filepath: StringProperty(  # type: ignore
        name="File",
        description="File path for the ensemble to import",
        subtype="FILE_PATH",
        maxlen=0,
    )
    node_setup: BoolProperty(  # type: ignore
        name="Setup Nodes",
        default=True,
        description="Create and set up a Geometry Nodes tree on import",
    )

    def draw(self, context):
        layout = self.layout
        assert layout
        layout.prop_tabs_enum(self, "ensemble_type")
        layout.prop(self, "filepath")
        layout.prop(self, "node_setup")

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)

    def execute(self, context):
        file_path = path_resolve(self.filepath)
        if self.ensemble_type == "cellpack":
            ensemble.load_cellpack(
                file_path=file_path,
                name=Path(self.filepath).name,
                node_setup=self.node_setup,
            )
        else:
            ensemble.load_starfile(
                file_path=file_path,
                node_setup=self.node_setup,
            )
        return {"FINISHED"}


DENSITY_STYLE_ITEMS = (
    (
        "density_surface",
        "Surface",
        "A mesh surface based on the specified threshold",
        0,
    ),
    (
        "density_iso_surface",
        "ISO Surface",
        "A mesh surface based on the specified iso value",
        1,
    ),
    (
        "density_wire",
        "Wire",
        "A wire mesh surface based on the specified threshold",
        2,
    ),
)


class MN_OT_Import_Map(bpy.types.Operator):
    bl_idname = "mn.import_density"
    bl_label = "Import Density"
    bl_description = "Import an EM density map into Blender"
    bl_options = {"REGISTER", "UNDO"}

    filepath: StringProperty(  # type: ignore
        name="File",
        description="File path for the map file.",
        subtype="FILE_PATH",
        maxlen=0,
    )
    invert: BoolProperty(  # type: ignore
        name="Invert Data",
        description="Invert the values in the map. Low becomes high, high becomes low.",
        default=False,
    )
    center: BoolProperty(  # type: ignore
        name="Center Density",
        description="Translate the density so that the center of the box is at the origin.",
        default=False,
    )
    overwrite: BoolProperty(  # type: ignore
        name="Overwrite Intermediate File",
        description="Overwrite generated intermediate .vdb file.",
        default=False,
    )
    setup_nodes: BoolProperty(  # type: ignore
        name="Setup Nodes",
        default=True,
        description="Create and set up a Geometry Nodes tree on import",
    )
    style: EnumProperty(  # type: ignore
        name="Style",
        items=DENSITY_STYLE_ITEMS,
    )

    def draw(self, context):
        layout = self.layout
        assert layout
        layout.prop(self, "filepath")
        layout.prop(self, "invert")
        layout.prop(self, "center")
        layout.prop(self, "overwrite")
        row = layout.row()
        row.prop(self, "setup_nodes", text="")
        col = row.column()
        col.prop(self, "style")
        col.enabled = self.setup_nodes

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)

    def execute(self, context):
        density.load(
            file_path=path_resolve(self.filepath),
            invert=self.invert,
            style=self.style if self.setup_nodes else None,
            center=self.center,
            overwrite=self.overwrite,
        )
        return {"FINISHED"}


class MN_OT_Reload_Trajectory(bpy.types.Operator):
    bl_idname = "mn.reload_trajectory"
    bl_label = "Reload Trajectory"
    bl_description = (
        "Reload the `mda.UNiverse` of the current Object to renable updating"
    )
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        loaded_trajectory = context.scene.MNSession.match(obj)
        return obj.mn.entity_type.startswith("md") and not loaded_trajectory

    def execute(self, context):
        obj = context.active_object
        path_topo = path_resolve(obj.mn.filepath_topology)
        path_traj = path_resolve(obj.mn.filepath_trajectory)

        if "oxdna" in obj.mn.entity_type:
            uni = mda.Universe(
                path_topo,
                path_traj,
                topology_format=trajectory.oxdna.OXDNAParser,
                format=trajectory.oxdna.OXDNAReader,
            )
            traj = trajectory.oxdna.OXDNA(uni, create_object=False)
        elif "streaming" in obj.mn.entity_type:
            traj = StreamingTrajectory.load(path_topo, path_traj, create_object=False)
        else:
            traj = Trajectory.load(
                path_topo, path_traj, style=None, create_object=False
            )

        traj.object = obj
        traj.set_frame(context.scene.frame_current)
        return {"FINISHED"}


class MN_OT_Import_Trajectory(bpy.types.Operator):
    bl_idname = "mn.import_trajectory"
    bl_label = "Import Trajectory"
    bl_description = "Load a molecular dynamics or oxDNA trajectory"
    bl_options = {"REGISTER", "UNDO"}

    format: EnumProperty(  # type: ignore
        name="Format",
        description="The kind of trajectory to import",
        default="md",
        items=(
            ("md", "MD", "A molecular dynamics trajectory (via MDAnalysis)"),
            ("oxdna", "oxDNA", "An oxDNA trajectory"),
        ),
    )
    topology: StringProperty(  # type: ignore
        name="Topology",
        description="File path for the toplogy file for the trajectory",
        subtype="FILE_PATH",
        maxlen=0,
    )
    trajectory: StringProperty(  # type: ignore
        name="Trajectory",
        description="File path for the trajectory file for the trajectory",
        subtype="FILE_PATH",
        maxlen=0,
    )
    name: StringProperty(  # type: ignore
        name="Name",
        description="Name of the molecule on import",
        default="NewTrajectory",
        maxlen=0,
    )
    style: EnumProperty(  # type: ignore
        name="Style",
        description="Default style for importing",
        items=STYLE_ITEMS,
        default="spheres",
    )
    setup_nodes: BoolProperty(  # type: ignore
        name="Setup Nodes",
        description="Add nodes to the scene to load the trajectory",
        default=True,
    )

    def draw(self, context):
        layout = self.layout
        assert layout
        layout.prop_tabs_enum(self, "format")
        layout.prop(self, "name")
        layout.prop(self, "topology")
        layout.prop(self, "trajectory")
        # oxDNA imports don't set up a style/node tree in the same way
        if self.format == "md":
            row = layout.row()
            row.prop(self, "setup_nodes", text="")
            col = row.column()
            col.prop(self, "style")
            col.enabled = self.setup_nodes

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)

    def execute(self, context):
        if self.format == "oxdna":
            trajectory.load_oxdna(
                top=path_resolve(self.topology),
                traj=path_resolve(self.trajectory),
                name=self.name,
            )
            return {"FINISHED"}

        topology = path_resolve(self.topology)
        coordinates = path_resolve(self.trajectory)

        if self.trajectory.startswith("imd://"):
            traj = StreamingTrajectory.load(
                topology=topology,
                coordinates=coordinates,
                name=self.name,
                style=self.style,
                selection="all",
            )
        else:
            traj = Trajectory.load(
                topology=topology,
                coordinates=coordinates,
                name=self.name,
                style=self.style if self.setup_nodes else None,
                selection="all",
            )

        context.view_layer.objects.active = traj.object
        context.scene.frame_start = 0

        if isinstance(traj, StreamingTrajectory):
            self.report(
                {"INFO"},
                message=f"Streaming trajectory '{traj.name}' from '{self.trajectory}'",
            )
        else:
            n_frames = int(traj.object.mn.n_frames - 1)
            context.scene.frame_end = n_frames
            self.report(
                {"INFO"},
                message=f"Imported '{self.topology}' as {traj.name} "
                f"with {n_frames} frames from '{self.trajectory}'.",
            )

        return {"FINISHED"}


class MN_OT_Add_Style(Operator):
    """
    Operator to add a new style to an entity
    """

    bl_idname = "mn.add_style"
    bl_label = "Add Style"
    bl_description = "Add new style to Fpointntity"

    uuid: StringProperty()  # type: ignore

    style: EnumProperty(  # type: ignore
        name="Style",
        default="spheres",
        description="Style type",
        items=STYLE_ITEMS,
    )

    use_uniform_color: BoolProperty(  # type: ignore
        name="Use uniform color",
        description="Use uniform color for style",
        default=False,
    )

    uniform_color: bpy.props.FloatVectorProperty(  # type: ignore
        name="Uniform color",
        description="Uniform color for style",
        subtype="COLOR",
        size=4,
        default=(0.162, 0.624, 0.196, 1),
        min=0.0,
        max=1.0,
    )

    color_scheme: StringProperty(  # type: ignore
        name="Coloring scheme",
        description="Coloring scheme for style",
        default="common",
    )

    selection: StringProperty(  # type: ignore
        name="Selection",
        description="Selection for which the style applies",
        default="all",
    )

    name: StringProperty(  # type: ignore
        name="Name",
        description="Label for the style",
        default="",
    )

    def draw(self, context):
        layout = self.layout
        assert layout
        layout.prop(self, "style")
        layout.prop(self, "use_uniform_color")
        if self.use_uniform_color:
            layout.prop(self, "uniform_color")
        else:
            layout.prop(self, "color_scheme")
        layout.prop(self, "selection")
        layout.prop(self, "name")
        # Note: Materials cannot be passed into operators

    def execute(self, context: Context):
        entity = get_session().get(self.uuid)
        if self.use_uniform_color:
            color = self.uniform_color
        else:
            color = self.color_scheme
        entity.add_style(
            style=self.style,
            color=color,
            selection=self.selection.strip() or None,
            name=self.name.strip() or None,
        )
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)


class MN_OT_Remove_Style(Operator):
    """
    Operator to remove a style from an entity
    """

    bl_idname = "mn.remove_style"
    bl_label = "Remove Style"
    bl_description = "Remove style from entity"

    uuid: StringProperty()  # type: ignore
    style_node_index: IntProperty()  # type: ignore

    def execute(self, context: Context):
        entity = get_session().get(self.uuid)
        node_group = entity.node_group
        style_node = node_group.nodes[self.style_node_index]
        remove_style_node(style_node)
        # set the active index in UI to the last style

        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_confirm(self, event, title="Remove Style?")


def _register_temp_annotation_add_op(entity):
    """Register a temporary annotation add operator with custom input properties"""

    def get_annotation_types(self, context):
        annotation_types = []
        for type in entity.annotations._classes.keys():
            annotation_types.append((type, type, f"Annotation of type {type}"))
        if not annotation_types:
            annotation_types = [("None", "None", "None")]
        return annotation_types

    registered_classes = []

    def register(cls):
        bpy.utils.register_class(cls)
        registered_classes.append(cls)

    attributes = {"__annotations__": {}}
    attributes["__annotations__"]["type"] = EnumProperty(items=get_annotation_types)
    for cls in entity.annotations._classes.values():
        AnnotationTypeInputs = create_annotation_type_inputs(cls)
        register(AnnotationTypeInputs)
        entity_annotation_type = (
            f"{entity._get_annotation_entity_type()}_{cls.annotation_type}"
        )
        attributes["__annotations__"][entity_annotation_type] = (
            bpy.props.PointerProperty(type=AnnotationTypeInputs)
        )
    AnnotationProps = type("AnnotationProps", (bpy.types.PropertyGroup,), attributes)
    register(AnnotationProps)

    # Temporary annotation add operator
    class TempAnnotationAddOperator(bpy.types.Operator):
        bl_idname = "mn.temp_annotation_add"
        bl_label = "Add annotation"
        bl_description = "Add a new annotation"

        props: bpy.props.PointerProperty(type=AnnotationProps)  # type: ignore

        def draw(self, context):
            layout = self.layout
            layout.prop(self.props, "type")
            entity_annotation_type = (
                f"{entity._get_annotation_entity_type()}_{self.props.type}"
            )
            inputs = getattr(self.props, entity_annotation_type, None)
            if inputs is not None:
                for prop_name in inputs.__annotations__.keys():
                    if prop_name == "uuid":
                        continue
                    layout.prop(inputs, prop_name)

        def invoke(self, context, event):
            return context.window_manager.invoke_props_dialog(self)

        def execute(self, context):
            if self.props.type == "None":
                return {"CANCELLED"}
            annotation_class = entity.annotations._classes[self.props.type]
            api_inputs = {}
            entity_annotation_type = (
                f"{entity._get_annotation_entity_type()}_{self.props.type}"
            )
            ui_inputs = getattr(self.props, entity_annotation_type, None)
            if ui_inputs is not None:
                for prop_name in ui_inputs.__annotations__.keys():
                    if prop_name in ui_inputs:
                        api_inputs[prop_name] = getattr(ui_inputs, prop_name)
                    else:
                        if hasattr(annotation_class, prop_name):
                            api_inputs[prop_name] = getattr(annotation_class, prop_name)
            # call the same method used by the APIs
            method_name = f"add_{annotation_class.annotation_type}"
            method = getattr(entity.annotations, method_name)
            method(**api_inputs)
            return {"FINISHED"}

    register(TempAnnotationAddOperator)
    return registered_classes


class MN_OT_Add_Annotation(Operator):
    """
    Operator to add a new annotation to an entity
    """

    bl_idname = "mn.add_annotation"
    bl_label = "Add annotation"
    bl_description = "Add annotation to entity"

    uuid: StringProperty()  # type: ignore

    _temp_classes = []

    def draw(self, context):
        layout = self.layout
        assert layout
        layout.prop(self, "type")

    def execute(self, context: Context):
        # unregister any temp classes from previous invocation
        # verify in: bpy.types.Operator.__subclasses__() and
        # bpy.types.PropertyGroup.__subclasses__()
        for cls in MN_OT_Add_Annotation._temp_classes:
            bpy.utils.unregister_class(cls)
        entity = get_session().get(self.uuid)
        # register the temporary operator with required type inputs
        MN_OT_Add_Annotation._temp_classes = _register_temp_annotation_add_op(entity)
        # invoke the temporary add operator
        bpy.ops.mn.temp_annotation_add("INVOKE_DEFAULT")
        return {"FINISHED"}


class MN_OT_Remove_Annotation(Operator):
    """
    Operator to remove an annotation from an entity
    """

    bl_idname = "mn.remove_annotation"
    bl_label = "Remove annotation"
    bl_description = "Remove annotation from entity"

    uuid: StringProperty()  # type: ignore
    annotation_uuid: StringProperty()  # type: ignore

    def execute(self, context: Context):
        entity = get_session().get(self.uuid)
        entity.annotations._remove_annotation_by_uuid(self.annotation_uuid)
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_confirm(
            self, event, title="Remove Annotation?"
        )


class MN_OT_Setup_Compositor(Operator):
    """
    Operator to setup compositor
    """

    bl_idname = "mn.setup_compositor"
    bl_label = "Setup Compositor"
    bl_description = "Setup Molecular Nodes Compositor"

    def execute(self, context: Context):
        setup_compositor(context.scene)
        return {"FINISHED"}


class MN_OT_DSSP_init(Operator):
    """
    Operator to initialize DSSP for trajectories
    """

    bl_idname = "mn.dssp_init"
    bl_label = "Initialize"
    bl_description = "Initialize DSSP analysis for trajectory"

    uuid: StringProperty()  # type: ignore

    def execute(self, context: Context):
        entity = get_session().get(self.uuid)
        if entity is None:
            return {"CANCELLED"}
        entity.dssp.init()
        return {"FINISHED"}


class MN_OT_DSSP_apply(Operator):
    """
    Operator to apply changed DSSP options
    """

    bl_idname = "mn.dssp_apply"
    bl_label = "Apply"
    bl_description = "Apply changed DSSP options"

    uuid: StringProperty()  # type: ignore
    apply_ta_threshold: BoolProperty()  # type: ignore
    ta_threshold: FloatProperty()  # type: ignore

    def execute(self, context: Context):
        entity = get_session().get(self.uuid)
        if entity is None:
            return {"CANCELLED"}
        props = entity.object.mn.dssp
        if props.display_option == "trajectory-average":
            if self.apply_ta_threshold:
                entity.dssp.show_trajectory_average(threshold=self.ta_threshold)
            else:
                entity.dssp.show_trajectory_average()
        return {"FINISHED"}


class MN_OT_DSSP_cancel(Operator):
    """
    Operator to cancel and restore current DSSP options
    """

    bl_idname = "mn.dssp_cancel"
    bl_label = "Cancel"
    bl_description = "Restore current DSSP options"

    uuid: StringProperty()  # type: ignore

    def execute(self, context: Context):
        entity = get_session().get(self.uuid)
        if entity is None:
            return {"CANCELLED"}
        props = entity.object.mn.dssp
        props.cancelling = True
        props.display_option = entity.dssp._display_option
        props.window_size = entity.dssp._window_size
        if entity.dssp._sw_threshold is not None:
            props.sw_threshold = entity.dssp._sw_threshold
            props.apply_sw_threshold = True
        else:
            props.apply_sw_threshold = False
        if entity.dssp._ta_threshold is not None:
            props.ta_threshold = entity.dssp._ta_threshold
            props.apply_ta_threshold = True
        else:
            props.apply_ta_threshold = False
        props.applied = True
        props.cancelling = False
        return {"FINISHED"}


CLASSES = [
    MN_OT_Import_Fetch,
    MN_OT_Import_Trajectory,
    MN_OT_Reload_Trajectory,
    MN_OT_Import_Map,
    MN_OT_Import_Ensemble,
    MN_OT_Import_Molecule,
    MN_FH_Import_Molecule,
    MN_OT_Add_Style,
    MN_OT_Remove_Style,
    MN_OT_Add_Annotation,
    MN_OT_Remove_Annotation,
    MN_OT_Setup_Compositor,
    MN_OT_DSSP_init,
    MN_OT_DSSP_apply,
    MN_OT_DSSP_cancel,
]
