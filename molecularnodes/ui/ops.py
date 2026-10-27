import gc
import os
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Annotated, Any, cast
import bpy
import MDAnalysis as mda
import nodebpy
from biotite import InvalidFileError
from bpy.props import (
    BoolProperty,
    CollectionProperty,
    EnumProperty,
    FloatProperty,
    IntProperty,
    StringProperty,
)
from bpy.types import Context, Operator
from ..annotations.props import create_annotation_type_inputs
from ..blender import coll
from ..blender.registered import object_props
from ..blender.utils import annotated_props, path_resolve, resolve_file_path
from ..download import CACHE_DIR, FileDownloadPDBError
from ..entities import (
    OXDNA,
    Molecule,
    StreamingTrajectory,
    density,
    ensemble,
    molecule,
)
from ..entities.molecule.pdbx import MultipleStructuresError
from ..handlers import update_entities
from ..nodes._utils import STYLE_LITERALS, remove_style_node, styles_mapping, swap
from ..scene.compositor import setup_compositor
from ..session import get_session
from ..utils import _increase_view_distance
from .pref import addon_preferences
from .style import STYLE_ITEMS

# MD file formats accepted when dropped into the viewport. Trajectory files hold
# only coordinates, so they are paired with the single other (topology) file
# dropped alongside them; topology files dropped on their own are read with
# MDAnalysis when they carry coordinates (e.g. a `.gro`)
MD_TOPOLOGY_EXTENSIONS = {".gro", ".tpr", ".psf", ".prmtop", ".parm7", ".data"}
MD_TRAJECTORY_EXTENSIONS = {
    ".xtc",
    ".trr",
    ".dcd",
    ".nc",
    ".ncdf",
    ".mdcrd",
    ".lammpstrj",
    ".h5md",
}


def _natural_key(path: Path) -> list:
    """Sort key ordering `md_2.xtc` before `md_10.xtc`."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", path.name)]


class MN_FH_Import_Molecule(bpy.types.FileHandler):
    bl_idname = "MN_FH_import_molecule"
    bl_label = "File handler for import molecular data files."
    bl_import_operator = "mn.import_molecule"
    bl_file_extensions = ";".join(
        [".pdb", ".cif", ".mmcif", ".bcif", ".pdbx", ".xyz"]
        + sorted(MD_TOPOLOGY_EXTENSIONS | MD_TRAJECTORY_EXTENSIONS)
    )

    @classmethod
    def poll_drop(cls, context):
        return context.area and context.area.type == "VIEW_3D"


DOWNLOAD_FORMATS = (
    ("bcif", ".bcif", "Binary compressed .cif file, fastest for downloading"),
    ("cif", ".cif", "The new standard of .cif / .mmcif"),
    ("pdb", ".pdb", "The classic (and depcrecated) PDB format"),
)

# structure-file extensions `Molecule.load()` can parse, used when scanning a
# folder; MD topology + trajectory pairs are excluded, as which trajectory
# belongs to which topology can't be matched up automatically
FOLDER_IMPORT_EXTENSIONS = {".pdb", ".cif", ".bcif", ".xyz", ".sdf", ".mol"}


def _load_structures(path: Path) -> list[Molecule]:
    """Load a structure file, as one Molecule per structure when a .cif/.bcif
    holds several structures in separate data blocks."""
    try:
        return [Molecule.load(path)]
    except MultipleStructuresError as e:
        return [
            Molecule.load(path, data_block=block, name=f"{path.name}:{block}")
            for block in e.blocks
        ]


# operator that is called by the 'button' press which calls the fetch function


@annotated_props
class MN_OT_Import_Molecule(bpy.types.Operator):
    bl_idname = "mn.import_molecule"
    bl_label = "Import Molecule"
    bl_description = (
        "Open a local structure file or MD trajectory, or download a structure "
        "from a database"
    )
    bl_options = {"REGISTER", "UNDO"}

    style: Annotated[
        STYLE_LITERALS,
        EnumProperty(
            name="Style",
            default="spheres",
            description="Starting style for the structure on import",
            items=STYLE_ITEMS,
        ),
    ]
    node_setup: Annotated[
        bool,
        BoolProperty(
            name="Node Setup",
            default=True,
            description="Whether to setup the starting default node tree on import",
        ),
    ]
    assembly: Annotated[
        bool,
        BoolProperty(
            default=False,
            name="Build Biological Assembly",
            description="Build the biological assembly for the structure on import",
        ),
    ]
    share_node_group: Annotated[
        bool,
        BoolProperty(
            default=False,
            name="Share Node Group",
            description=(
                "When importing multiple structures, style them all with the same "
                "node group instead of creating one per structure, so tweaks apply "
                "to every structure at once (useful for conformations of one protein)"
            ),
        ),
    ]
    # filled in by the file handler when structure files are dropped into the 3D
    # viewport, in place of the dialog's own filepath field
    directory: Annotated[
        str, StringProperty(subtype="FILE_PATH", options={"SKIP_SAVE", "HIDDEN"})
    ]
    files: Annotated[
        "bpy.types.bpy_prop_collection_idprop[bpy.types.OperatorFileListElement]",
        CollectionProperty(
            type=bpy.types.OperatorFileListElement, options={"SKIP_SAVE", "HIDDEN"}
        ),
    ]

    method: Annotated[
        str,
        EnumProperty(
            name="Database",
            default="fetch",
            items=(
                (
                    "fetch",
                    "Fetch",
                    "Download a structure from a database",
                ),
                (
                    "local",
                    "Local File",
                    "Open a structure file or MD topology + trajectory already on disk",
                ),
            ),
        ),
    ]
    database: Annotated[
        str,
        EnumProperty(
            name="Database",
            default="wwpdb",
            items=(
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
                (
                    "pubchem",
                    "PubChem",
                    "3D conformers of small molecules from PubChem, by CID or name",
                ),
            ),
        ),
    ]
    code: Annotated[
        str,
        StringProperty(
            name="PDB",
            description=(
                "The code to download (4-character e.g. '1abc' or 12-character e.g."
                " 'pdb_00001abc'). Separate multiple codes with commas to download"
                " several structures at once"
            ),
            options={"TEXTEDIT_UPDATE"},
        ),
    ]
    filepath: Annotated[
        str,
        StringProperty(
            name="File",
            description=(
                "Path of the local structure (or MD topology) file to open, or a "
                "folder to recursively import every structure file inside of it"
            ),
            subtype="FILE_PATH",
        ),
    ]
    recursive: Annotated[
        bool,
        BoolProperty(
            name="Recursive",
            default=True,
            description=(
                "Import structure files from all subfolders too, rather than just "
                "the top level of the folder"
            ),
        ),
    ]
    objects_only: Annotated[
        bool,
        BoolProperty(
            name="Objects Only",
            default=False,
            description=(
                "Discard the Python `Molecule` (and its `mda.Universe`) after each "
                "import, keeping only the created Blender object - saves memory when "
                "importing many structures. Objects can be relinked later with the "
                "session Reload operator"
            ),
        ),
    ]
    trajectory: Annotated[
        str,
        StringProperty(
            name="Trajectory",
            description=(
                "Optional trajectory file to load alongside the topology, "
                "or an imd:// URL to stream frames from a running simulation"
            ),
            subtype="FILE_PATH",
        ),
    ]
    additional_arguments: Annotated[
        str,
        StringProperty(
            name="Arguments",
            description="Additional arguments to pass to the `mda.Universe(topology, trajectory, **kwargs)` constructor",
            default="",
        ),
    ]
    file_format: Annotated[
        str,
        EnumProperty(
            name="Format",
            description="Format to download as from the PDB",
            default="bcif",
            items=DOWNLOAD_FORMATS,
        ),
    ]

    cache_dir: Annotated[
        str,
        StringProperty(
            name="Cache Directory",
            description="Where to store the structures downloaded from the Protein Data Bank",
            default=str(CACHE_DIR),
            subtype="DIR_PATH",
        ),
    ]

    def apply_import_options(self, mol: Molecule) -> None:
        """Set up the default node tree on a freshly imported molecule."""
        if not self.node_setup:
            return
        mol.create_asset_nodes()
        kwargs = {}
        if self.style == "spheres":
            # the API picks instanced spheres for non-Cycles engines so renders
            # look right, but a GUI import has to stay interactive on large
            # systems, so it always keeps the point cloud (#1220)
            kwargs["sphere"] = "Point"
        mol.add_style(
            style=self.style, color="common", assembly=self.assembly, **kwargs
        )

    def draw(self, context):
        layout = self.layout
        assert layout
        show_assembly = True
        if self.files:
            topologies, trajectories = self._split_dropped_files()
            if trajectories:
                show_assembly = False
                layout.label(
                    text=f"Importing {len(trajectories)} trajectory files",
                    icon="FILE_MOVIE",
                )
                if len(topologies) != 1:
                    layout.label(
                        text="Drop exactly one topology file with the trajectory",
                        icon="ERROR",
                    )
                layout.prop(self, "additional_arguments")
            else:
                layout.label(text=f"Importing {len(self.files)} molecules")
            if len(self.files) > 1 and not trajectories:
                layout.prop(self, "share_node_group")
                layout.prop(self, "objects_only")
        else:
            layout.prop_tabs_enum(self, "method")
            if self.method == "local":
                layout.prop(self, "filepath")
                folder = self._folder_path()
                if folder is not None:
                    layout.prop(self, "recursive")
                    files, complete = self._scan_folder(folder, max_entries=10_000)
                    if not files and complete:
                        layout.label(
                            text="No structure files found in the folder",
                            icon="ERROR",
                        )
                    else:
                        count = f"{len(files)}" if complete else f"{len(files)}+"
                        layout.label(
                            text=f"Importing {count} structure files from the folder",
                            icon="FILE_FOLDER",
                        )
                        layout.prop(self, "share_node_group")
                        layout.prop(self, "objects_only")
                else:
                    layout.prop(self, "trajectory")
                    layout.prop(self, "additional_arguments")
            else:
                layout.prop_tabs_enum(self, "database")
                row = layout.row().split(factor=0.7)
                row.prop(self, "code")
                # file format only applies to wwPDB downloads; others pick their own
                if self.method == "fetch" and self.database != "pubchem":
                    row.prop(self, "file_format", text="")
                if "," in self.code:
                    layout.prop(self, "share_node_group")
        row = layout.row()
        row.prop(self, "node_setup", text="")
        col = row.column()
        col.prop(self, "style")
        col.enabled = self.node_setup
        if show_assembly and not self.trajectory:
            layout.prop(self, "assembly")

    def invoke(self, context, event):
        prefs = addon_preferences()
        self.cache_dir = str(prefs.cache_dir) if prefs is not None else bpy.app.tempdir
        return context.window_manager.invoke_props_dialog(self)

    def _folder_path(self) -> Path | None:
        """The dialog's filepath as a directory, when it points at one."""
        if self.method != "local" or not self.filepath:
            return None
        path = path_resolve(self.filepath)
        return path if path.is_dir() else None

    def _scan_folder(
        self, folder: Path, max_entries: int | None = None
    ) -> tuple[list[Path], bool]:
        """The structure files under `folder`, and whether the scan was complete.

        `max_entries` bounds how many directory entries are visited, so the
        dialog can preview the count without walking an unexpectedly huge tree
        (e.g. a half-typed path that resolves to the home directory); a scan
        that hits the bound is returned as incomplete. Results are cached per
        folder (and recursive setting), as `draw` runs on every UI redraw."""
        cached = getattr(self, "_folder_scan", None)
        if cached is not None and cached[0] == (folder, self.recursive):
            _, files, complete = cached
            # a cached complete scan answers everything; an incomplete one is
            # only enough when a bounded preview is being asked for again
            if complete or max_entries is not None:
                return files, complete
        if self.recursive:
            walker = os.walk(folder)
        else:
            walker = [
                (
                    str(folder),
                    [],
                    [entry.name for entry in os.scandir(folder) if entry.is_file()],
                )
            ]
        files = []
        complete = True
        visited = 0
        for root, _dirs, names in walker:
            for name in names:
                visited += 1
                if max_entries is not None and visited > max_entries:
                    complete = False
                    break
                if Path(name).suffix.lower() in FOLDER_IMPORT_EXTENSIONS:
                    files.append(Path(root, name))
            if not complete:
                break
        files.sort()
        self._folder_scan = ((folder, self.recursive), files, complete)
        return files, complete

    def _universe_kwargs(self) -> dict:
        """Parse the additional arguments string into `mda.Universe` kwargs."""
        if self.additional_arguments == "":
            return {}
        try:
            return {
                key.strip(): arg.strip()
                for key, arg in [
                    kv.split("=") for kv in self.additional_arguments.split(",")
                ]
            }
        except Exception as e:
            self.report(
                {"WARNING"},
                message=f"Failed to parse additional arguments: {e}",
            )
            return {}

    def _setup_molecule(self, mol, shared_tree=None):
        """Apply the import options to a freshly imported molecule, re-using
        `shared_tree` for the styling instead when one is given. Returns the
        tree to share with subsequently imported molecules."""
        if shared_tree is not None:
            # re-use the first structure's node group rather than creating a
            # separate styling tree per structure
            if self.node_setup:
                mol.create_asset_nodes()
            mol.tree = shared_tree
            return shared_tree
        self.apply_import_options(mol)
        return mol.tree if self.share_node_group else None

    def execute(self, context):
        if self.files:
            return self._execute_dropped_files(context)

        folder = self._folder_path()
        if folder is not None:
            return self._execute_folder(context, folder)

        mols = []
        try:
            if self.method == "local":
                topology = resolve_file_path(self.filepath, "Structure file")
                if self.trajectory.startswith("imd://"):
                    mol = StreamingTrajectory.load(
                        topology=topology,
                        coordinates=self.trajectory,
                        name=topology.name,
                        style=None,
                    )
                    message = (
                        f"Streaming trajectory '{mol.name}' from '{self.trajectory}'"
                    )
                elif self.trajectory:
                    mol = Molecule.load(
                        topology,
                        resolve_file_path(self.trajectory, "Trajectory file"),
                        **self._universe_kwargs(),
                    )
                    message = (
                        f"Imported '{self.filepath}' as {mol.name} with "
                        f"{int(mol.props.n_frames)} frames from '{self.trajectory}'."
                    )
                else:
                    mols = _load_structures(topology)
                    mol = mols[-1]
                    if len(mols) == 1:
                        message = f"Imported '{self.filepath}' as {mol.name}"
                    else:
                        message = (
                            f"Imported {len(mols)} structures from '{self.filepath}'"
                        )
                if not mols:
                    mols = [mol]
            else:
                codes = [code.strip() for code in self.code.split(",") if code.strip()]
                if not codes:
                    self.report({"ERROR"}, "No code given to fetch")
                    return {"CANCELLED"}
                mols = [
                    Molecule.fetch(
                        code=code,
                        cache=self.cache_dir,
                        format=self.file_format,
                        database=self.database,
                    )
                    for code in codes
                ]
                mol = mols[-1]
                if len(mols) == 1:
                    message = f"Downloaded {codes[0]} as {mol.name}"
                else:
                    message = f"Downloaded {len(mols)} structures: {', '.join(codes)}"
        except FileDownloadPDBError as e:
            self.report({"ERROR"}, str(e))
            if self.file_format == "pdb" and self.database != "pubchem":
                self.report(
                    {"ERROR"},
                    "There may not be a `.pdb` formatted file available - try a different download format.",
                )
            return {"CANCELLED"}
        except (InvalidFileError, ValueError, OSError) as e:
            # a missing file, or unreadable file contents (e.g. a structure
            # factor file with no coordinates); show the message instead of a
            # console traceback
            self.report({"ERROR"}, str(e))
            return {"CANCELLED"}

        shared_tree = None
        for imported in mols:
            shared_tree = self._setup_molecule(imported, shared_tree)

        if isinstance(mol, StreamingTrajectory):
            context.scene.frame_start = 0
        else:
            n_frames = int(mol.props.n_frames)
            if n_frames > 1:
                context.scene.frame_start = 0
                context.scene.frame_end = n_frames

        try:
            context.view_layer.objects.active = mol.object
        except RuntimeError:
            message += " - Molecular Nodes collection is disabled"

        self.report({"INFO"}, message=message)

        _increase_view_distance()
        return {"FINISHED"}

    def _split_dropped_files(self) -> tuple[list[Path], list[Path]]:
        """The dropped files split into topology/structure files and MD
        trajectory files, the trajectories in natural filename order."""
        paths = [Path(self.directory, file.name) for file in self.files]
        trajectories = sorted(
            (p for p in paths if p.suffix.lower() in MD_TRAJECTORY_EXTENSIONS),
            key=_natural_key,
        )
        topologies = [p for p in paths if p not in trajectories]
        return topologies, trajectories

    def _execute_dropped_trajectory(
        self, context, topology: Path, trajectories: list[Path]
    ):
        """Import a dropped topology with one or more trajectory files, which
        MDAnalysis chains together into a single trajectory."""
        try:
            mol = Molecule.load(
                topology,
                trajectories[0] if len(trajectories) == 1 else trajectories,
                **self._universe_kwargs(),
            )
        except Exception as e:
            self.report({"ERROR"}, f"Failed importing {topology.name}: {e}")
            return {"CANCELLED"}
        self.apply_import_options(mol)

        n_frames = int(mol.props.n_frames)
        if n_frames > 1:
            context.scene.frame_start = 0
            context.scene.frame_end = n_frames
        self.report(
            {"INFO"},
            message=(
                f"Imported '{topology.name}' as {mol.name} with {n_frames} frames "
                f"from {len(trajectories)} trajectory files"
            ),
        )
        _increase_view_distance()
        return {"FINISHED"}

    def _execute_dropped_files(self, context):
        """Import each structure file dropped into the viewport, or a topology
        with its trajectory files when any are dropped."""
        topologies, trajectories = self._split_dropped_files()
        if trajectories:
            if len(topologies) != 1:
                names = ", ".join(p.name for p in topologies) or "none"
                self.report(
                    {"ERROR"},
                    "Trajectory files need exactly one topology file dropped "
                    f"with them (got: {names})",
                )
                return {"CANCELLED"}
            return self._execute_dropped_trajectory(
                context, topologies[0], trajectories
            )

        imported = 0
        shared_tree = None
        session = get_session(context)
        for path in topologies:
            try:
                if path.suffix.lower() in MD_TOPOLOGY_EXTENSIONS:
                    # an MD topology with its own coordinates (e.g. a `.gro`)
                    mols = [Molecule(mda.Universe(path), name=path.name)]
                else:
                    mols = _load_structures(path)
                for mol in mols:
                    shared_tree = self._setup_molecule(mol, shared_tree)
                    if self.objects_only:
                        session.remove_entity(mol.uuid)
                    imported += 1
            except Exception as e:
                self.report({"WARNING"}, message=f"Failed importing {path.name}: {e}")

        if imported == 0:
            return {"CANCELLED"}

        self.report({"INFO"}, message=f"Imported {imported} molecules")
        _increase_view_distance()
        return {"FINISHED"}

    def _execute_folder(self, context, folder: Path):
        """Import every structure file under `folder` (recursively unless turned
        off), mirroring its subfolder layout as nested collections that the
        objects are placed in."""
        files, _ = self._scan_folder(folder)
        if not files:
            self.report({"ERROR"}, f"No structure files found under '{folder}'")
            return {"CANCELLED"}

        collections: dict[Path, bpy.types.Collection] = {}

        def collection_for(directory: Path) -> bpy.types.Collection:
            if directory in collections:
                return collections[directory]
            parent = (
                coll.mn() if directory == folder else collection_for(directory.parent)
            )
            # reuse a matching collection already under this parent (a repeated
            # import), but never one from elsewhere in the scene - collection
            # names are global in Blender, so same-named subfolders in different
            # branches get a .001-suffixed collection rather than being merged
            child = next((c for c in parent.children if c.name == directory.name), None)
            if child is None:
                child = bpy.data.collections.new(directory.name)
                parent.children.link(child)
            collections[directory] = child
            return child

        session = get_session(context)
        window_manager = context.window_manager
        window_manager.progress_begin(0, len(files))
        imported = 0
        failed = 0
        shared_tree = None
        try:
            for i, path in enumerate(files):
                window_manager.progress_update(i)
                try:
                    mols = _load_structures(path)
                    for mol in mols:
                        shared_tree = self._setup_molecule(mol, shared_tree)
                except Exception as e:
                    failed += 1
                    self.report(
                        {"WARNING"},
                        f"Failed importing {path.relative_to(folder)}: {e}",
                    )
                    continue
                for mol in mols:
                    obj = mol.object
                    for c in obj.users_collection:
                        c.objects.unlink(obj)
                    collection_for(path.parent).objects.link(obj)
                    imported += 1
                    if self.objects_only:
                        session.remove_entity(mol.uuid)
                        del mol
                        # the entity and its managers reference each other, so the
                        # universes only free on a cycle collection - don't let
                        # hundreds of them pile up while the import runs
                        if imported % 50 == 0:
                            gc.collect()
                del mols
        finally:
            window_manager.progress_end()
            if self.objects_only:
                gc.collect()

        if imported == 0:
            return {"CANCELLED"}
        message = f"Imported {imported} structures from '{folder.name}'"
        if failed:
            message += f", {failed} failed"
        self.report({"INFO"}, message)
        _increase_view_distance()
        return {"FINISHED"}


ENSEMBLE_TYPES = (
    (
        "starfile",
        "Starfile",
        "Import a .star mapback file or CryoET .ndjson point annotations",
    ),
    ("cellpack", "CellPack", "Import a CellPack .cif / .bcif model"),
)


@annotated_props
class MN_OT_Import_Ensemble(bpy.types.Operator):
    bl_idname = "mn.import_ensemble"
    bl_label = "Import Ensemble"
    bl_description = "Import an ensemble as a Starfile or CellPack model"
    bl_options = {"REGISTER", "UNDO"}

    ensemble_type: Annotated[
        str,
        EnumProperty(
            name="Type",
            description="The kind of ensemble to import",
            default="starfile",
            items=ENSEMBLE_TYPES,
        ),
    ]
    filepath: Annotated[
        str,
        StringProperty(
            name="File",
            description="File path for the ensemble to import",
            subtype="FILE_PATH",
            maxlen=0,
        ),
    ]
    node_setup: Annotated[
        bool,
        BoolProperty(
            name="Setup Nodes",
            default=True,
            description="Create and set up a Geometry Nodes tree on import",
        ),
    ]

    def draw(self, context):
        layout = self.layout
        assert layout
        layout.prop_tabs_enum(self, "ensemble_type")
        layout.prop(self, "filepath")
        layout.prop(self, "node_setup")

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)

    def execute(self, context):
        try:
            file_path = resolve_file_path(self.filepath, "Ensemble file")
            if self.ensemble_type == "cellpack":
                ensemble.CellPack.load(
                    file_path=file_path,
                    name=file_path.name,
                    node_setup=self.node_setup,
                )
            else:
                ensemble.StarFile.load(
                    file_path=file_path,
                    node_setup=self.node_setup,
                )
        except (ValueError, OSError) as e:
            self.report({"ERROR"}, str(e))
            return {"CANCELLED"}

        _increase_view_distance()
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


@annotated_props
class MN_OT_Import_Map(bpy.types.Operator):
    bl_idname = "mn.import_density"
    bl_label = "Import Density"
    bl_description = "Import an EM density map into Blender"
    bl_options = {"REGISTER", "UNDO"}

    filepath: Annotated[
        str,
        StringProperty(
            name="File",
            description="File path for the map file.",
            subtype="FILE_PATH",
            maxlen=0,
        ),
    ]
    invert: Annotated[
        bool,
        BoolProperty(
            name="Invert Data",
            description="Invert the values in the map. Low becomes high, high becomes low.",
            default=False,
        ),
    ]
    center: Annotated[
        bool,
        BoolProperty(
            name="Center Density",
            description="Translate the density so that the center of the box is at the origin.",
            default=False,
        ),
    ]
    overwrite: Annotated[
        bool,
        BoolProperty(
            name="Overwrite Intermediate File",
            description="Overwrite generated intermediate .vdb file.",
            default=False,
        ),
    ]
    setup_nodes: Annotated[
        bool,
        BoolProperty(
            name="Setup Nodes",
            default=True,
            description="Create and set up a Geometry Nodes tree on import",
        ),
    ]
    style: Annotated[
        str,
        EnumProperty(
            name="Style",
            items=DENSITY_STYLE_ITEMS,
        ),
    ]

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
        try:
            density.Grids.load(
                file_path=resolve_file_path(self.filepath, "Map file"),
                invert=self.invert,
                style=self.style if self.setup_nodes else None,
                center=self.center,
                overwrite=self.overwrite,
            )
        except (ValueError, OSError) as e:
            self.report({"ERROR"}, str(e))
            return {"CANCELLED"}
        _increase_view_distance()
        return {"FINISHED"}


@annotated_props
class MN_OT_Frames_To_Collection(bpy.types.Operator):
    bl_idname = "mn.frames_to_collection"
    bl_label = "Bake Frames to Collection"
    bl_description = (
        "Bake a range of trajectory frames into a collection of objects (one per frame) "
        "with the positions and per-frame attributes, which can be read inside geometry "
        "nodes with the Animate Frames node"
    )
    bl_options = {"REGISTER", "UNDO"}

    start: Annotated[
        int,
        IntProperty(
            name="Start",
            default=0,
            min=0,
            description="First trajectory frame to bake",
        ),
    ]
    stop: Annotated[
        int,
        IntProperty(
            name="Stop",
            default=0,
            min=0,
            description="One past the last frame to bake (0 bakes through the final frame)",
        ),
    ]
    step: Annotated[
        int,
        IntProperty(
            name="Step",
            default=1,
            min=1,
            description="Stride between baked frames",
        ),
    ]

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        if obj is None:
            return False
        entity = context.scene.MNSession.match(obj)
        return isinstance(entity, Molecule) and obj.mn.n_frames > 1

    def invoke(self, context, event):
        # default the stop to the trajectory length for convenience
        if self.stop == 0:
            self.stop = context.active_object.mn.n_frames
        return context.window_manager.invoke_props_dialog(self)

    def execute(self, context):
        obj = context.active_object
        entity = context.scene.MNSession.match(obj)
        if entity is None:
            self.report({"ERROR"}, "Active object is not linked to a molecule")
            return {"CANCELLED"}

        stop = self.stop if self.stop > 0 else None
        try:
            frames = entity.frames_to_collection(
                start=self.start, stop=stop, step=self.step
            )
        except ValueError as e:
            self.report({"ERROR"}, str(e))
            return {"CANCELLED"}

        self.report(
            {"INFO"},
            f"Baked {len(frames.objects)} frames into collection '{frames.name}'",
        )
        return {"FINISHED"}


@annotated_props
class MN_OT_Import_OxDNA(bpy.types.Operator):
    bl_idname = "mn.import_oxdna"
    bl_label = "Import oxDNA"
    bl_description = "Load an oxDNA topology and trajectory"
    bl_options = {"REGISTER", "UNDO"}

    topology: Annotated[
        str,
        StringProperty(
            name="Topology",
            description="File path for the topology file for the trajectory",
            subtype="FILE_PATH",
            maxlen=0,
        ),
    ]
    trajectory: Annotated[
        str,
        StringProperty(
            name="Trajectory",
            description="File path for the trajectory file for the trajectory",
            subtype="FILE_PATH",
            maxlen=0,
        ),
    ]
    name: Annotated[
        str,
        StringProperty(
            name="Name",
            description="Name of the molecule on import",
            default="oxDNA",
            maxlen=0,
        ),
    ]

    def draw(self, context):
        layout = self.layout
        assert layout
        layout.prop(self, "name")
        layout.prop(self, "topology")
        layout.prop(self, "trajectory")

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)

    def execute(self, context):
        try:
            OXDNA.load(
                topology=resolve_file_path(self.topology, "Topology file"),
                coordinates=resolve_file_path(self.trajectory, "Trajectory file"),
                name=self.name,
            )
        except (ValueError, OSError) as e:
            self.report({"ERROR"}, str(e))
            return {"CANCELLED"}

        _increase_view_distance()
        return {"FINISHED"}


@annotated_props
class MN_OT_add_selection_to_style(Operator):
    """
    Create a new selection and add the corresponding `Named Attribute` node to the node tree's style
    """

    bl_idname = "mn.add_selection_to_style"
    bl_label = "Add Style Selection"
    bl_description = ""

    node_tree: Annotated[str, StringProperty()]
    node_name: Annotated[str, StringProperty()]

    def execute(self, context: bpy.types.Context):
        obj = context.active_object
        assert obj
        session = get_session(context)
        mol = cast(Molecule, session.match(obj))

        if mol is None:
            self.report({"ERROR"}, message="No molecule found for this object")
            return {"CANCELLED"}

        node_group = bpy.data.node_groups.get(self.node_tree)
        if node_group is None:
            self.report({"ERROR"}, message=f"Node tree '{self.node_tree}' not found")
            return {"CANCELLED"}

        node = node_group.nodes.get(self.node_name)
        if node is None:
            self.report(
                {"ERROR"},
                message=f"Node '{self.node_name}' not found in '{self.node_tree}'",
            )
            return {"CANCELLED"}

        selection = mol.selections.from_string("all")

        with nodebpy.TreeBuilder(node_group):
            style = nodebpy.geometry.Group._from_node(node)
            selection.node() >> style.i["Selection"]

        return {"FINISHED"}


@annotated_props
class MN_OT_Add_Style(Operator):
    """
    Operator to add a new style to an entity
    """

    bl_idname = "mn.add_style"
    bl_label = "Add Style"
    bl_description = "Add a new style to the entity"
    bl_options = {"REGISTER", "UNDO"}

    uuid: Annotated[str, StringProperty()]
    # fallback identifier for objects not linked to a session entity, whose
    # style branch is built directly in the object's node tree
    name_object: Annotated[str, StringProperty()]

    style: Annotated[
        STYLE_LITERALS,
        EnumProperty(
            name="Style",
            default="spheres",
            description="Style type",
            items=STYLE_ITEMS,
        ),
    ]

    use_uniform_color: Annotated[
        bool,
        BoolProperty(
            name="Use uniform color",
            description="Use uniform color for style",
            default=False,
        ),
    ]

    # a bpy_prop_array at runtime, which the stubs don't type as a sequence
    uniform_color: Annotated[
        Sequence[float],
        bpy.props.FloatVectorProperty(
            name="Uniform color",
            description="Uniform color for style",
            subtype="COLOR",
            size=4,
            default=(0.162, 0.624, 0.196, 1),
            min=0.0,
            max=1.0,
        ),
    ]

    color_scheme: Annotated[
        str,
        StringProperty(
            name="Coloring scheme",
            description="Coloring scheme for style",
            default="common",
        ),
    ]

    selection: Annotated[
        str,
        StringProperty(
            name="Selection",
            description="Selection for which the style applies",
            default="all",
        ),
    ]

    name: Annotated[
        str,
        StringProperty(
            name="Name",
            description="Label for the style",
            default="",
        ),
    ]

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
        selection = self.selection.strip() or None
        name = self.name.strip() or None
        if entity is not None:
            # the Styles panel only enables adding for molecules
            cast(Molecule, entity).add_style(
                style=self.style,
                color=self.uniform_color
                if self.use_uniform_color
                else self.color_scheme,
                selection=selection,
                name=name,
            )
            return {"FINISHED"}

        # without a session entity the style branch is built directly in the
        # object's node tree — no universe means no MDAnalysis selections and
        # no color schemes, so those degrade with a warning
        obj = bpy.data.objects.get(self.name_object)
        mod = obj.modifiers.get("Molecular Nodes") if obj is not None else None
        tree = getattr(mod, "node_group", None)
        if obj is None or tree is None:
            self.report({"ERROR"}, "Object has no Molecular Nodes node tree")
            return {"CANCELLED"}

        attribute = None
        if selection is not None and selection != "all":
            if obj.data is not None and selection in obj.data.attributes:
                attribute = selection
            else:
                self.report(
                    {"WARNING"},
                    f"'{selection}' is not an attribute on the object; selections "
                    "need the entity reloaded (linked) to be evaluated. Style "
                    "added to all atoms",
                )
        if not self.use_uniform_color and self.color_scheme not in (
            "common",
            "default",
        ):
            # baked attribute colors already reflect the common scheme; anything
            # else needs the entity to compute
            self.report(
                {"WARNING"},
                f"Color scheme '{self.color_scheme}' needs a linked entity; "
                "keeping the existing colors",
            )
        molecule.base.add_style_to_tree(
            tree,
            style=self.style,
            color=tuple(self.uniform_color) if self.use_uniform_color else None,
            selection_attribute=attribute,
            name=name,
        )
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)


@annotated_props
class MN_OT_Remove_Style(Operator):
    """
    Operator to remove a style from an entity
    """

    bl_idname = "mn.remove_style"
    bl_label = "Remove Style"
    bl_description = (
        "Remove the selected style, along with its selection and color nodes"
    )
    bl_options = {"REGISTER", "UNDO"}

    # the node is identified by tree + node name, as operators can't take
    # pointers; working on Blender data directly (rather than through the
    # session) keeps removal available for unlinked objects
    name_tree: Annotated[str, StringProperty()]
    name_node: Annotated[str, StringProperty()]

    def execute(self, context: Context):
        tree = bpy.data.node_groups.get(self.name_tree)
        if tree is None or self.name_node not in tree.nodes:
            self.report({"ERROR"}, "Style node to remove was not found")
            return {"CANCELLED"}
        remove_style_node(tree.nodes[self.name_node])
        # keep the styles list selection valid by pointing it at the last
        # remaining style node, on every object rendered by this tree
        last = max(
            (i for i, node in enumerate(tree.nodes) if "Style" in node.name),
            default=-1,
        )
        for obj in bpy.data.objects:
            mod = obj.modifiers.get("Molecular Nodes")
            if mod is not None and getattr(mod, "node_group", None) == tree:
                object_props(obj).styles_active_index = last
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_confirm(self, event, title="Remove Style?")


def _swap_style_items(self, context):
    # module-level tuples, so the references stay alive (avoids the dynamic-enum
    # garbage-collection crash)
    obj = context.active_object
    if obj is not None and obj.mn.entity_type == "density":
        return DENSITY_STYLE_ITEMS
    return STYLE_ITEMS


@annotated_props
class MN_OT_Swap_Style(Operator):
    """
    Operator to swap a style node for a different style, keeping its connections.
    """

    bl_idname = "mn.swap_style"
    bl_label = "Swap Style"
    bl_description = "Swap the selected style node for a different style"
    bl_options = {"REGISTER", "UNDO"}

    # the node is identified by tree + node name, as operators can't take pointers
    name_tree: Annotated[str, StringProperty()]
    name_node: Annotated[str, StringProperty()]
    style: Annotated[str, EnumProperty(items=_swap_style_items)]

    def execute(self, context: Context):
        tree = bpy.data.node_groups.get(self.name_tree)
        if tree is None or self.name_node not in tree.nodes:
            self.report({"ERROR"}, "Style node to swap was not found")
            return {"CANCELLED"}
        swap(tree.nodes[self.name_node], styles_mapping[self.style])
        return {"FINISHED"}


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
    @annotated_props
    class TempAnnotationAddOperator(bpy.types.Operator):
        bl_idname = "mn.temp_annotation_add"
        bl_label = "Add annotation"
        bl_description = "Add a new annotation"

        props: Annotated[Any, bpy.props.PointerProperty(type=AnnotationProps)]

        def draw(self, context):
            layout = self.layout
            assert layout is not None
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


@annotated_props
class MN_OT_Add_Annotation(Operator):
    """
    Operator to add a new annotation to an entity
    """

    bl_idname = "mn.add_annotation"
    bl_label = "Add annotation"
    bl_description = "Add annotation to entity"

    uuid: Annotated[str, StringProperty()]

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
        # operators registered at runtime aren't in the bpy stubs
        bpy.ops.mn.temp_annotation_add("INVOKE_DEFAULT")  # ty: ignore[unresolved-attribute]
        return {"FINISHED"}


@annotated_props
class MN_OT_Remove_Annotation(Operator):
    """
    Operator to remove an annotation from an entity
    """

    bl_idname = "mn.remove_annotation"
    bl_label = "Remove annotation"
    bl_description = "Remove annotation from entity"

    uuid: Annotated[str, StringProperty()]
    annotation_uuid: Annotated[str, StringProperty()]

    def execute(self, context: Context):
        entity = get_session().get(self.uuid)
        if not isinstance(entity, (Molecule, density.Density)):
            return {"CANCELLED"}
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
        assert context.scene is not None
        setup_compositor(context.scene)
        return {"FINISHED"}


@annotated_props
class MN_OT_DSSP_init(Operator):
    """
    Operator to add DSSP as the `sec_struct` calculation, or change its selection
    """

    bl_idname = "mn.dssp_init"
    bl_label = "Add DSSP"
    bl_description = (
        "Calculate secondary structure with DSSP on every frame, or update the "
        "selection DSSP is run on if already added"
    )

    uuid: Annotated[str, StringProperty()]
    selection: Annotated[str, StringProperty(default="protein")]

    def execute(self, context: Context):
        entity = get_session().get(self.uuid)
        if not isinstance(entity, Molecule):
            return {"CANCELLED"}
        try:
            entity.calculations.dssp(selection=self.selection)
        except Exception as e:
            self.report({"ERROR"}, f"DSSP failed for selection '{self.selection}': {e}")
            return {"CANCELLED"}
        update_entities(context.scene)
        return {"FINISHED"}


@annotated_props
class MN_OT_HBonds_init(Operator):
    """
    Operator to add hydrogen bonds as the `hbonds` calculation, or change its selection
    """

    bl_idname = "mn.hbonds_init"
    bl_label = "Add Hydrogen Bonds"
    bl_description = (
        "Calculate hydrogen bonds on every frame, for the Style HBonds node, or "
        "update the selection they are calculated for if already added"
    )

    uuid: Annotated[str, StringProperty()]
    selection: Annotated[str, StringProperty(default="all")]

    def execute(self, context: Context):
        entity = get_session().get(self.uuid)
        if not isinstance(entity, Molecule):
            return {"CANCELLED"}
        try:
            entity.calculations.hbonds(selection=self.selection)
        except Exception as e:
            self.report(
                {"ERROR"},
                f"Hydrogen bonds failed for selection '{self.selection}': {e}",
            )
            return {"CANCELLED"}
        update_entities(context.scene)
        return {"FINISHED"}


@annotated_props
class MN_OT_FrameAttribute_Remove(Operator):
    """
    Operator to remove a selection or calculation from the per-frame attributes
    """

    bl_idname = "mn.frame_attribute_remove"
    bl_label = "Remove"
    bl_description = (
        "Remove this selection (and its attribute) or calculation (its attribute "
        "keeps the values it last had)"
    )
    bl_options = {"REGISTER", "UNDO"}

    uuid: Annotated[str, StringProperty()]
    name: Annotated[str, StringProperty()]

    def execute(self, context: Context):
        entity = get_session().get(self.uuid)
        if not isinstance(entity, Molecule):
            return {"CANCELLED"}
        item = entity.props.frame_attributes.get(self.name)
        if item is None:
            return {"CANCELLED"}
        if item.source == "SELECTION":
            entity.selections.remove(self.name)
        elif item.source == "CALCULATION" and self.name in entity.calculations:
            del entity.calculations[self.name]
        else:
            self.report({"ERROR"}, f"'{self.name}' can't be removed")
            return {"CANCELLED"}
        return {"FINISHED"}


@annotated_props
class MN_OT_DSSP_apply(Operator):
    """
    Operator to apply changed DSSP options
    """

    bl_idname = "mn.dssp_apply"
    bl_label = "Apply"
    bl_description = "Apply changed DSSP options"

    uuid: Annotated[str, StringProperty()]
    apply_ta_threshold: Annotated[bool, BoolProperty()]
    ta_threshold: Annotated[float, FloatProperty()]

    def execute(self, context: Context):
        entity = get_session().get(self.uuid)
        if not isinstance(entity, Molecule):
            return {"CANCELLED"}
        props = entity.props.dssp
        if props.display_option == "trajectory-average":
            if self.apply_ta_threshold:
                entity.calculations._dssp.show_trajectory_average(
                    threshold=self.ta_threshold
                )
            else:
                entity.calculations._dssp.show_trajectory_average()
        return {"FINISHED"}


@annotated_props
class MN_OT_DSSP_cancel(Operator):
    """
    Operator to cancel and restore current DSSP options
    """

    bl_idname = "mn.dssp_cancel"
    bl_label = "Cancel"
    bl_description = "Restore current DSSP options"

    uuid: Annotated[str, StringProperty()]

    def execute(self, context: Context):
        entity = get_session().get(self.uuid)
        if not isinstance(entity, Molecule):
            return {"CANCELLED"}
        props = entity.props.dssp
        props.cancelling = True
        props.display_option = entity.calculations._dssp._display_option
        props.window_size = entity.calculations._dssp._window_size
        if entity.calculations._dssp._sw_threshold is not None:
            props.sw_threshold = entity.calculations._dssp._sw_threshold
            props.apply_sw_threshold = True
        else:
            props.apply_sw_threshold = False
        if entity.calculations._dssp._ta_threshold is not None:
            props.ta_threshold = entity.calculations._dssp._ta_threshold
            props.apply_ta_threshold = True
        else:
            props.apply_ta_threshold = False
        props.applied = True
        props.cancelling = False
        return {"FINISHED"}


@annotated_props
class MN_OT_Align(Operator):
    bl_idname = "mn.align"
    bl_label = "Align Molecules"
    bl_description = (
        "Superimpose the selected molecules onto the active one by setting their"
        " transforms. Residues are paired by sequence or structure, so the"
        " sequences don't need to match"
    )
    bl_options = {"REGISTER", "UNDO"}

    method: Annotated[
        str,
        EnumProperty(
            name="Method",
            items=(
                (
                    "sequence",
                    "Sequence",
                    "Pair residues by sequence alignment and fit on their CA / P atoms,"
                    " dropping outliers",
                ),
                (
                    "structure",
                    "Structure",
                    "Pair residues by structural similarity (3Di), for remote homologs"
                    " with low sequence identity. Proteins only",
                ),
            ),
            default="sequence",
        ),
    ]
    selection: Annotated[
        str,
        StringProperty(
            name="Selection",
            description="MDAnalysis selection of the atoms to fit on in the selected"
            " molecules, e.g. 'chainID H L'. Empty uses all protein and nucleic acid"
            " chains",
            default="",
        ),
    ]
    fixed_selection: Annotated[
        str,
        StringProperty(
            name="Selection",
            description="MDAnalysis selection of the atoms to fit on in the active"
            " molecule, e.g. 'chainID A'. Empty uses all protein and nucleic acid"
            " chains",
            default="",
        ),
    ]

    @staticmethod
    def _molecules(context):
        session = context.scene.MNSession
        fixed = session.match(context.active_object) if context.active_object else None
        mobiles = [
            session.match(obj)
            for obj in context.selected_objects
            if obj != context.active_object
        ]
        return fixed, [m for m in mobiles if isinstance(m, Molecule)]

    @classmethod
    def poll(cls, context):
        fixed, mobiles = cls._molecules(context)
        if not isinstance(fixed, Molecule):
            cls.poll_message_set("The active object must be a molecule")
            return False
        if not mobiles:
            cls.poll_message_set("Select molecules to align, then the reference last")
            return False
        return True

    def invoke(self, context, event):
        # open a dialog first: the chains often need picking before the fit
        return context.window_manager.invoke_props_dialog(
            self, width=450, confirm_text="Align"
        )

    @staticmethod
    def _chain_summary(mol: Molecule) -> str:
        """The chains of `mol` with their residue counts, e.g. 'A 102, H 117'."""
        polymer = mol.universe.select_atoms("protein or nucleic")
        chains = getattr(polymer, "chainIDs", polymer.segids)
        counts = {
            chain: len(polymer[chains == chain].residues)
            for chain in dict.fromkeys(chains)
        }
        return ", ".join(f"{chain} {n}" for chain, n in counts.items()) or "none"

    def _draw_side(self, layout, label: str, mols: list, prop: str) -> None:
        box = layout.box()
        box.label(text=label)
        for mol in mols:
            col = box.column(align=True)
            col.label(text=mol.name, icon="OBJECT_DATA")
            col.label(text=f"Chains (residues): {self._chain_summary(mol)}")
        box.prop(self, prop)

    def draw(self, context):
        layout = self.layout
        assert layout is not None
        fixed, mobiles = self._molecules(context)
        layout.prop(self, "method", expand=True)
        self._draw_side(layout, "Reference (active)", [fixed], "fixed_selection")
        self._draw_side(layout, "Move", mobiles, "selection")
        layout.label(
            text="Chains are paired in order; if the counts differ, each side is"
            " joined into one chain",
            icon="INFO",
        )

    def execute(self, context):
        fixed, mobiles = self._molecules(context)
        results = []
        for mobile in mobiles:
            try:
                result = mobile.align_to(
                    fixed,
                    method=self.method,
                    selection=self.selection or None,
                    fixed_selection=self.fixed_selection or None,
                )
            except (ValueError, IndexError) as e:
                self.report({"ERROR"}, f"Failed aligning {mobile.name}: {e}")
                return {"CANCELLED"}
            results.append(
                f"{mobile.name}: RMSD {result.rmsd:.2f} Å over {result.n_anchors} atoms"
                + (" (chains joined)" if result.merged_chains else "")
            )
        self.report({"INFO"}, f"Aligned to {fixed.name}. " + "; ".join(results))
        return {"FINISHED"}


CLASSES = [
    MN_OT_add_selection_to_style,
    MN_OT_Import_Molecule,
    MN_OT_Import_OxDNA,
    MN_OT_Frames_To_Collection,
    MN_OT_Import_Map,
    MN_OT_Import_Ensemble,
    MN_FH_Import_Molecule,
    MN_OT_Add_Style,
    MN_OT_Remove_Style,
    MN_OT_Swap_Style,
    MN_OT_Add_Annotation,
    MN_OT_Remove_Annotation,
    MN_OT_Setup_Compositor,
    MN_OT_DSSP_init,
    MN_OT_HBonds_init,
    MN_OT_FrameAttribute_Remove,
    MN_OT_DSSP_apply,
    MN_OT_DSSP_cancel,
    MN_OT_Align,
]
