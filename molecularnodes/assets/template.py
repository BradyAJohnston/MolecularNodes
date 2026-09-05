import itertools
import os
import shutil
from pathlib import Path
import bpy

SUBFOLDER = "Molecular Nodes"


def list_templates():
    paths = bpy.utils.app_template_paths()
    t_paths = list(paths) if paths is not None else []
    names = list(
        itertools.chain.from_iterable(
            [[folder.stem for folder in Path(p).rglob("")] for p in t_paths]
        )
    )
    return [name for name in names if not name.startswith("bl_app")]


def is_installed():
    base_path = Path(
        bpy.utils.user_resource(
            "SCRIPTS", path=str(Path("startup") / "bl_app_templates_user"), create=False
        )
    )
    molecular_nodes_path = base_path / "Molecular Nodes"
    return molecular_nodes_path.exists()


def install():
    base_path = Path(
        bpy.utils.user_resource(
            "SCRIPTS", path=str(Path("startup") / "bl_app_templates_user"), create=True
        )
    )

    path_app_templates = base_path / SUBFOLDER
    startup_file = Path(__file__).parent / "template/startup.blend"
    destination = path_app_templates / startup_file.name

    # copy to a process-unique name and atomically replace, so concurrent
    # installs (e.g. parallel test workers) never collide writing the same
    # file, which fails on Windows; losing the race is fine as long as one
    # of the installs got the template in place. Retry once in case a
    # concurrent uninstall removed the folder mid-copy.
    temporary = path_app_templates / f"startup-{os.getpid()}.blend.tmp"
    last_error: OSError | None = None
    for _ in range(2):
        path_app_templates.mkdir(parents=True, exist_ok=True)
        try:
            shutil.copy(startup_file, temporary)
            os.replace(temporary, destination)
            last_error = None
            break
        except OSError as e:
            last_error = e
        finally:
            temporary.unlink(missing_ok=True)
    if last_error is not None and not destination.exists():
        raise last_error
    bpy.utils.refresh_script_paths()


def uninstall():
    for folder in bpy.utils.app_template_paths():
        path = Path(folder).absolute() / SUBFOLDER
        if "Molecular Nodes" not in str(path):
            continue

        if path.exists():
            try:
                shutil.rmtree(path)
            except OSError:
                # a concurrent install may be re-populating the folder
                # (e.g. parallel test workers); retry, then give up quietly
                shutil.rmtree(path, ignore_errors=True)
    bpy.utils.refresh_script_paths()
