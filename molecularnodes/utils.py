import os
import traceback
import zipfile
from pathlib import Path

import bpy
import numpy as np
from bpy.app.translations import pgettext_tip as tip_
from mathutils import Matrix

ADDON_DIR = Path(__file__).resolve().parent


def lerp(a: np.ndarray, b: np.ndarray, t: float = 0.5) -> np.ndarray:
    """
    Linearly interpolate between two values.

    Parameters
    ----------
    a : array_like
        The starting value.
    b : array_like
        The ending value.
    t : float, optional
        The interpolation parameter. Default is 0.5.

    Returns
    -------
    array_like
        The interpolated value(s).

    Notes
    -----
    This function performs linear interpolation between `a` and `b` using the
    interpolation parameter `t` such that the result lies between `a` and `b`.

    Examples
    --------
    >>> lerp(1, 2, 0.5)
    1.5

    >>> lerp(3, 7, 0.2)
    3.8

    >>> lerp([1, 2, 3], [4, 5, 6], 0.5)
    array([2.5, 3.5, 4.5])

    """
    return np.add(a, np.multiply(np.subtract(b, a), t))


def _module_filesystem_remove(path_base, module_name):
    # taken from the bpy.ops.preferences.app_template_install() operator source code
    # Remove all Python modules with `module_name` in `base_path`.
    # The `module_name` is expected to be a result from `_zipfile_root_namelist`.
    import os
    import shutil

    module_name = os.path.splitext(module_name)[0]
    for f in os.listdir(path_base):
        f_base = os.path.splitext(f)[0]
        if f_base == module_name:
            f_full = os.path.join(path_base, f)
            if os.path.isdir(f_full):
                shutil.rmtree(f_full)
            else:
                os.remove(f_full)


def _zipfile_root_namelist(file_to_extract):
    # taken from the bpy.ops.preferences.app_template_install() operator source code
    # Return a list of root paths from zipfile.ZipFile.namelist.
    import os

    root_paths = []
    for f in file_to_extract.namelist():
        # Python's `zipfile` API always adds a separate at the end of directories.
        # use `os.path.normpath` instead of `f.removesuffix(os.sep)`
        # since paths could be stored as `./paths/./`.
        #
        # Note that `..` prefixed paths can exist in ZIP files but they don't write to parent directory when extracting.
        # Nor do they pass the `os.sep not in f` test, this is important,
        # otherwise `shutil.rmtree` below could made to remove directories outside the installation directory.
        f = os.path.normpath(f)
        if os.sep not in f:
            root_paths.append(f)
    return root_paths


def template_install():
    print(os.path.abspath(ADDON_DIR))
    template = os.path.join(
        os.path.abspath(ADDON_DIR), "assets", "template", "Molecular Nodes.zip"
    )
    _install_template(template)
    bpy.utils.refresh_script_paths()


def template_uninstall():
    import shutil

    for folder in bpy.utils.app_template_paths():
        path = os.path.join(os.path.abspath(folder), "MolecularNodes")
        if os.path.exists(path):
            shutil.rmtree(path)
    bpy.utils.refresh_script_paths()


def _install_template(filepath, subfolder="", overwrite=True):
    # taken from the bpy.ops.preferences.app_template_install() operator source code

    path_app_templates = bpy.utils.user_resource(
        "SCRIPTS",
        path=os.path.join("startup", "bl_app_templates_user", subfolder),
        create=True,
    )

    if not os.path.isdir(path_app_templates):
        try:
            os.makedirs(path_app_templates, exist_ok=True)
        except:
            traceback.print_exc()

    app_templates_old = set(os.listdir(path_app_templates))

    # check to see if the file is in compressed format (.zip)
    if zipfile.is_zipfile(filepath):
        try:
            file_to_extract = zipfile.ZipFile(filepath, "r")
        except:
            traceback.print_exc()
            return {"CANCELLED"}

        file_to_extract_root = _zipfile_root_namelist(file_to_extract)
        if overwrite:
            for f in file_to_extract_root:
                _module_filesystem_remove(path_app_templates, f)
        else:
            for f in file_to_extract_root:
                path_dest = os.path.join(path_app_templates, os.path.basename(f))
                if os.path.exists(path_dest):
                    # self.report({'WARNING'}, tip_("File already installed to %r\n") % path_dest)
                    return {"CANCELLED"}

        try:  # extract the file to "bl_app_templates_user"
            file_to_extract.extractall(path_app_templates)
        except:
            traceback.print_exc()
            return {"CANCELLED"}

    else:
        # Only support installing zipfiles
        print("no zipfile")
        return {"CANCELLED"}

    app_templates_new = set(os.listdir(path_app_templates)) - app_templates_old

    # in case a new module path was created to install this addon.
    bpy.utils.refresh_script_paths()

    # print message
    msg = tip_("Template Installed (%s) from %r into %r") % (
        ", ".join(sorted(app_templates_new)),
        filepath,
        path_app_templates,
    )
    print(msg)


# data types for the np.array that will store per-chain symmetry operations
dtype = [
    ("assembly_id", int),
    ("transform_id", int),
    ("chain_id", "U10"),
    ("rotation", float, 4),  # quaternion form
    ("translation", float, 3),
]


def array_quaternions_from_dict(transforms_dict):
    n_transforms = 0
    for assembly in transforms_dict.values():
        for transform in assembly:
            n_transforms += len(transform[0])

    arr = np.array((n_transforms), dtype=dtype)

    transforms = []
    for i, assembly in enumerate(transforms_dict.values()):
        for j, transform in enumerate(assembly):
            chains = transform[0]
            matrix = transform[1]
            arr = np.zeros((len(chains)), dtype=dtype)
            translation, rotation, scale = Matrix(matrix).decompose()
            arr["assembly_id"] = i + 1
            arr["transform_id"] = j
            arr["chain_id"] = chains
            arr["rotation"] = rotation
            arr["translation"] = translation
            transforms.append(arr)

    return np.hstack(transforms)
