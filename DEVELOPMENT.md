

## Environment

- [uv](https://docs.astral.sh/uv/)
- Blender
- [blext](https://pypi.org/project/blext/)


## Starting


```python
# uvx blext show spec
BLExtSpec(
    path_proj_root=PosixPath('/Users/zcpowers/Documents/Projects/MolecularNodes'),
    req_python_version='~=3.11',
    is_universal_blext=True,
    bl_platform_pypa_tags={},
    init_settings_filename='init_settings.toml',
    manifest_filename='blender_manifest.toml',
    init_schema_version='0.1.0',
    use_path_local=False,
    use_log_file=False,
    log_file_path=PosixPath('addon.log'),
    log_file_level=<StrLogLevel.Debug: 'DEBUG'>,
    use_log_console=True,
    log_console_level=<StrLogLevel.Debug: 'DEBUG'>,
    manifest_schema_version='1.0.0',
    id='molecularnodes',
    name='BLExt Simple Example',
    version='4.4.0',
    tagline='Toolbox for molecular animations with Blender and Geometry Nodes.',
    maintainer='Brady Johnston <brady.johnston@me.com>',
    type='add-on',
    blender_version_min='4.2.0',
    blender_version_max='4.3.10',
    permissions={},
    tags=('Development',),
    license=('SPDX:GLP',),
    copyright=('2024 blext Contributors',),
    website=HttpUrl('https://bradyajohnston.github.io/MolecularNodes'),
    bl_platforms=frozenset(),
    wheels=()
)
```



```sh
# you need to use the head
# see https://codeberg.org/so-rose/blext/issues/57
uvx --from git+https://codeberg.org/so-rose/blext blext build
uvx --from git+https://codeberg.org/so-rose/blext blext check # not implemented
uvx --from git+https://codeberg.org/so-rose/blext blext run
uvx --from git+https://codeberg.org/so-rose/blext blext show blender_manifest # broken
uvx --from git+https://codeberg.org/so-rose/blext blext show deps
uvx --from git+https://codeberg.org/so-rose/blext blext show global_config
uvx --from git+https://codeberg.org/so-rose/blext blext show profile
uvx --from git+https://codeberg.org/so-rose/blext blext show path

uv add  git+https://codeberg.org/so-rose/blext


# show
uvx blext show manifest
uvx blext show init_settings
uvx blext show spec   # see above
uvx blext show path_blender # works locally

# build
uvx blext build

# dev
uvx blext dev
```
