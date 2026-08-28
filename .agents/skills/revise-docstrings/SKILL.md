---
name: revise-docstrings
description: >
  Review and improve Python docstrings for Great Docs API reference
  generation. Covers NumPy and Google style conventions, parameter
  documentation, return types, examples, cross-references, and
  Great Docs directives (%seealso, %nodoc). Use when auditing,
  writing, or fixing docstrings in a Python package.
license: MIT
compatibility: Requires Great Docs >=0.8, Quarto CLI installed.
metadata:
  author: rich-iannone
  version: "1.0"
  tags:
    - documentation
    - docstrings
    - python
    - api-reference
---

# Revise Docstrings

Skill for reviewing and improving Python docstrings so they render
correctly and completely in a Great Docs API reference site.

## Quick start

```bash
# Preview what Great Docs will document
great-docs scan --verbose

# Build and check the rendered reference
great-docs build && great-docs preview
```

## Skill directory structure

```
skills/revise-docstrings/
├── SKILL.md
└── references/
    ├── docstring-checklist.md
    └── style-examples.md
```

## When to use this skill

| Need                               | Action                                  |
| ---------------------------------- | --------------------------------------- |
| Audit docstring completeness       | Run checklist on each public symbol     |
| Convert Google→NumPy (or reverse)  | Reformat following style-examples.md    |
| Add missing parameter descriptions | Fill in Parameters section              |
| Add a live example                 | Use Examples section with `>>>` prompts |
| Cross-reference related symbols    | Add `%seealso` directive                |
| Hide an internal symbol            | Add `%nodoc` directive                  |
| Fix rendering issues               | Check against common pitfalls below     |

## Core concepts

### Docstring formats

Great Docs supports two formats. Set the parser in config:

```yaml
# great-docs.yml
parser: numpy # or "google"
```

Both produce identical rendered output. Use whichever your codebase
already uses.

### NumPy style

````python
def connect(host: str, port: int = 5432) -> Connection:
    """
    Open a connection to the database server.

    Establishes a TCP connection to the specified host and port, authenticates
    with default credentials, and returns an active connection handle.

    Parameters
    ----------
    host
        Hostname or IP address of the database server.
    port
        TCP port number. Defaults to `5432`.

    Returns
    -------
    Connection
        An authenticated connection ready for queries.

    Raises
    ------
    ConnectionError
        If the server is unreachable.

    Examples
    --------
    Connect to a local server and run a simple query:

    ```{python}
    conn = connect("localhost")
    conn.execute("SELECT 1")
    ```

    See Also
    --------
    disconnect : Close a connection.

    Notes
    -----
    The connection uses TLS by default when port 5433 is specified.
    """
````

### Google style

````python
def connect(host: str, port: int = 5432) -> Connection:
    """Open a connection to the database server.

    Establishes a TCP connection to the specified host and port, authenticates
    with default credentials, and returns an active connection handle.

    Args:
        host: Hostname or IP address of the database server.
        port: TCP port number. Defaults to `5432`.

    Returns:
        An authenticated connection ready for queries.

    Raises:
        ConnectionError: If the server is unreachable.

    Examples:
        Connect and run a query:

        ```{python}
        conn = connect("localhost")
        conn.execute("SELECT 1")
        ```
    """
````

### Sections recognized by Great Docs

| Section          | NumPy header | Google header | Purpose                                |
| ---------------- | ------------ | ------------- | -------------------------------------- |
| Summary          | first line   | first line    | One-line description                   |
| Extended summary | body text    | body text     | Multi-paragraph explanation            |
| Parameters       | `Parameters` | `Args:`       | Function/method arguments              |
| Returns          | `Returns`    | `Returns:`    | Return value description               |
| Raises           | `Raises`     | `Raises:`     | Exceptions that may be raised          |
| Examples         | `Examples`   | `Examples:`   | Usage examples (Quarto cells or `>>>`) |
| Notes            | `Notes`      | `Notes:`      | Implementation details                 |
| See Also         | `See Also`   | —             | Related symbols                        |
| Warns            | `Warns`      | —             | Warnings issued                        |
| References       | `References` | `References:` | Citations or links                     |

### Great Docs directives

Special inline directives using `%` prefix:

```python
def my_function():
    """
    Description.

    %seealso other_func, SomeClass: related utilities
    %nodoc
    """
```

| Directive  | Effect                                              |
| ---------- | --------------------------------------------------- |
| `%seealso` | Renders a "See Also" box with cross-reference links |
| `%nodoc`   | Excludes this symbol from the API reference         |

### Type annotations vs docstring types

- **Prefer type annotations** in the function signature.
- Great Docs reads annotations automatically — no need to duplicate
  types in the docstring. Write bare parameter names (e.g., `host`
  not `host : str`) and let the signature annotation render on the
  reference page.
- If you _do_ add a type in the docstring (e.g., `host : str`), it
  **overwrites** the annotation in the rendered output. This can be
  useful to show a simplified form (e.g., `str or Path` instead of
  `str | pathlib.Path`), but it creates a maintenance risk: the
  docstring type and the annotation can drift apart silently.
- **Rule of thumb:** omit docstring types unless the annotation is
  confusing to readers. Keep one source of truth.

## Workflows

### Auditing a package's docstrings

```
Task Progress:
- [ ] Step 1: List public API
- [ ] Step 2: Check each symbol
- [ ] Step 3: Fix issues
- [ ] Step 4: Rebuild and verify
```

**Step 1**: Run `great-docs scan --verbose` to see every symbol
Great Docs will document.

**Step 2**: For each symbol, run through the checklist in
[references/docstring-checklist.md](references/docstring-checklist.md).

**Step 3**: Fix missing or incorrect sections. Use
[references/style-examples.md](references/style-examples.md) as
a template.

**Step 4**: Rebuild with `great-docs build` and check the rendered
pages in the browser.

### Writing a docstring from scratch

1. Start with a one-line summary (imperative mood: "Connect to...",
   "Return the...", "Parse the...").
2. Add an extended summary if the one-liner is insufficient.
3. Document every parameter with name, type, and description.
4. Document the return value.
5. Add `Raises` if the function can raise exceptions.
6. Add an `Examples` section with Quarto code cells (preferred)
   or `>>>` prompts. Use `{python}` for executable cells,
   `{.python}` for display-only. Wrap cells with short prose.
7. Add `%seealso` for closely related symbols.

### Fixing a rendering issue

Common rendering problems and their fixes:

| Problem                         | Cause                          | Fix                                         |
| ------------------------------- | ------------------------------ | ------------------------------------------- |
| Parameter not showing           | Wrong indentation              | Align to 4 spaces under section header      |
| Code block not rendering        | Missing blank line before code | Add blank line above cell or `>>>`          |
| Type mismatch warning           | Annotation ≠ docstring type    | Remove type from docstring, keep annotation |
| Symbol missing from reference   | Not exported in `__init__.py`  | Add import to `__init__.py`                 |
| Entire docstring shown as prose | Wrong parser setting           | Check `parser:` in great-docs.yml           |

## Gotchas

1. **One-line summary is required.** Without it, the API reference
   page shows no description at all.
2. **Blank line after summary.** NumPy style requires a blank line
   between the summary and the extended summary.
3. **Indentation matters.** Parameters must be indented consistently
   (4 spaces for NumPy, nested under `Args:` for Google).
4. **Don't mix styles.** All docstrings in a package must use the
   same format. Mixing causes parsing failures.
5. **`%nodoc` hides completely.** Use `exclude` in config for
   selective hiding without modifying source code.
6. **Prefer Quarto cells over `>>>` prompts.** `{python}` cells
   render output automatically and support prose between steps.
   Use `{.python}` for non-executable illustration.
7. **Class docstrings go on the class, not `__init__`.** Great Docs
   reads the class-level docstring for the class page.
