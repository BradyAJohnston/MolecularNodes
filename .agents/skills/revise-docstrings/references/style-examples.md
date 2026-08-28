# Style Examples — Docstring Templates

## Function — NumPy style

```python
def read_csv(
    path: str | Path,
    delimiter: str = ",",
    header: bool = True,
) -> DataFrame:
    """
    Read a CSV file into a DataFrame.

    Parses the file at *path* using the specified delimiter and returns a
    DataFrame. The first row is treated as column headers by default.

    Parameters
    ----------
    path
        Path to the CSV file. Supports local paths and `s3://` URIs.
    delimiter
        Column separator character. Defaults to `","`.
    header
        If `True`, the first row is used as column names. If `False`, columns
        are numbered `0, 1, 2, ...`.

    Returns
    -------
    DataFrame
        Parsed tabular data with one row per record.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If the file contains no data rows.

    Examples
    --------
    Read a simple CSV file and inspect its shape:

    ```{python}
    df = read_csv("data/sales.csv")
    df.shape
    ```

    Use a tab delimiter for TSV files:

    ```{.python}
    df = read_csv("data/tsv_file.tsv", delimiter="\t")
    ```

    See Also
    --------
    write_csv : Write a DataFrame to CSV.
    read_parquet : Read Parquet files.

    Notes
    -----
    Large files (>1 GB) are read in chunks internally. Memory usage stays
    bounded regardless of file size.
    """
```

### Parameter types — when to include them

In NumPy-style docstrings the parameter name can optionally carry a
type after a colon (`path : str or Path`). When the function
signature already has a type annotation, **omit the type from the
docstring** — Great Docs reads the annotation automatically and
renders it on the reference page.

If you *do* add a type in the docstring, it **overwrites** the
annotation in the rendered output. This is occasionally useful (e.g.,
to show a simplified union like `str or Path` instead of
`str | pathlib.Path`), but it creates a maintenance burden: the
docstring type and the annotation can drift out of sync silently.

**Rule of thumb:** leave parameter lines as bare names and let the
annotations speak for themselves. Add a docstring type only when the
rendered annotation would confuse readers.

## Examples — Quarto code cells

Great Docs renders docstrings through Quarto, so the Examples section
can use Quarto code cells instead of `>>>` prompts. This gives you
executable output, syntax highlighting, and the ability to weave
prose between cells.

### Executable cell

Use `{python}` (curly braces, no dot) to make a cell that runs
during the build. The output is captured and rendered inline:

```
    Examples
    --------
    Read a CSV and check the row count:

    ```{python}
    df = read_csv("data/sales.csv")
    df.shape
    ```
```

Quarto executes the cell, so readers see both the code and its
result — no need to hardcode output.

### Non-executable cell

Use `{.python}` (dot prefix) when the code is illustrative but
should not run (e.g., it depends on external resources or is
expensive):

```
    Examples
    --------
    Read a remote Parquet file:

    ```{.python}
    df = read_parquet("s3://bucket/warehouse/events.parquet")
    ```
```

The cell renders with syntax highlighting but no execution.

### Mixing cells with prose

Wrap each cell in a short sentence or two of context. This produces
a mini-tutorial on the reference page:

```
    Examples
    --------
    Open a connection to the default local server:

    ```{python}
    conn = connect("localhost")
    ```

    Run a query and fetch results as a list of rows:

    ```{python}
    rows = conn.execute("SELECT id, name FROM users LIMIT 5")
    rows
    ```

    Always close the connection when finished:

    ```{.python}
    conn.close()
    ```
```

### When to use `>>>` prompts

The traditional `>>>` style still works and is fine for simple,
self-contained snippets. Prefer Quarto cells when:

- The example benefits from rendered output (tables, plots)
- You want prose between steps
- The example has multiple logical steps

## Function — Google style

```python
def read_csv(
    path: str | Path,
    delimiter: str = ",",
    header: bool = True,
) -> DataFrame:
    """Read a CSV file into a DataFrame.

    Parses the file at *path* using the specified delimiter and returns a
    DataFrame. The first row is treated as column headers by default.

    Args:
        path: Path to the CSV file. Supports local paths and `s3://` URIs.
        delimiter: Column separator character. Defaults to `","`.
        header: If `True`, the first row is used as column names. If `False`,
            columns are numbered `0, 1, 2, ...`.

    Returns:
        Parsed tabular data with one row per record.

    Raises:
        FileNotFoundError: If *path* does not exist.
        ValueError: If the file contains no data rows.

    Examples:
        Read a CSV and check its shape:

        ```{python}
        df = read_csv("data/sales.csv")
        df.shape
        ```
    """
```

## Class — NumPy style

```python
class Connection:
    """
    A database connection handle.

    Wraps a TCP socket to the database server and provides methods for
    executing queries, managing transactions, and streaming result sets.

    Parameters
    ----------
    host
        Server hostname or IP address.
    port
        TCP port. Defaults to `5432`.
    timeout
        Connection timeout in seconds. Defaults to `30.0`.

    Examples
    --------
    Create a connection and run a query:

    ```{python}
    conn = Connection("localhost")
    result = conn.execute("SELECT 1")
    conn.close()
    ```

    See Also
    --------
    connect : Factory function for creating connections.
    ConnectionPool : Manage a pool of reusable connections.
    """

    def __init__(self, host: str, port: int = 5432, timeout: float = 30.0):
        ...
```

## Property — NumPy style

```python
@property
def is_connected(self) -> bool:
    """
    Check whether the connection is still alive.

    Returns
    -------
    bool
        `True` if the connection is open and responsive.
    """
```

## Directive usage

```python
def internal_helper():
    """
    Normalize column names.

    %nodoc
    """

def public_function():
    """
    Transform the dataset.

    %seealso other_transform, Pipeline: related transformation tools
    """
```
