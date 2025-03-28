#
# Core commands for Developing with Blender
#
# try `just -l` in your shell
#


# Builds the MN Extension
build:
    echo 'Building Python Packages for Molecular Nodes Extension'
    uv run build.py


# Creates the Website. NOte: requires `quarto`
document: update-deps
    echo 'Creating the UV docs'
    cd docs && \
    uv run generate.py && \
    uv run quartodoc build ** \
    cd .. && \
    uv run quarto render docs

# Installs MN as a local Blender Extension
install-in-blender:
    echo 'Creating the UV docs'
    blender -b -P tests/python.py -- -m pip install ".[test]"
    blender -b -P tests/run.py -- -vv tests --cov --cov-report=xml

# run Blender with with a Jupyter connecttion
run-jupyter:
    uv run jupyter lab # needs tweaking to work


update-deps:
    echo 'Updating Molecular Nodes Dependencies'
    uv sync --all-extras

# Runs the testsuite
test:
  echo 'This is another recipe.'
  uv run pytest
