docs-build:
	cd docs && quartodoc build
	quarto render docs

test:
	pip install .
	pytest -vv

version := $(shell grep -o -E "\b[0-9]+\.[0-9]+\.[0-9]+\b" pyproject.toml)

release:
	git clean -dfX
	zip -r MolecularNodes_$(version).zip MolecularNodes -x *pycache* *.blend1
