docs-build:
	cd docs && quartodoc build
	quarto render docs

test:
	pip install .
	pytest -v

version := $(shell grep version pyproject.toml | grep -o -E "\b[0-9]+\.[0-9]+\.[0-9]+\b")

template:
	cd MolecularNodes/assets/template && zip -r MolecularNodes.zip MolecularNodes -x *blend1

release:
	git clean -dfX
	make template
	zip -r MolecularNodes_$(version).zip MolecularNodes -x *pycache* *.blend1 MolecularNodes/assets/template/Molecular_Nodes/*
