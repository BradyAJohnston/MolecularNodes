docs-build:
	# cd docs && quartodoc build
	quarto render docs

test:
	pip install -e .
	pytest -v

version := $(shell grep version pyproject.toml | grep -o -E "\b[0-9]+\.[0-9]+\.[0-9]+\b")

template:
	cd 'molecularnodes/assets/template/' && zip -r 'Molecular Nodes.zip' 'Molecular Nodes'

release:
	git clean -dfX
	make template
	zip -r molecularnodes_$(version).zip molecularnodes -x *pycache* *.blend1 "molecularnodes/assets/template/MolecularNodes/*"
