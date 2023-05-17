docs-build:
	cd docs && quartodoc build
	quarto render docs

test:
	pip3 install .
	pytest -vv
