# MolViewSpec example files

Official example states vendored for the import tests:

- `landing_*.mvsj` — the landing-page examples at
  <https://molstar.org/mol-view-spec/>, from
  `landing/public/examples/*/state.mvsj` in the
  [molstar/mol-view-spec](https://github.com/molstar/mol-view-spec) repository.
- `colab_*.mvsj` — from `test-data/colab_examples/` in the same repository.

The files reference structures by wwPDB URL, so the tests that import them
download on first run (cached under `MolecularNodesCache/mvs/`).
