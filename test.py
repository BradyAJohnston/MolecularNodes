import molecularnodes as mn

m = mn.io.parse.bcif.BBCIF(mn.io.download('4ozs', 'bcif', cache='/Users/brady/.MolecularNodes'))

m.create_model()
