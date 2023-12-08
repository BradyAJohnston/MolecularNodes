import molecularnodes as mn

def test_name_versions():
    for name in mn.pkg.get_pkgs().keys():
        print(f"{name}")
        print(f"{mn.pkg.is_current(name)}")
        assert mn.pkg.is_current(name)

def test_is_current():
    assert mn.pkg.is_current('biotite')

def test_get_pkgs():
    names = ['biotite', 'MDAnalysis', 'mrcfile', 'starfile', 'msgpack']
    for name in mn.pkg.get_pkgs().keys():
        assert name in names