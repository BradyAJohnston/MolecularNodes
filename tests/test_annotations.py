import bpy
import MDAnalysis as mda
import pytest
import molecularnodes as mn
from .constants import data_dir

bpy.utils.expose_bundled_modules()


class TestAnnotations:
    @pytest.fixture(scope="module")
    def universe(self):
        topo = data_dir / "md_ppr/box.gro"
        traj = data_dir / "md_ppr/first_5_frames.xtc"
        return mda.Universe(topo, traj)

    @pytest.fixture(scope="module")
    def session(self):
        return mn.session.get_session()

    @pytest.fixture(scope="module")
    def density_file(self):
        file = data_dir / "emd_24805.map.gz"
        vdb_file = data_dir / "emd_24805.vdb"
        vdb_file.unlink(missing_ok=True)
        # Make all densities are removed
        for o in bpy.data.objects:
            if o.mn.entity_type == "density":
                bpy.data.objects.remove(o, do_unlink=True)
        return file

    def test_registered_trajectory_annotations(self):
        # test if trajectory annotations are correctly auto registered
        manager = mn.entities.trajectory.TrajectoryAnnotationManager
        # atom_info
        assert "atom_info" in manager._classes
        assert hasattr(manager, "add_atom_info")
        assert callable(getattr(manager, "add_atom_info"))
        # com
        assert "com" in manager._classes
        assert hasattr(manager, "add_com")
        assert callable(getattr(manager, "add_com"))
        # com_distance
        assert "com_distance" in manager._classes
        assert hasattr(manager, "add_com_distance")
        assert callable(getattr(manager, "add_com_distance"))
        # canonical_dihedrals
        assert "canonical_dihedrals" in manager._classes
        assert hasattr(manager, "add_canonical_dihedrals")
        assert callable(getattr(manager, "add_canonical_dihedrals"))
        # universe_info
        assert "universe_info" in manager._classes
        assert hasattr(manager, "add_universe_info")
        assert callable(getattr(manager, "add_universe_info"))
        # label_2d
        assert "label_2d" in manager._classes
        assert hasattr(manager, "add_label_2d")
        assert callable(getattr(manager, "add_label_2d"))
        # label_3d
        assert "label_3d" in manager._classes
        assert hasattr(manager, "add_label_3d")
        assert callable(getattr(manager, "add_label_3d"))

    def test_registered_molecule_annotations(self):
        # test if molecule annotations are correctly auto registered
        manager = mn.entities.molecule.annotations.MoleculeAnnotationManager
        # molecule_info
        assert "molecule_info" in manager._classes
        assert hasattr(manager, "add_molecule_info")
        assert callable(getattr(manager, "add_molecule_info"))
        # label_2d
        assert "label_2d" in manager._classes
        assert hasattr(manager, "add_label_2d")
        assert callable(getattr(manager, "add_label_2d"))
        # label_3d
        assert "label_3d" in manager._classes
        assert hasattr(manager, "add_label_3d")
        assert callable(getattr(manager, "add_label_3d"))

    def test_registered_density_annotations(self):
        # test if density annotations are correctly auto registered
        manager = mn.entities.density.annotations.DensityAnnotationManager
        # density_info
        assert "density_info" in manager._classes
        assert hasattr(manager, "add_density_info")
        assert callable(getattr(manager, "add_density_info"))
        # grid_axes
        assert "grid_axes" in manager._classes
        assert hasattr(manager, "add_grid_axes")
        assert callable(getattr(manager, "add_grid_axes"))
        # label_2d
        assert "label_2d" in manager._classes
        assert hasattr(manager, "add_label_2d")
        assert callable(getattr(manager, "add_label_2d"))
        # label_3d
        assert "label_3d" in manager._classes
        assert hasattr(manager, "add_label_3d")
        assert callable(getattr(manager, "add_label_3d"))

    def test_trajectory_annotations_registration(self, universe, session):
        manager = mn.entities.trajectory.TrajectoryAnnotationManager

        # test register exceptions
        class TestAnnotation:
            pass

        # test presence of annotation_type
        with pytest.raises(ValueError):
            manager.register(TestAnnotation)

        class TestAnnotation:
            annotation_type = "test_annotation"

        # test not an annotation class
        with pytest.raises(ValueError):
            manager.register(TestAnnotation)

        # test no draw method
        with pytest.raises(ValueError):

            class TestAnnotation(mn.entities.trajectory.TrajectoryAnnotation):
                annotation_type = "test_annotation"

        class TestAnnotation:
            annotation_type = "test_annotation"

        # test unregister exceptions
        with pytest.raises(ValueError):
            manager.unregister(TestAnnotation)

        class TestAnnotation(mn.entities.trajectory.TrajectoryAnnotation):
            annotation_type = "test_annotation"

            selection: str  # required param
            param1: int = 0
            param2: float = 0.0

            def draw(self):
                print(
                    self.trajectory,
                    self.interface.selection,
                    self.interface.param1,
                    self.interface.param2,
                )

        # test auto register
        assert "test_annotation" in manager._classes
        assert hasattr(manager, "add_test_annotation")
        # test unregister
        manager.unregister(TestAnnotation)
        assert "test_annotation" not in manager._classes
        assert not hasattr(manager, "add_test_annotation")
        # test register
        manager.register(TestAnnotation)
        assert "test_annotation" in manager._classes
        assert hasattr(manager, "add_test_annotation")
        # test add
        t1 = mn.Trajectory(universe)
        # test add with required param missing
        with pytest.raises(ValueError):
            a1 = t1.annotations.add_test_annotation()
        # test add with invalid param
        with pytest.raises(ValueError):
            a1 = t1.annotations.add_test_annotation(invalid_param="value")
        a1 = t1.annotations.add_test_annotation(selection="all")
        # test draw
        a1._instance.draw()

    def test_trajectory_annotation_lifecycle(self, universe, session):
        t1 = mn.Trajectory(universe)
        # test annotation atom_info add
        assert len(t1.annotations._interfaces) == 0
        a1 = t1.annotations.add_atom_info(selection="all")
        assert len(t1.annotations._interfaces) == 1
        # test len and iterable access
        assert len(t1.annotations) == 1
        for ant in t1.annotations:
            assert ant == a1
        # test subscriptable access
        assert t1.annotations["Annotation"] == a1
        # test index based access
        assert t1.annotations[0] == a1
        # test invalid index
        with pytest.raises(ValueError):
            a1 = t1.annotations[100]
        # test annotation get by name
        ga1 = t1.annotations.get("Annotation")
        assert ga1 == a1
        # test invalid get by name
        with pytest.raises(ValueError):
            ga1 = t1.annotations.get("InvalidAnnotation")
        # test key completions
        list1 = t1.annotations._ipython_key_completions_()
        assert list1 == ["Annotation"]
        # test annotation inputs in interface
        assert hasattr(a1, "selection")
        assert hasattr(a1, "show_resid")
        assert hasattr(a1, "show_segid")
        # test annotation common params in interface
        assert hasattr(a1, "visible")
        # test annotation instance having trajectory entity
        assert a1._instance.trajectory == t1
        # test annotation instance interface
        assert a1 == a1._instance.interface
        # test all annotations visibility
        assert t1.annotations.visible
        t1.annotations.visible = False
        assert not t1.annotations.visible
        # test remove annotation
        a2 = t1.annotations.add_atom_info(selection="resid 2", name="A2")
        a3 = t1.annotations.add_atom_info(selection="resid 3")
        assert a2 in t1.annotations
        assert a3 in t1.annotations
        # remove annotation by name
        t1.annotations.remove("A2")
        assert a2 not in t1.annotations
        # remove annotation by instance
        t1.annotations.remove(a3)
        assert a3 not in t1.annotations
        # test invalid remove by name
        with pytest.raises(ValueError):
            t1.annotations.remove("InvalidAnnotation")
        # test invalid remove by instance
        with pytest.raises(ValueError):
            t1.annotations.remove(int)
        # test invalid remove by instance uuid
        orig_uuid = a1._uuid
        a1._uuid = "InvalidUUID"
        with pytest.raises(ValueError):
            t1.annotations.remove(a1)
        a1._uuid = orig_uuid
        # test clear annotations (remove all)
        assert len(t1.annotations) != 0
        t1.annotations.clear()
        assert len(t1.annotations) == 0
        # test invalid remove by uuid
        with pytest.raises(ValueError):
            t1.annotations._remove_annotation_by_uuid("InvalidUUID")

    def test_annotation_ops(self, universe, session):
        t1 = mn.Trajectory(universe)
        assert len(t1.annotations) == 0
        # add annotation operator
        # first annotation type due to no invoke - must have no required inputs
        bpy.ops.mn.add_annotation("EXEC_DEFAULT", uuid=t1.uuid)
        assert len(t1.annotations) == 1
        a1 = t1.annotations[0]
        # remove annotation operator
        bpy.ops.mn.remove_annotation(
            "EXEC_DEFAULT", uuid=t1.uuid, annotation_uuid=a1._uuid
        )
        assert len(t1.annotations) == 0

    def test_trajectory_annotation_atom_info(self, universe, session):
        t1 = mn.Trajectory(universe)
        assert len(t1.annotations) == 0
        # test defaults
        t1.annotations.add_atom_info()
        assert len(t1.annotations) == 1
        t1.annotations.clear()
        # test selection string
        phrase = "all"
        a1 = t1.annotations.add_atom_info(selection=phrase)
        assert len(t1.annotations) == 1
        assert phrase == a1.selection
        t1.annotations.clear()
        # test selection atom group
        ag = universe.select_atoms("all")
        a1 = t1.annotations.add_atom_info(selection=ag)
        assert len(t1.annotations) == 1
        assert ag == a1.selection
        t1.annotations.clear()
        # test invalid selection type (not str or AtomGroup)
        with pytest.raises(ValueError):
            t1.annotations.add_atom_info(selection=1)

    def test_trajectory_annotation_com(self, universe, session):
        t1 = mn.Trajectory(universe)
        assert len(t1.annotations) == 0
        # test defaults
        t1.annotations.add_com()
        assert len(t1.annotations) == 1
        t1.annotations.clear()
        # test selection string
        phrase = "all"
        a1 = t1.annotations.add_com(selection=phrase)
        assert len(t1.annotations) == 1
        assert phrase == a1.selection
        t1.annotations.clear()
        # test selection atom group
        ag = universe.select_atoms("all")
        a1 = t1.annotations.add_com(selection=ag)
        assert len(t1.annotations) == 1
        assert ag == a1.selection
        t1.annotations.clear()
        # test invalid selection type (not str or AtomGroup)
        with pytest.raises(ValueError):
            t1.annotations.add_com(selection=1)

    def test_trajectory_annotation_com_distance(self, universe, session):
        t1 = mn.Trajectory(universe)
        assert len(t1.annotations) == 0
        # test defaults
        with pytest.raises(ValueError):
            t1.annotations.add_com_distance()
        assert len(t1.annotations) == 0
        # test both required selections
        with pytest.raises(ValueError):
            t1.annotations.add_com_distance(selection1="all")
        assert len(t1.annotations) == 0
        # test selection strings
        phrase1 = "all"
        a1 = t1.annotations.add_com_distance(selection1=phrase1, selection2="protein")
        assert len(t1.annotations) == 1
        assert phrase1 == a1.selection1
        t1.annotations.clear()
        # test selection atom groups
        ag1 = universe.select_atoms("all")
        a1 = t1.annotations.add_com_distance(selection1=ag1, selection2="protein")
        assert len(t1.annotations) == 1
        assert ag1 == a1.selection1
        t1.annotations.clear()
        # test invalid selection type (not str or AtomGroup)
        with pytest.raises(ValueError):
            t1.annotations.add_com_distance(selection1=1, selection2=2)

    def test_trajectory_annotation_canonical_dihedrals(self, universe, session):
        t1 = mn.Trajectory(universe)
        assert len(t1.annotations) == 0
        # test defaults - needs resid input
        with pytest.raises(ValueError):
            t1.annotations.add_canonical_dihedrals()
        assert len(t1.annotations) == 0
        # test invalid resid
        with pytest.raises(IndexError):
            t1.annotations.add_canonical_dihedrals(resid=1000)
        assert len(t1.annotations) == 0
        # test valid resid
        a1 = t1.annotations.add_canonical_dihedrals(resid=1)
        assert len(t1.annotations) == 1
        # test change of resid through API
        a1.resid = 2

    def test_trajectory_annotation_universe_info(self, universe, session):
        t1 = mn.Trajectory(universe)
        assert len(t1.annotations) == 0
        t1.annotations.add_universe_info()
        assert len(t1.annotations) == 1

    def test_molecule_annotation_molecule_info(self):
        mol = mn.Molecule.load(data_dir / "1cd3.cif")
        assert len(mol.annotations) == 0
        mol.annotations.add_molecule_info()
        assert len(mol.annotations) == 1

    def test_density_annotation_density_info(self, density_file):
        d1 = mn.entities.density.load(density_file)
        assert len(d1.annotations) == 0
        d1.annotations.add_density_info()
        assert len(d1.annotations) == 1

    def test_density_annotation_grid_axes(self, density_file):
        d1 = mn.entities.density.load(density_file)
        assert len(d1.annotations) == 0
        d1.annotations.add_grid_axes()
        assert len(d1.annotations) == 1

    def test_common_annotation_label_2d(self, universe, session):
        t1 = mn.Trajectory(universe)
        assert len(t1.annotations) == 0
        t1.annotations.add_label_2d(text="2D Text", location=(0.25, 0.75))
        assert len(t1.annotations) == 1

    def test_common_annotation_label_3d(self, universe, session):
        t1 = mn.Trajectory(universe)
        assert len(t1.annotations) == 0
        t1.annotations.add_label_3d(text="3D Text", location=(0.25, 0.5, 0.75))
        assert len(t1.annotations) == 1

    @pytest.mark.skip(reason="This currently fails on MacOS")
    def test_annotations_render_image(self, universe, session):
        assert "mn_annotations" not in bpy.data.images
        canvas = mn.Canvas(resolution=(192, 108))
        canvas.engine = "CYCLES"  # Only works for this
        canvas.engine.samples = 1
        t1 = mn.Trajectory(universe.select_atoms("resid 1"))
        t1.annotations.add_com(selection="resid 1")
        bpy.ops.render.render()
        assert mn.scene.compositor.annotations_image in bpy.data.images
        scene = bpy.context.scene
        assert scene.mn.auto_setup_compositor
        assert scene.node_tree
        nodes = scene.node_tree.nodes
        mn_compositor_node_name = mn.scene.compositor.mn_compositor_node_name
        assert mn_compositor_node_name in nodes
        mn_compositor_node = nodes[mn_compositor_node_name]
        assert (
            mn_compositor_node.inputs["Image"].links[0].from_node.name
            == "Render Layers"
        )
        assert mn_compositor_node.outputs["Image"].links[0].to_node.name == "Composite"
