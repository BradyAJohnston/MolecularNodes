import bpy
import MDAnalysis as mda
import pytest
import molecularnodes as mn
from .constants import data_dir


class TestAnnotations:
    @pytest.fixture(scope="module")
    def universe(self):
        topo = data_dir / "md_ppr/box.gro"
        traj = data_dir / "md_ppr/first_5_frames.xtc"
        return mda.Universe(topo, traj)

    @pytest.fixture(scope="module")
    def session(self):
        return mn.session.get_session()

    def test_registered_trajectory_annotations(self):
        # test if trajectory annotations are correctly auto registered
        manager = mn.entities.trajectory.TrajectoryAnnotationManager
        # atom_info
        assert "atom_info" in manager._classes
        assert hasattr(manager, "add_atom_info")
        assert callable(getattr(manager, "add_atom_info"))

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
        t1 = session.add_trajectory(universe)
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
        t1 = session.add_trajectory(universe)
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
        t1 = session.add_trajectory(universe)
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
