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

    def test_trajectory_annotations_registration(self):
        class TestAnnotation(mn.entities.trajectory.TrajectoryAnnotation):
            annotation_type = "test_annotation"

            def draw(self):
                pass

        manager = mn.entities.trajectory.TrajectoryAnnotationManager
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
        # test annotation get by name
        ga1 = t1.annotations.get("Annotation")
        assert ga1 == a1
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
        # test clear annotations (remove all)
        assert len(t1.annotations) != 0
        t1.annotations.clear()
        assert len(t1.annotations) == 0

    def test_annotation_ops(self, universe, session):
        t1 = session.add_trajectory(universe)
        assert len(t1.annotations) == 0
        # add annotation operator
        bpy.ops.mn.add_annotation("EXEC_DEFAULT", uuid=t1.uuid, type="atom_info")
        assert len(t1.annotations) == 1
        a1 = t1.annotations[0]
        # remove annotation operator
        bpy.ops.mn.remove_annotation(
            "EXEC_DEFAULT", uuid=t1.uuid, annotation_uuid=a1._uuid
        )
        assert len(t1.annotations) == 0
