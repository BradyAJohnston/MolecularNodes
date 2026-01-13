"""
Tests for the node builder API.

Tests the TreeBuilder, NodeBuilder, and the >> operator chaining system.
"""

from math import pi
import bpy
import pytest
from numpy.testing import assert_allclose
from molecularnodes.nodes import generated as n
from molecularnodes.nodes import sockets
from molecularnodes.nodes.builder import TreeBuilder
from molecularnodes.nodes.generated.input import Boolean
from molecularnodes.nodes.generated.manually_specified import BooleanMath


class TestTreeBuilder:
    """Tests for TreeBuilder basic functionality."""

    def test_tree_creation(self):
        """Test creating a basic node tree."""
        tree = TreeBuilder("TestTree")
        assert tree.tree is not None
        assert tree.tree.name == "TestTree"
        assert isinstance(tree.tree, bpy.types.GeometryNodeTree)

    def test_interface_definition(self):
        """Test defining tree interface with socket types."""
        tree = TreeBuilder("InterfaceTest")
        tree.interface(
            inputs=[
                sockets.SocketGeometry(name="Geometry"),
                sockets.SocketBoolean(name="Selection", default=True),
                sockets.SocketVector(name="Offset", default=(1.0, 2.0, 3.0)),
            ],
            outputs=[
                sockets.SocketGeometry(name="Geometry"),
                sockets.SocketInt(name="Count"),
            ],
        )

        # Check inputs were created
        input_names = [
            socket.name
            for socket in tree.tree.interface.items_tree
            if socket.in_out == "INPUT"
        ]
        assert "Geometry" in input_names
        assert "Selection" in input_names
        assert "Offset" in input_names

        # Check outputs were created
        output_names = [
            socket.name
            for socket in tree.tree.interface.items_tree
            if socket.in_out == "OUTPUT"
        ]
        assert "Geometry" in output_names
        assert "Count" in output_names

    def test_socket_defaults(self):
        """Test that socket defaults are set correctly."""
        tree = TreeBuilder("DefaultsTest")
        tree.interface(
            inputs=[
                sockets.SocketBoolean(name="Selection", default=True),
                sockets.SocketVector(name="Offset", default=(1.0, 2.0, 3.0)),
                sockets.SocketFloat(
                    name="Scale", default=2.5, min_value=0.0, max_value=10.0
                ),
            ],
            outputs=[sockets.SocketGeometry(name="Geometry")],
        )

        # Find the selection input
        selection_socket = None
        for item in tree.tree.interface.items_tree:
            if item.name == "Selection" and item.in_out == "INPUT":
                selection_socket = item
                break

        assert selection_socket is not None
        assert selection_socket.default_value is True


class TestContextManager:
    """Tests for the context manager functionality."""

    def test_context_manager_basic(self):
        """Test using the tree as a context manager."""
        tree = TreeBuilder("ContextTest")
        tree.interface(
            inputs=[sockets.SocketGeometry(name="Geometry")],
            outputs=[sockets.SocketGeometry(name="Geometry")],
        )

        with tree:
            # Should be able to create nodes without passing tree
            pos = n.Position()
            assert pos.node is not None
            assert pos.tree == tree

    def test_context_manager_node_creation(self):
        """Test that nodes created in context use the active tree."""
        tree = TreeBuilder("NodeCreationTest")
        tree.interface(
            inputs=[sockets.SocketGeometry(name="Geometry")],
            outputs=[sockets.SocketGeometry(name="Geometry")],
        )

        with tree:
            node1 = n.Position()
            node2 = n.SetPosition()

            # Both nodes should be in the same tree
            assert node1.tree == tree
            assert node2.tree == tree
            assert node1.node.id_data == tree.tree
            assert node2.node.id_data == tree.tree


class TestOperatorChaining:
    """Tests for the >> operator chaining."""

    def test_basic_chaining(self):
        """Test basic node chaining with >> operator."""
        tree = TreeBuilder("ChainingTest")
        tree.interface(
            inputs=[sockets.SocketGeometry(name="Geometry")],
            outputs=[sockets.SocketGeometry(name="Geometry")],
        )

        with tree:
            pos = n.Position()
            set_pos = n.SetPosition()

            # Chain with >> operator
            result = pos >> set_pos

            # Should return the right-hand node
            assert result == set_pos

            # Should create a link between the nodes
            links = tree.tree.links
            assert len(links) > 0

    def test_multi_node_chaining(self):
        """Test chaining multiple nodes together."""
        tree = TreeBuilder("MultiChainTest")
        tree.interface(
            inputs=[sockets.SocketGeometry(name="Geometry")],
            outputs=[sockets.SocketGeometry(name="Geometry")],
        )

        with tree:
            # Chain multiple nodes
            _ = (
                tree.inputs.geometry
                >> n.SetPosition()
                >> n.TransformGeometry(translation=(0, 0, 1))
                >> tree.outputs.geometry
            )

        # Check that links were created
        assert len(tree.tree.links) >= 3


class TestNamedSocketAccess:
    """Tests for named socket access (tree.inputs.socket_name)."""

    def test_input_socket_access(self):
        """Test accessing input sockets by name."""
        tree = TreeBuilder("InputAccessTest")
        tree.interface(
            inputs=[
                sockets.SocketGeometry(name="Geometry"),
                sockets.SocketBoolean(name="Selection"),
                sockets.SocketVector(name="My Offset"),
            ],
            outputs=[sockets.SocketGeometry(name="Geometry")],
        )

        with tree:
            # Access sockets by normalized name
            geo = tree.inputs.geometry
            sel = tree.inputs.selection
            offset = tree.inputs.my_offset

            assert geo is not None
            assert sel is not None
            assert offset is not None

    def test_output_socket_access(self):
        """Test accessing output sockets by name."""
        tree = TreeBuilder("OutputAccessTest")
        tree.interface(
            inputs=[sockets.SocketGeometry(name="Geometry")],
            outputs=[
                sockets.SocketGeometry(name="Geometry"),
                sockets.SocketInt(name="Count"),
                sockets.SocketFloat(name="Total Area"),
            ],
        )

        with tree:
            # Access output sockets
            geo = tree.outputs.geometry
            count = tree.outputs.count
            area = tree.outputs.total_area

            assert geo is not None
            assert count is not None
            assert area is not None


class TestExamples:
    """Tests that replicate the original example functions."""

    def test_example_basic(self):
        """Test the basic example from example.py."""
        tree = TreeBuilder("ExampleTree")

        tree.interface(
            inputs=[
                sockets.SocketGeometry(name="Geometry"),
                sockets.SocketBoolean(name="Selection", default=True),
                sockets.SocketVector(name="Offset", default=(1.0, 2.0, 3.0)),
            ],
            outputs=[
                sockets.SocketGeometry(name="Geometry"),
            ],
        )

        with tree:
            _ = (
                tree.inputs.geometry
                >> n.SetPosition(position=n.Position())
                >> n.TransformGeometry(translation=(0, 0, 1))
                >> tree.outputs.geometry
            )

        # Verify tree was created correctly
        assert tree.tree is not None
        assert len(tree.tree.nodes) > 0
        assert len(tree.tree.links) > 0

    def test_example_multi_socket(self):
        """Test the multi-socket example from example.py."""
        tree = TreeBuilder("MultiSocketExample")

        tree.interface(
            inputs=[
                sockets.SocketGeometry(name="Geometry"),
                sockets.SocketBoolean(name="Selection", default=True),
            ],
            outputs=[
                sockets.SocketGeometry(name="Geometry"),
                sockets.SocketInt(name="Count"),
            ],
        )

        with tree:
            # Access multiple named sockets
            _ = (
                tree.inputs.geometry
                >> n.SetPosition(selection=tree.inputs.selection)
                >> tree.outputs.geometry
            )

        # Verify tree structure
        assert tree.tree is not None
        assert len(tree.tree.nodes) > 0

        # Check that both outputs were created
        output_names = [
            socket.name
            for socket in tree.tree.interface.items_tree
            if socket.in_out == "OUTPUT"
        ]
        assert "Geometry" in output_names
        assert "Count" in output_names


class TestGeneratedNodes:
    """Tests for generated node classes."""

    def test_position_node(self):
        """Test the Position input node."""
        tree = TreeBuilder("PositionTest")
        tree.interface(
            inputs=[sockets.SocketGeometry(name="Geometry")],
            outputs=[sockets.SocketGeometry(name="Geometry")],
        )

        with tree:
            pos = n.Position()
            assert pos.node is not None
            assert pos.node.bl_idname == "GeometryNodeInputPosition"

    def test_set_position_node(self):
        """Test the SetPosition node with parameters."""
        tree = TreeBuilder("SetPositionTest")
        tree.interface(
            inputs=[sockets.SocketGeometry(name="Geometry")],
            outputs=[sockets.SocketGeometry(name="Geometry")],
        )

        with tree:
            pos = n.Position()
            set_pos = n.SetPosition(position=pos)

            assert set_pos.node is not None
            assert set_pos.node.bl_idname == "GeometryNodeSetPosition"

            # Check that the position input was linked
            assert len(set_pos.node.inputs["Position"].links) > 0

    def test_transform_geometry_node(self):
        """Test the TransformGeometry node."""
        tree = TreeBuilder("TransformTest")
        tree.interface(
            inputs=[sockets.SocketGeometry(name="Geometry")],
            outputs=[sockets.SocketGeometry(name="Geometry")],
        )

        with tree:
            transform = n.TransformGeometry(translation=(1, 2, 3))

            assert transform.node is not None
            assert transform.node.bl_idname == "GeometryNodeTransform"

    def test_node_output_properties(self):
        """Test that output properties are accessible."""
        tree = TreeBuilder("OutputPropsTest")
        tree.interface(
            inputs=[sockets.SocketGeometry(name="Geometry")],
            outputs=[sockets.SocketGeometry(name="Geometry")],
        )

        with tree:
            bbox = n.BoundingBox()

            # Test output property accessors
            assert hasattr(bbox, "o_bounding_box")
            assert hasattr(bbox, "o_min")
            assert hasattr(bbox, "o_max")

            # They should return sockets
            assert bbox.o_bounding_box is not None
            assert bbox.o_min is not None
            assert bbox.o_max is not None


class TestComplexWorkflow:
    """Test more complex node tree workflows."""

    def test_branching_workflow(self):
        """Test a workflow with branching node connections."""
        tree = TreeBuilder("BranchingTest")
        tree.interface(
            inputs=[sockets.SocketGeometry(name="Geometry")],
            outputs=[sockets.SocketGeometry(name="Geometry")],
        )

        with tree:
            pos = n.Position()

            # Use the same position node in multiple places
            _set_pos1 = n.SetPosition(position=pos, offset=(1, 0, 0))
            _set_pos2 = n.SetPosition(position=pos, offset=(0, 1, 0))

            # Both should reference the same position node
            assert len(pos.node.outputs[0].links) == 2

    def test_multiple_inputs_workflow(self):
        """Test using multiple tree inputs in a workflow."""
        tree = TreeBuilder("MultiInputTest")
        tree.interface(
            inputs=[
                sockets.SocketGeometry(name="Geometry"),
                sockets.SocketBoolean(name="Selection"),
                sockets.SocketVector(name="Translation"),
            ],
            outputs=[sockets.SocketGeometry(name="Geometry")],
        )

        with tree:
            _ = (
                tree.inputs.geometry
                >> n.SetPosition(selection=tree.inputs.selection)
                >> n.TransformGeometry(translation=tree.inputs.translation)
                >> tree.outputs.geometry
            )

        # Verify all inputs are used
        group_input_node = tree.tree.nodes.get("Group Input")
        assert group_input_node is not None

        # Check that inputs have outgoing links
        assert len(group_input_node.outputs["Geometry"].links) > 0
        assert len(group_input_node.outputs["Selection"].links) > 0
        assert len(group_input_node.outputs["Translation"].links) > 0


def create_tree_chain():
    tree = TreeBuilder("MathTest")
    tree.interface(
        inputs=[sockets.SocketFloat(name="Value")],
        outputs=[sockets.SocketFloat(name="Result")],
    )

    with tree:
        _ = (
            tree.inputs.value
            >> n.Math.add(..., 0.1)
            >> n.VectorMath.multiply(..., (2.0, 2.0, 2.0))
            >> tree.outputs.result
        )

    return tree


def create_tree():
    tree = TreeBuilder("MathTest")
    tree.interface(
        inputs=[sockets.SocketFloat(name="Value")],
        outputs=[sockets.SocketFloat(name="Result")],
    )

    with tree:
        final = n.VectorMath.multiply(
            n.Math.add(tree.inputs.value, 0.1), (2.0, 2.0, 2.0)
        )

        final >> tree.outputs.result

    return tree


@pytest.mark.parametrize("maker", [create_tree_chain, create_tree])
def test_math_nodes(maker):
    """Test math nodes."""
    tree = maker()
    # Verify all inputs are used
    node_input = tree.tree.nodes.get("Group Input")
    assert node_input is not None

    # check the default values have been property set
    assert_allclose(tree.tree.nodes["Math"].inputs[0].default_value, 0.5)
    assert_allclose(tree.tree.nodes["Math"].inputs[1].default_value, 0.1)
    assert_allclose(tree.nodes["Vector Math"].inputs[0].default_value, (0.0, 0.0, 0.0))
    assert_allclose(tree.nodes["Vector Math"].inputs[1].default_value, (2.0, 2.0, 2.0))

    # Check that inputs have outgoing links
    assert len(node_input.outputs["Value"].links) == 1
    assert len(tree.tree.nodes.get("Group Output").inputs["Result"].links) == 1

    assert (
        tree.tree.nodes["Math"].inputs[0].links[0].from_node == tree.inputs.value.node
    )


def test_nodes():
    tree = TreeBuilder()
    tree.interface(outputs=[sockets.SocketGeometry("Geometry")])

    with tree:
        _ = (
            n.Points(1_000, position=n.RandomValue.vector())
            >> n.PointsToCurves(curve_group_id=n.RandomValue.integer(min=0, max=10))
            >> n.CurveToMesh(profile_curve=n.CurveCircle(12, radius=0.1))
            >> tree.outputs.geometry
        )


def test_mix_node():
    tree = TreeBuilder()
    tree.interface(
        inputs=[sockets.SocketInt("Count", min_value=0, max_value=100, default=50)],
        outputs=[sockets.SocketGeometry("Instances")],
    )

    with tree:
        rotation = n.Mix.rotation(
            n.RandomValue.vector((-pi, -pi, -pi), (pi, pi, pi)),
            (0, 0, 1),
            factor=n.RandomValue.float(seed=n.Index()),
        )

        selection = (
            n.RandomValue.boolean(probability=0.3)
            >> n.BooleanMath.l_not()
            >> n.BooleanMath.l_and(n.RandomValue.boolean(probability=0.8))
            >> n.BooleanMath.l_or(n.RandomValue.boolean(probability=0.5))
            >> n.BooleanMath.l_equal(n.RandomValue.boolean(probability=0.4))
            >> n.BooleanMath.l_not()
        )

        _ = (
            n.Points(tree.inputs.count, position=n.RandomValue.vector())
            >> n.InstanceOnPoints(
                selection=selection,
                instance=n.Cube(),
                rotation=rotation,
            )
            >> n.TranslateInstances(translation=(0.0, 0.1, 0.0))
            >> tree.outputs.instances
        )

    # some nodes with different data types have a different output for each data type
    # so for rotation the socket is the 4th output - this will change in the future
    # with raibow sockets eventually
    assert len(rotation.node.outputs[3].links) == 1
    assert len(tree.nodes) == 17


def test_warning_innactive_socket():
    "Raises an error because we want to not let a user silently link sockets that won't do anything"
    with TreeBuilder():
        pos = n.Position()
        mix = n.Mix.vector()
        # this works because by default we link to the currently active vector sockets
        pos >> mix
        # this now fails because we try to link to the innactive float sockets
        mix._default_input_id = "A_Float"
        with pytest.raises(RuntimeError):
            pos >> mix
