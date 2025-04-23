import typing
from collections import Counter, deque
import bpy
from mathutils import Vector


def contains_geo_socket(sockets: bpy.types.NodeInputs | bpy.types.NodeOutputs) -> bool:
    """
    Check if any socket in the collection is a geometry socket

    Parameters
    ----------
    sockets : bpy.types.NodeInputs | bpy.types.NodeOutputs
        Collection of node sockets to check

    Returns
    -------
    bool
        True if any socket is a geometry socket
    """
    return any([s.bl_idname == "NodeSocketGeometry" for s in sockets])


def node_has_geo_socket(node: bpy.types.GeometryNode) -> bool:
    """
    Check if a node has any geometry sockets in its inputs or outputs

    Parameters
    ----------
    node : bpy.types.GeometryNode
        The node to check

    Returns
    -------
    bool
        True if the node has at least one geometry socket
    """
    return any([contains_geo_socket(x) for x in [node.inputs, node.outputs]])


def build_dependency_graph(tree: bpy.types.NodeTree) -> tuple[dict, Counter]:
    """
    Build a graph representing node dependencies and count input connections

    Parameters
    ----------
    tree : bpy.types.NodeTree
        The node tree to analyze

    Returns
    -------
    tuple
        Contains:
            dict
                Mapping of nodes to their dependent nodes
            Counter
                Count of connections for each socket
    """
    dependency_graph = {node: set() for node in tree.nodes}
    socket_input_connection_count = Counter()

    # populate the graph based on node connections
    for link in tree.links:
        dependency_graph[link.from_node].add(link.to_node)
        socket_input_connection_count[link.to_socket] += 1

    return dependency_graph, socket_input_connection_count


def topological_sort(dependency_graph: dict) -> list:
    """
    Sort nodes by their dependencies using a topological sort algorithm

    Parameters
    ----------
    dependency_graph : dict
        Mapping of nodes to their dependent nodes

    Returns
    -------
    list
        Nodes sorted in topological order
    """
    # count incoming connections for each node
    incoming_connection_count = {node: 0 for node in dependency_graph}
    for source_node in dependency_graph:
        for target_node in dependency_graph[source_node]:
            incoming_connection_count[target_node] += 1

    # start with nodes that have no dependencies
    processing_queue = deque(
        [
            node
            for node in incoming_connection_count
            if incoming_connection_count[node] == 0
        ]
    )
    sorted_node_order = []

    # process nodes in breadth-first order
    while processing_queue:
        current_node = processing_queue.popleft()
        sorted_node_order.append(current_node)

        # update counts for nodes that depend on the current node
        for dependent_node in dependency_graph[current_node]:
            incoming_connection_count[dependent_node] -= 1
            # if all dependencies are processed, add to queue
            if incoming_connection_count[dependent_node] == 0:
                processing_queue.append(dependent_node)

    return sorted_node_order


def organize_into_columns(nodes_in_order: list, dependency_graph: dict) -> list:
    """
    Organize nodes into columns based on their dependencies

    Parameters
    ----------
    nodes_in_order : list
        Nodes sorted in topological order
    dependency_graph : dict
        Mapping of nodes to their dependent nodes

    Returns
    -------
    list
        Columns of nodes, where each column is a list of nodes
    """
    node_columns = []
    node_column_assignment = {}

    for node in reversed(nodes_in_order):
        # node goes in column after its furthest dependent
        node_column_assignment[node] = (
            max(
                [
                    node_column_assignment[dependent_node]
                    for dependent_node in dependency_graph[node]
                ],
                default=-1,
            )
            + 1
        )

        # add node to its assigned column
        if node_column_assignment[node] == len(node_columns):
            node_columns.append([node])
        else:
            node_columns[node_column_assignment[node]].append(node)

    # reverse columns to get left-to-right flow
    return list(reversed(node_columns))


def calculate_node_dimensions(
    node: bpy.types.Node, socket_input_connection_count: Counter, interface_scale: float
) -> tuple[float, float]:
    """
    Calculate the visual dimensions of a node

    Parameters
    ----------
    node : bpy.types.Node
        The node to calculate dimensions for
    socket_input_connection_count : Counter
        Counter of connections for each socket
    interface_scale : float
        UI scale factor

    Returns
    -------
    tuple[float, float]
        Width and height of the node
    """
    # height constants for different node elements
    node_header_height = 20
    node_socket_height = 28
    node_property_row_height = 28
    node_vector_input_height = 84

    # count enabled inputs and outputs
    enabled_input_count = len(
        list(filter(lambda input_socket: input_socket.enabled, node.inputs))
    )
    enabled_output_count = len(
        list(filter(lambda output_socket: output_socket.enabled, node.outputs))
    )

    # get properties specific to this node type (not inherited)
    inherited_property_ids = [
        property.identifier
        for base_class in type(node).__bases__
        for property in base_class.bl_rna.properties
    ]

    node_specific_property_count = len(
        [
            property
            for property in node.bl_rna.properties
            if property.identifier not in inherited_property_ids
        ]
    )

    # count vector inputs that need UI widgets (not connected)
    unconnected_vector_input_count = len(
        list(
            filter(
                lambda input_socket: input_socket.enabled
                and input_socket.type == "VECTOR"
                and socket_input_connection_count[input_socket] == 0,
                node.inputs,
            )
        )
    )

    # calculate total node height based on components
    total_node_height = (
        node_header_height
        + (enabled_output_count * node_socket_height)
        + (node_specific_property_count * node_property_row_height)
        + (enabled_input_count * node_socket_height)
        + (unconnected_vector_input_count * node_vector_input_height)
    ) * interface_scale

    return node.width, total_node_height


def position_nodes_in_columns(
    node_columns: list,
    socket_input_connection_count: Counter,
    spacing: typing.Tuple[float, float] = (50, 25),
) -> None:
    """Position nodes in columns with appropriate spacing

    Parameters
    ----------
    node_columns : list
        List of columns, where each column is a list of nodes
    socket_input_connection_count : Counter
        Counter of connections for each socket
    spacing : tuple of float, optional
        Tuple of (horizontal, vertical) spacing between nodes, by default (50, 25)
    """
    interface_scale = bpy.context.preferences.view.ui_scale
    non_geo_offset = 20 + 28 * 2  # header + 2 socket heights

    # position nodes column by column
    position_x = 0
    for column in node_columns:
        widest_node_in_column = 0
        position_y = 0

        for node in column:
            node.update()

            width, height = calculate_node_dimensions(
                node, socket_input_connection_count, interface_scale
            )

            # track widest node for column spacing
            if width > widest_node_in_column:
                widest_node_in_column = width

            # position node
            node.location = (position_x, position_y)

            # adjust position for non-geometry nodes
            if not node_has_geo_socket(node):
                node.location -= Vector((0, non_geo_offset))

            # move down for next node with spacing
            position_y -= height + spacing[1]

        # move right for next column with spacing
        position_x += widest_node_in_column + spacing[0]


def position_special_nodes(
    tree: bpy.types.NodeTree, vertical_offset: float = 100
) -> None:
    """Position special nodes like Group Input and Group Output at the top

    Parameters
    ----------
    tree : bpy.types.NodeTree
        The node tree to modify
    vertical_offset : float, optional
        Vertical offset above the highest node, by default 100

    Returns
    -------
    None
    """
    highest_y_position = max([node.location[1] for node in tree.nodes])
    for special_node_name in ["Group Input", "Group Output"]:
        if special_node_name in tree.nodes:
            special_node = tree.nodes[special_node_name]
            special_node.location = (
                special_node.location[0],
                highest_y_position + vertical_offset,
            )


def cleanup_orphaned_nodes(tree: bpy.types.NodeTree, max_iter: int = 100) -> None:
    """Remove nodes that are not connected to anything"""
    to_remove = []
    for _ in range(max_iter):
        for node in tree.nodes:
            if len(node.outputs) == 0:
                continue
            if not any([s.is_linked for s in node.outputs]):
                to_remove.append(node)

        if len(to_remove) == 0:
            break

        for node in to_remove:
            tree.nodes.remove(node)

        to_remove = []

    if "Group Input" not in tree.nodes:
        n_input = tree.nodes.new("NodeGroupInput")
        if "Join Geometry" in tree.nodes:
            tree.links.new(n_input.outputs[0], tree.nodes["Join Geometry"].inputs[0])


def arrange_tree(
    tree: bpy.types.NodeTree,
    spacing: typing.Tuple[float, float] = (50, 25),
) -> None:
    """Arrange nodes in a node tree based on their dependencies

    Parameters
    ----------
    tree : bpy.types.GeometryNodeTree
        The node tree to arrange
    spacing : tuple of float
        Tuple of (horizontal, vertical) spacing between nodes

    Returns
    -------
    None
        This function modifies the node tree in place

    Notes
    -----
    This function organizes nodes into columns and positions them with appropriate spacing.
    Nodes are arranged from left to right based on their dependencies, with special
    handling for geometry nodes and group input/output nodes.
    """
    cleanup_orphaned_nodes(tree)
    dependency_graph, socket_input_connection_count = build_dependency_graph(tree)
    nodes_in_dependency_order = topological_sort(dependency_graph)
    node_columns = organize_into_columns(nodes_in_dependency_order, dependency_graph)
    position_nodes_in_columns(node_columns, socket_input_connection_count, spacing)
    position_special_nodes(tree)
