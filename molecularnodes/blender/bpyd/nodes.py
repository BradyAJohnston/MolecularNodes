import bpy
from typing import List
import re
import time
import warnings


NODE_DUP_SUFFIX = r"\.\d{3}$"


class NodeGroupCreationError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


def deduplicate_node_trees(node_trees: List[str]):
    # Compile the regex pattern for matching a suffix of a dot followed by 3 numbers
    node_duplicate_pattern = re.compile(r"\.\d{3}$")
    to_remove: List[bpy.types.GeometryNodeTree] = []

    for node_tree in node_trees:
        # Check if the node tree's name matches the duplicate pattern and is not a "NodeGroup"
        for node in node_tree.nodes:
            if not (
                hasattr(node, "node_tree")
                and node_duplicate_pattern.search(node.node_tree.name)
                and "NodeGroup" not in node.node_tree.name
            ):
                continue

            old_name = node.node_tree.name
            # Remove the numeric suffix to get the original name
            name_sans = old_name.rsplit(".", 1)[0]
            replacement = bpy.data.node_groups.get(name_sans)
            if not replacement:
                continue

            # print(f"matched {old_name} with {name_sans}")
            node.node_tree = replacement
            to_remove.append(bpy.data.node_groups[old_name])

    for tree in to_remove:
        try:
            # remove the data from the blend file
            bpy.data.node_groups.remove(tree)
        except ReferenceError:
            pass


def cleanup_duplicates(purge: bool = False):
    # Collect all node trees from node groups, excluding "NodeGroup" named ones
    node_trees = [tree for tree in bpy.data.node_groups if "NodeGroup" not in tree.name]

    # Call the deduplication function with the collected node trees
    deduplicate_node_trees(node_trees)

    if purge:
        # Purge orphan data blocks from the file
        bpy.ops.outliner.orphans_purge()


class DuplicatePrevention:
    "Context manager to cleanup duplicated node trees when appending node groups"

    def __init__(self, timing=False):
        self.current_names: List[str] = []
        self.start_time = None
        self.timing = timing

    def __enter__(self):
        self.current_names = [tree.name for tree in bpy.data.node_groups]
        if self.timing:
            self.start_time = time.time()

    def __exit__(self, type, value, traceback):
        new_trees = [
            tree for tree in bpy.data.node_groups if tree.name not in self.current_names
        ]
        deduplicate_node_trees(new_trees)
        if self.timing:
            end_time = time.time()
            print(f"De-duplication time: {end_time - self.start_time:.2f} seconds")


class MaintainConnections:
    # capture input and output links, so we can rebuild the links based on name
    # and the sockets they were connected to
    # as we collect them, remove the links so they aren't automatically connected
    # when we change the node_tree for the group

    def __init__(self, node: bpy.types.GeometryNode) -> None:
        self.node = node
        self.input_links = []
        self.output_links = []

    def __enter__(self):
        "Store all the connections in and out of this node for rebuilding on exit."
        self.node_tree = self.node.id_data

        for input in self.node.inputs:
            for input_link in input.links:
                self.input_links.append((input_link.from_socket, input.name))
                self.node_tree.links.remove(input_link)

        for output in self.node.outputs:
            for output_link in output.links:
                self.output_links.append((output.name, output_link.to_socket))
                self.node_tree.links.remove(output_link)

        try:
            self.material = self.node.inputs["Material"].default_value
        except KeyError:
            self.material = None

    def __exit__(self, type, value, traceback):
        "Rebuild the connections in and out of this node that were stored on entry."
        # rebuild the links based on names of the sockets, not their identifiers
        link = self.node_tree.links.new
        for input_link in self.input_links:
            try:
                link(input_link[0], self.node.inputs[input_link[1]])
            except KeyError:
                pass
        for output_link in self.output_links:
            try:
                link(self.node.outputs[output_link[0]], output_link[1])
            except KeyError:
                pass

        # reset all values to tree defaults
        tree = self.node.node_tree
        for item in tree.interface.items_tree:
            if item.item_type == "PANEL":
                continue
            if item.in_out == "INPUT":
                if hasattr(item, "default_value"):
                    self.node.inputs[item.identifier].default_value = item.default_value

        if self.material:
            try:
                self.node.inputs["Material"].default_value = self.material
            except KeyError:
                # the new node doesn't contain a material slot
                pass


def swap_tree(node: bpy.types.GeometryNode, tree: bpy.types.GeometryNodeTree) -> None:
    with MaintainConnections(node):
        node.node_tree = tree
        node.name = tree.name


def append_from_blend(
    name: str, filepath: str, link: bool = False
) -> bpy.types.GeometryNodeTree:
    "Append a Geometry Nodes node tree from the given .blend file"
    try:
        return bpy.data.node_groups[name]
    except KeyError:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with DuplicatePrevention():
                bpy.ops.wm.append(
                    "EXEC_DEFAULT",
                    directory=filepath,
                    filename=name,
                    link=link,
                    use_recursive=True,
                )
        return bpy.data.node_groups[name]
