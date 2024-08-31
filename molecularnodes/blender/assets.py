import time
from typing import List
import bpy
import re


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
    def __init__(self, timing=False):
        self.current_names: List[str] = []
        self.start_time = None
        self.timing = timing

    @property
    def trees(self) -> List:
        node_groups = list(bpy.data.node_groups)
        material_trees = [
            mat.node_tree for mat in bpy.data.materials if mat.node_tree is not None
        ]
        return node_groups + material_trees

    def tree_snapshot(self) -> None:
        self.current_names = [tree.name for tree in self.trees]

    def new_trees(self) -> List:
        return [tree for tree in self.trees if tree.name not in self.current_names]

    def __enter__(self):
        if self.timing:
            self.start_time = time.time()
        self.tree_snapshot()

    def __exit__(self, type, value, traceback):
        deduplicate_node_trees(self.new_trees())
        if self.timing:
            end_time = time.time()
            print(f"De-duplication time: {end_time - self.start_time:.2f} seconds")
