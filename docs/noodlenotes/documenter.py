import bpy

from .interface import InterfaceGroup, InterfaceItem
from . import markdown


class TreeDocumenter:
    def __init__(
        self, tree: bpy.types.NodeTree, description: str = None, image: str = None
    ) -> None:
        self.tree = tree
        self.items = [InterfaceItem(x) for x in tree.interface.items_tree]
        self.inputs = InterfaceGroup([x for x in self.items if x.is_input])
        self.outputs = InterfaceGroup([x for x in self.items if x.is_output])
        self._description = description
        self.image = markdown.Video(image)
        self.level = 2

    @property
    def name(self) -> str:
        return self.tree.name

    def title(self) -> str:
        return f"## {self.tree.name.removesuffix('_')}"

    def description(self) -> str:
        if self._description:
            return self._description + "\n\n" + self.tree.description

        return self.tree.description

    def collect_items(self):
        items = [
            self.title(),
            self.description(),
            self.image.as_markdown(),
            "### Inputs",
            self.inputs.as_markdown(),
            "### Outputs",
            self.outputs.as_markdown(),
        ]
        return [item for item in items if item is not None]

    def as_markdown(self) -> str:
        text = "\n"
        text += "\n\n".join(self.collect_items())
        return text


def MenuItemDocumenter(menu_item) -> TreeDocumenter:
    if menu_item.backup:
        tree = bpy.data.node_groups[menu_item.backup]
    else:
        tree = bpy.data.node_groups[menu_item.name]

    doc = TreeDocumenter(
        tree=tree, description=menu_item.description, image=menu_item.video_url
    )

    return doc
