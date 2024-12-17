import bpy
import pathlib
import sys

from .interface import InterfaceGroup, InterfaceItem
from . import markdown
from typing import List

TOP_FOLDER = pathlib.Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(TOP_FOLDER))

from molecularnodes.ui.menu import Menu, MenuItem, CustomItem, Item, Break


class Documenter:
    def __init__(self, tree: bpy.types.NodeTree, menu_item: MenuItem = None) -> None:
        self.tree = tree
        self.items = [InterfaceItem(x) for x in tree.interface.items_tree]
        self.inputs = InterfaceGroup([x for x in self.items if x.is_input])
        self.outputs = InterfaceGroup(
            [x for x in self.items if x.is_output], is_output=True
        )
        self.menu_item = menu_item
        self.level = 2

    @property
    def name(self) -> str:
        return self.tree.name

    def title(self) -> str:
        return f"## {self.tree.name.removesuffix('_')}"

    def description(self) -> str:
        return self.menu_item.description + "\n\n" + self.tree.description

    def videos(self) -> List[str]:
        links = self.menu_item.videos

        if links is None:
            return None

        for x in links:
            if x is None:
                return None

        if isinstance(links, str):
            links = [links]

        if not all([isinstance(x, str) for x in links]):
            raise ValueError(f"All url values must be strings: {links=}")

        videos = "\n\n".join(
            [markdown.Video(x).as_markdown() for x in links if x is not None]
        )

        return "\n\n" + videos + "\n\n"

    def collect_items(self):
        items = [
            self.title(),
            self.description(),
            self.videos(),
            self.outputs.as_markdown("Outputs"),
            self.inputs.as_markdown("Inputs"),
        ]
        return [item for item in items if item is not None]

    def as_markdown(self) -> str:
        text = "\n"
        text += "\n\n".join(self.collect_items())
        return text


class MenuItemDocummenter(Documenter):
    def __init__(self, menu_item: MenuItem) -> None:
        super().__init__(tree=menu_item.tree, menu_item=menu_item)


class TreeDocumenter(Documenter):
    def __init__(self, tree: bpy.types.NodeTree) -> None:
        super().__init__(tree=tree)
