import bpy

from .interface import InterfaceGroup, InterfaceItem
from . import markdown
from typing import List


class TreeDocumenter:
    def __init__(
        self,
        tree: bpy.types.NodeTree,
        description: str = None,
        videos: str | List[str] = None,
    ) -> None:
        self.tree = tree
        self.items = [InterfaceItem(x) for x in tree.interface.items_tree]
        self.inputs = InterfaceGroup([x for x in self.items if x.is_input])
        self.outputs = InterfaceGroup([x for x in self.items if x.is_output])
        self._description = description

        if videos:
            if isinstance(videos, str):
                videos = [videos]
            if not all([isinstance(x, str) for x in videos]):
                raise TypeError(f"All images must be strings: {videos=}")
            self.videos = [markdown.Video(link) for link in videos]
        else:
            self.videos = None
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

    def images(self) -> str:
        if not self.videos:
            return None
        images_string = "\n\n"
        images_string += "\n\n".join(
            [video.as_markdown() for video in self.videos if video]
        )
        images_string += "\n\n"
        return images_string

    def collect_items(self):
        items = [
            self.title(),
            self.description(),
            self.images(),
            self.inputs.as_markdown("Inputs"),
            self.outputs.as_markdown("Outputs"),
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
        tree=tree, description=menu_item.description, videos=menu_item.video_url
    )

    return doc
