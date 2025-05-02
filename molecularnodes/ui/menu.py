from typing import List, Union
import bpy


class Item:
    def __init__(self) -> None:
        self.is_break = False
        self.is_custom = False
        self.backup: str = None

    @classmethod
    def menu(
        cls, layout: bpy.types.UILayout, context: bpy.types.Context = None
    ) -> None:
        pass

    @property
    def node_name(self) -> str:
        if hasattr(self, "backup"):
            return self.backup
        else:
            return self.name

    @property
    def tree(self) -> bpy.types.GeometryNodeGroup:
        if self.name.startswith("mn."):
            return bpy.data.node_groups[self.backup]
        return bpy.data.node_groups[self.name]


class MenuItem(Item):
    def __init__(
        self,
        name: str,
        label: str | None = None,
        description: str = "",
        videos: list[str] | None = None,
        backup: str | None = None,
    ) -> None:
        super().__init__()
        self.name = name
        self.label = label
        if self.label is None:
            self.label = self.name
        self.description = description
        self._videos = videos
        self.backup = backup

    def short_description(self):
        return self.description.split("\n")[0].removesuffix(".")

    @property
    def videos(self) -> list[str]:
        if self._videos is None:
            return []
        if isinstance(self._videos, str):
            return [self._videos]
        return self._videos

    def menu(
        self, layout: bpy.types.UILayout, context: bpy.types.Context = None
    ) -> None:
        if self.label is None:
            self.label = self.name

        if self.name.startswith("mn."):
            layout.operator(self.name)
            return None

        op = layout.operator("mn.add_custom_node_group", text=self.label)
        op.node_label = self.label
        op.node_name = self.name
        op.node_description = self.description
        op.node_link = False

    def to_dict(self) -> dict:
        "Return a dictionary with key: node_name and value: node_data (name, label, description, links, backup)"
        return {
            self.name: {
                "name": self.name,
                "label": self.label,
                "description": self.description,
                "videos": self.videos,
                "backup": self.backup,
            }
        }


class CustomItem(Item):
    def __init__(
        self,
        label: str,
        field: str,
        dtype: str,
        name: str,
        prefix: str,
        property_id: str,
        description: str,
        videos: str = None,
        type: str = None,
    ) -> None:
        super().__init__()
        self.label = label
        self.field = field
        self.dtype = dtype
        self.name = name
        self.prefix = prefix
        self.property_id = property_id
        self.description = description
        self.videos = videos
        self.is_custom = True

    def menu(
        self, layout: bpy.types.UILayout, context: bpy.types.Context = None
    ) -> None:
        row = layout.row()
        op = row.operator("mn.iswitch_custom", text=self.label)
        op.field = self.field
        op.dtype = self.dtype
        op.prefix = self.prefix
        op.node_property = self.property_id
        op.node_name = self.label

        if self.dtype == "RGBA":
            op.description = f"Choose custom colors for {self.label}"
        elif self.dtype == "BOOLEAN":
            op.description = f"Choose custom selections for {self.label}"
        else:
            raise ValueError(f"Data type currently not supported: {self.dtype}")
        # test if the object has the currently tested property to enable operator
        obj = context.active_object
        if obj is None:
            row.enabled = False
        else:
            row.enabled = bool(obj.get(self.property_id))


class Break:
    def __init__(self, text: str = None) -> None:
        super().__init__()
        self.is_break = True
        self.text = text

    def menu(
        self, layout: bpy.types.UILayout, context: bpy.types.Context = None
    ) -> None:
        layout.separator()
        # optionally we can add a subtitle for the next section of nodes
        if self.text and self.text.strip() != "":
            layout.label(text=self.text)


class Submenu:
    def __init__(self, name, items, title: str = None, description: str = None) -> None:
        self.name: str = name
        self.items: List[Union[MenuItem, Break, CustomItem]] = items
        self.title = title
        self.description = description

    def node_names(self):
        return [item.name for item in self.items if not isinstance(item, Break)]

    def menu(self, layout: bpy.types.UILayout, context: bpy.types.Context):
        for item in self.items:
            item.menu(layout=layout, context=context)


class Menu:
    def __init__(self, submenus: List[Submenu]) -> None:
        self.submenus = submenus

    def get_submenu(self, name: str) -> Submenu:
        for sub in self.submenus:
            if sub.name == name:
                return sub
