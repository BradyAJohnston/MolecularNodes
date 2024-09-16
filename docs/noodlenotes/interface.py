from typing import List

import bpy


class InterfaceItem:
    def __init__(self, item: bpy.types.NodeTreeInterface) -> None:
        self.item = item

    @property
    def is_socket(self) -> bool:
        return self.item.item_type == "SOCKET"

    @property
    def is_panel(self) -> bool:
        return self.item.item_type == "PANEL"

    @property
    def is_input(self) -> bool:
        if self.is_panel:
            return False
        return self.item.in_out == "INPUT"

    @property
    def is_output(self) -> bool:
        if self.is_panel:
            return False
        return self.item.in_out == "OUTPUT"

    @property
    def type(self) -> str:
        if self.is_panel:
            return "PANEL"
        return "`{}`".format(self.item.socket_type.replace("NodeSocket", ""))

    @property
    def default(self, round_length: int = 3) -> str:
        try:
            default = self.item.default_value

            if isinstance(self.item, bpy.types.NodeTreeInterfaceSocketVector):
                if self.item.default_input in ["NORMAL", "POSITION"]:
                    default = self.item.default_input.title()
                else:
                    default = [round(x, round_length) for x in default]
            if isinstance(self.item, bpy.types.NodeTreeInterfaceSocketColor):
                default = [round(x, round_length) for x in default]

            if isinstance(self.item, bpy.types.NodeTreeInterfaceSocketFloat):
                default = round(default, round_length)

            if isinstance(self.item, bpy.types.NodeTreeInterfaceSocketInt):
                if self.item.default_input in ["INDEX", "ID"]:
                    default = self.item.default_input.title()
                else:
                    default = int(self.item.default)

            if default == "":
                default = "_None_"

            return "`{}`".format(default)

        except AttributeError:
            return "_None_"

    @property
    def name(self):
        return self.item.name

    @property
    def min(self, round_length=4):
        try:
            return "`{}`".format(round(self.item.min_value, round_length))
        except AttributeError:
            return "_None_"

    @property
    def max(self, round_length=4):
        try:
            return "`{}`".format(round(self.item.max_value, round_length))
        except AttributeError:
            return "_None_"

    @property
    def description(self):
        try:
            return self.item.description
        except AttributeError:
            return ""

    def max_length(self):
        info_to_test = [self.description, self.min, self.max, self.default]
        return max([len(x) for x in info_to_test if x is not None])

    def formatted(self):
        text = f"Default Value: {self.default}\n"
        try:
            text += f"Min: {self.min}\n"
            text += f"Max: {self.max}\n"
        except AttributeError:
            pass

        return text


class InterfaceGroup:
    def __init__(self, items: List[InterfaceItem]) -> None:
        self.items = items
        self.attributes = ["name", "type", "description", "default"]  # , "min", "max"]
        self.lengths = {attr: self.get_length(attr) for attr in self.attributes}

    def sep(self) -> int:
        text = ""
        for length in self.lengths.values():
            text += "|" + "-" * length

        return text + ":|\n"

    def get_length(self, name: str) -> int:
        strings = [getattr(x, name) for x in self.items]
        return max([len(x) for x in strings + [name]])

    def __len__(self) -> int:
        return len(self.items)

    def get_padded_attr(self, item: InterfaceItem, attribute: str) -> str:
        return f"{getattr(item, attribute).ljust(self.lengths[attribute])}"

    def item_to_line(self, item):
        joined = "|".join(
            [self.get_padded_attr(item, attr) for attr in self.attributes]
        )
        return "|" + joined + "|"

    def top_line(self):
        joined = "|".join(
            [attr.title().ljust(self.lengths[attr]) for attr in self.attributes]
        )
        return "|" + joined + "|\n"

    def body(self) -> str:
        return "\n".join([self.item_to_line(x) for x in self.items])

    def tail(self) -> str:
        return '\n\n: {tbl-colwidths="[15, 10, 55, 20]"}\n\n'

    def as_markdown(self):
        return self.top_line() + self.sep() + self.body() + self.tail() + "\n"

    def __repr__(self) -> str:
        return self.as_markdown()
