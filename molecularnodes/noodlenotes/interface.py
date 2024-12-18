from typing import List, Union

import bpy


class InterfaceItem:
    def __init__(self, item: bpy.types.NodeTreeInterface) -> None:
        self.item: bpy.types.NodeTreeInterface = item

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
        return self.item.socket_type.replace("NodeSocket", "")

    @property
    def is_vector(self) -> bool:
        return self.type in ["Vector", "Color", "Rotation", "Matrix"]

    def __len__(self) -> int:
        if self.type == "PANEL":
            return 0
        elif self.type in ["Vector", "Rotation"]:
            return 3
        elif self.type == "Color":
            return 4
        elif self.type == "Matrix":
            return 16
        else:
            return 1

    @property
    def default(self) -> str:
        round_length: int = 3
        try:
            default_value: Union[str, float, int, list] = self.item.default_value

            if isinstance(self.item, bpy.types.NodeTreeInterfaceSocketVector):
                if self.item.default_input in ["NORMAL", "POSITION"]:
                    default_value = self.item.default_input.title()
                else:
                    default_value = [round(x, round_length) for x in default_value]
            if isinstance(self.item, bpy.types.NodeTreeInterfaceSocketColor):
                default_value = [round(x, round_length) for x in default_value]

            if isinstance(self.item, bpy.types.NodeTreeInterfaceSocketFloat):
                default_value = round(default_value, round_length)

            if isinstance(self.item, bpy.types.NodeTreeInterfaceSocketInt):
                if self.item.default_input in ["INDEX", "ID"]:
                    default_value = self.item.default_input.title()
                else:
                    default_value = int(self.item.default_value)

            if default_value == "":
                default_value = '""'

            if default_value in ["Index", "Position", "ID", "Normal"]:
                default_value = "{}::Input".format(default_value)

            return str(default_value)

        except AttributeError:
            return "_None_"

    @property
    def socket(self) -> str:
        return "{}::{}".format(self.name, self.type)

    @property
    def name(self) -> str:
        return self.item.name

    @property
    def min(self) -> str:
        round_length: int = 4
        try:
            return str(round(self.item.min_value, round_length))
        except AttributeError:
            return "_None_"

    @property
    def max(self) -> str:
        round_length: int = 4
        try:
            return str(round(self.item.max_value, round_length))
        except AttributeError:
            return "_None_"

    @property
    def description(self) -> str:
        try:
            return self.item.description
        except AttributeError:
            return ""

    def max_length(self) -> int:
        info_to_test: List[str] = [
            self.description,
            self.min,
            self.max,
            self.default,
            self.type,
        ]
        return max([len(x) for x in info_to_test if x is not None])


class InterfaceGroup:
    def __init__(self, items: List[InterfaceItem], is_output: bool = False) -> None:
        self.items: List[InterfaceItem] = items
        self._is_output: bool = is_output
        self.lengths: dict = {attr: self.get_length(attr) for attr in self.attributes}

    @property
    def attributes(self) -> List[str]:
        if self._is_output:
            return ["description", "socket"]
        else:
            return ["socket", "default", "description"]

    def sep(self) -> str:
        string = ""
        for attribute in self.attributes:
            lines = "-" * self.lengths[attribute]
            if attribute == "default":
                string += ":{}:|".format(lines)
            elif attribute == "socket" and self._is_output:
                string += "{}:|".format(lines)
            else:
                string += "{}|".format(lines)

        return "|" + string + "\n"

    def get_length(self, name: str) -> int:
        strings: List[str] = [getattr(x, name) for x in self.items]
        return max([len(x) for x in strings + [name]])

    def __len__(self) -> int:
        return len(self.items)

    def get_padded_attr(self, item: InterfaceItem, attribute: str) -> str:
        string = getattr(item, attribute)
        if attribute in ["type", "name"]:
            string = "{}::{}".format(string, item.type)

        if "::" in string:
            string = "`{}`".format(string)

        return string.ljust(self.lengths[attribute])

    def item_to_line(self, item: InterfaceItem) -> str:
        joined: str = "|".join(
            [self.get_padded_attr(item, attr) for attr in self.attributes]
        )
        return "|" + joined + "|"

    def top_line(self) -> str:
        joined: str = "|".join(
            [attr.title().ljust(self.lengths[attr]) for attr in self.attributes]
        )
        return "|" + joined + "|\n"

    def body(self) -> str:
        return "\n".join([self.item_to_line(x) for x in self.items])

    def tail(self) -> str:
        col_widths: list[int] = [10, 15, 80]
        if self._is_output:
            col_widths = [90, 10]

        widths_str = (
            "{"
            + 'tbl-colwidths="[{}]"'.format(", ".join([str(x) for x in col_widths]))
            + "}"
        )

        return "\n\n: {}\n\n".format(widths_str)

    def as_markdown(self, title: str = "", level: int = 3) -> str:
        body: str = self.body()
        if not body:
            return ""
        hashes: str = "#" * level
        lines: str = f"{hashes} {title}\n\n"
        for x in [
            self.top_line(),
            self.sep(),
            self.body(),
            self.tail(),
            "\n",
        ]:
            lines += x

        return lines

    def __repr__(self) -> str:
        return self.as_markdown()
