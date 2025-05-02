import yaml
from ..assets import ASSET_DIR
from .menu import Break, CustomItem, Menu, MenuItem, Submenu

with open(ASSET_DIR / "nodes.yml", "r") as f:
    data = yaml.safe_load(f)

submenus = []

for key, value in data.items():
    items = []
    for item in value["items"]:
        if "type" in item:
            if item["type"] == "break":
                items.append(Break())
            else:
                items.append(CustomItem(**item))
        else:
            items.append(MenuItem(**item))

    submenus.append(Submenu(name=key, title=value["title"], items=items))

menu_items = Menu(submenus=submenus)
