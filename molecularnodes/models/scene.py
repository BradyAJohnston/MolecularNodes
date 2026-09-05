from bpy.types import Scene
from ..ui import props


class MNScene(Scene):
    def __init__(self):
        super().__init__()
        self.mn = props.MolecularNodesSceneProperties
