from abc import ABCMeta
import bpy
from pathlib import Path
from io import BytesIO, StringIO
from uuid import uuid1
from .. import blender as bl
from bpy.types import Object
from ..bpyd import (
    BlenderObject,
)


class MolecularEntity(
    BlenderObject,
    metaclass=ABCMeta,
):
    def __init__(self, file_path: str | Path) -> None:
        if type(file_path) in [BytesIO, StringIO]:
            self.file = self._read(file_path=file_path)
            self.file_path = None
        else:
            self.file_path = bl.path_resolve(file_path)
            self.file = self._read(self.file_path)
        self.uuid: str = str(uuid1())
        self._object: Object | None

    @classmethod
    def _read(self, file_path: str | Path | BytesIO):
        """
        Initially open the file, ready to extract the required data.

        Parameters
        ----------
        file_path : Union[Path, io.BytesIO]
            The path to the file which stores the atomic data.
        """
        pass
