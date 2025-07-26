from ..annotations.base import BaseAnnotation


class Label2D(BaseAnnotation):
    """
    Show a label in 2D space (viewport / render)

    Attributes
    ----------
    text: str
        Text to be displayed

    location: tuple[float, float]
        Normalized coordinates (0.0 - 1.0) of the postion in viewport / render

    """

    annotation_type = "label_2d"

    text: str
    location: tuple[float, float] = (0.5, 0.5)

    def validate(self) -> bool:
        params = self.interface
        x, y = params.location
        if (not 0 <= x <= 1) or (not 0 <= y <= 1):
            raise ValueError("Normalized coordinates should lie between 0 and 1")
        return True

    def draw(self) -> None:
        params = self.interface
        self.draw_text_2d_norm(params.location, params.text)


class Label3D(BaseAnnotation):
    """
    Show a label in 3D space (world space)

    Attributes
    ----------
    text: str
        Text to be displayed

    location: tuple[float, float, float]
        3D Coordinates of the location to display text

    """

    annotation_type = "label_3d"

    text: str
    location: tuple[float, float, float] = (0.0, 0.0, 0.0)

    def draw(self) -> None:
        params = self.interface
        self.draw_text_3d(params.location, params.text)
