# Type aliases for node inputs
LINKABLE = "NodeSocket | NodeBuilder | Any"
TYPE_INPUT_VECTOR = "tuple[float, float, float] | NodeSocket | NodeBuilder | None"
TYPE_INPUT_ROTATION = (
    "tuple[float, float, float, float] | NodeSocket | NodeBuilder | None"
)
TYPE_INPUT_BOOLEAN = "bool | NodeSocket | NodeBuilder | None"


class DataTypes:
    FLOAT = "FLOAT"
    INT = "INT"
    BOOL = "BOOL"
    VECTOR = "FLOAT_VECTOR"
    ROTATION = "ROTATION"
    COLOR = "RGBA"
    RGBA = "RGBA"
