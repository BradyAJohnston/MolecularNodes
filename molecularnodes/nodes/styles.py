"""Constant style definitions for molecular visualization in Blender.

This module contains predefined style settings for different molecular visualization
modes including ball-and-stick, cartoon, ribbon, spheres, sticks and surface
representations.

"""

from dataclasses import dataclass, replace, field, fields
from typing import List, Tuple, Union


@dataclass(frozen=True)
class StyleBase:
    def get_by_key(self, original_key: str):
        for f in fields(self):
            if f.metadata.get('key') == original_key:
                return getattr(self, f.name)
        return None

    def update_properties(self, **changes):
        return replace(self, **changes)

@dataclass(frozen=True)
class BallStickStyle(StyleBase):
    style: str = field(default="ball+stick", metadata={"key": "Style"})
    quality: int = field(default=2, metadata={"key": "Quality"})
    as_mesh: bool = field(default=True, metadata={"key": "As Mesh"})
    sphere_radii: float = field(default=0.3, metadata={"key": "Sphere Radii"})
    bond_split: bool = field(default=False, metadata={"key": "Bond Split"})
    bond_find: bool = field(default=True, metadata={"key": "Bond Find"})
    bond_radius: float = field(default=0.3, metadata={"key": "Bond Radius"})
    color_blur: bool = field(default=False, metadata={"key": "Color Blur"})
    shade_smooth: bool = field(default=True, metadata={"key": "Shade Smooth"})


@dataclass(frozen=True)
class CartoonStyle(StyleBase):
    style: str = field(default="cartoon", metadata={"key": "Style"})
    quality: int = field(default=2, metadata={"key": "Quality"})
    dssp: bool = field(default=False, metadata={"key": "DSSP"})
    cylinders: bool = field(default=False, metadata={"key": "Cylinders"})
    arrows: bool = field(default=True, metadata={"key": "Arrows"})
    rounded: bool = field(default=False, metadata={"key": "Rounded"})
    thickness: float = field(default=0.6, metadata={"key": "Thickness"})
    width: float = field(default=2.2, metadata={"key": "Width"})
    loop_radius: float = field(default=0.3, metadata={"key": "Loop Radius"})
    smoothing: float = field(default=0.5, metadata={"key": "Smoothing"})
    color_blur: bool = field(default=False, metadata={"key": "Color Blur"})
    shade_smooth: bool = field(default=True, metadata={"key": "Shade Smooth"})


@dataclass(frozen=True)
class RibbonStyle(StyleBase):
    style: str = field(default="ribbon", metadata={"key": "Style"})
    quality: int = field(default=3, metadata={"key": "Quality"})
    radius: float = field(default=1.6, metadata={"key": "Radius"})
    smoothing: float = field(default=0.6, metadata={"key": "Smoothing"})
    color_blur: bool = field(default=False, metadata={"key": "Color Blur"})
    shade_smooth: bool = field(default=False, metadata={"key": "Shade Smooth"})

@dataclass(frozen=True)
class SpheresStyle(StyleBase):
    style: str = field(default="sphere", metadata={"key": "Style"})
    as_mesh: bool = field(default=True, metadata={"key": "Sphere As Mesh"})
    radii: float = field(default=0.8, metadata={"key": "Sphere Radii"})
    subdivisions: int = field(default=2, metadata={"key": "Sphere Subdivisions"})
    shade_smooth: bool = field(default=False, metadata={"key": "Shade Smooth"})


@dataclass(frozen=True)
class SticksStyle(StyleBase):
    style: str = field(default="sticks", metadata={"key": "Style"})
    quality: int = field(default=2, metadata={"key": "Quality"})
    radius: float = field(default=0.2, metadata={"key": "Radius"})
    color_blur: bool = field(default=False, metadata={"key": "Color Blur"})
    shade_smooth: bool = field(default=False, metadata={"key": "Shade Smooth"})


@dataclass(frozen=True)
class SurfaceStyle(StyleBase):
    style: str = field(default="surface", metadata={"key": "Style"})
    quality: int = field(default=3, metadata={"key": "Quality"})
    separate: bool = field(default=True, metadata={"key": "Separate"})
    attribute: str = field(default="chain_id", metadata={"key": "Attribute"})
    scale_radii: float = field(default=1.5, metadata={"key": "Scale Radii"})
    probe_size: float = field(default=1.0, metadata={"key": "Probe Size"})
    triangulate: bool = field(default=False, metadata={"key": "Triangulate"})
    relaxation_steps: int = field(default=10, metadata={"key": "Relaxation Steps"})
    by_ca: bool = field(default=False, metadata={"key": "by CA"})
    blur: int = field(default=2, metadata={"key": "Blur"})
    shade_smooth: bool = field(default=True, metadata={"key": "Shade Smooth"})


StyleType = Union[BallStickStyle, CartoonStyle, RibbonStyle, SpheresStyle, SticksStyle, SurfaceStyle]

class StyleCreator():
    def new() -> StyleType:
        return CartoonStyle()

    def ball_stick() -> BallStickStyle:
        return BallStickStyle()

    def cartoon() -> CartoonStyle:
        return CartoonStyle()

    def ribbon() -> RibbonStyle:
        return RibbonStyle()

    def spheres() -> SpheresStyle:
        return SpheresStyle()

    def sticks() -> SticksStyle:
        return SticksStyle()

    def surface() -> SurfaceStyle:
        return SurfaceStyle()
