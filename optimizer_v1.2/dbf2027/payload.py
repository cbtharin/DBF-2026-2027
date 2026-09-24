"""Payload geometry: sensor envelope, shipping container, and bay packing.

The shipping container is a SINGLE FIXED ASSUMPTION - one box, the same for
every design - set by `Rules.container_*` in config.py.  Its size comes from
the drop-test calculation for a roughly 1 kg sensor: a 5 ft drop reaches
5.47 m/s, holding the sensor under 75 g needs a 20 mm crush stroke, and
closed-cell foam crushes usefully through about 70 percent of its thickness,
so 29 mm of foam per face around a 165 x 83 x 83 mm sensor.  Rounded to a
buildable 230 x 150 x 150 mm box at 0.5 lb tare.

What still varies per design is how MANY boxes there are, and therefore how
big the payload bay is - which is the part that actually drives fuselage
size and drag.  Simulators share the container envelope, because that is
what they simulate.

The sensor's own envelope is still estimated at an assumed packing density,
but only so `PayloadGeometry.fits()` can reject a sensor too big for the box.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple

from .config import Rules


@dataclass
class PayloadGeometry:
    sensor_l: float
    sensor_w: float
    sensor_h: float
    pad_m: float                 # foam thickness per face
    box_l: float                 # container external dimensions
    box_w: float
    box_h: float
    container_tare_kg: float

    @property
    def box_volume_m3(self) -> float:
        return self.box_l * self.box_w * self.box_h

    @property
    def box_frontal_m2(self) -> float:
        return self.box_w * self.box_h

    @property
    def sensor_volume_m3(self) -> float:
        return self.sensor_l * self.sensor_w * self.sensor_h

    def fits(self, fill_fraction: float) -> bool:
        """Would the sensor fit at its NOMINAL packing density?

        Informational only.  A sensor heavier than this just has to be
        denser - see `required_density`.  Nothing rejects a design for
        failing this.
        """
        return self.sensor_volume_m3 <= fill_fraction * self.box_volume_m3

    def required_density(self, sensor_mass_kg: float,
                         fill_fraction: float) -> float:
        """Density the sensor must actually reach to fit the fixed box."""
        usable = max(fill_fraction * self.box_volume_m3, 1e-12)
        return sensor_mass_kg / usable


def drop_impact_speed_mps(height_in: float) -> float:
    """Speed at the end of a free fall from `height_in` inches.

    v = sqrt(2 g h).  The Ground Mission drops the container onto concrete
    from a height the team declares, up to 60 inches.
    """
    return math.sqrt(2.0 * 9.80665 * max(height_in, 0.0) * 0.0254)


def crush_stroke_m(height_in: float, survivable_g: float = 75.0) -> float:
    """Foam stroke needed to hold the sensor under `survivable_g`.

    Constant-force crush: s = v^2 / (2 a), with a = g * survivable_g.
    Real foam crushes usefully through about 70 percent of its thickness,
    so the pad is thicker than the stroke.
    """
    v = drop_impact_speed_mps(height_in)
    return v * v / (2.0 * 9.80665 * max(survivable_g, 1.0))


# RULE 3.1.3b, 3.1.3c - one rectangular container envelope, used for
# the real box and for every simulator, which must match it in size.
def size_payload(sensor_mass_kg: float, rules: Rules,
                 drop_height_in: Optional[float] = None) -> PayloadGeometry:
    """One fixed container for every design.

    The box is NOT resized per sensor: the external dimensions and tare in
    `Rules` are a single assumption.  What DOES respond to the design is
    the foam thickness, because the Ground Mission lets the team declare
    its own drop height up to 60 inches and the score rises with it - so
    the pad has to be sized for the height actually declared, not for a
    fixed 5 ft.

    The sensor envelope is still reported so you can see whether it fits,
    and `fits()` is what the feasibility check uses.
    """
    vol = max(sensor_mass_kg, 1e-9) / rules.sensor_packing_density_kg_m3
    w = (0.5 * vol) ** (1.0 / 3.0)          # 2:1:1 sensor envelope
    l, h = 2.0 * w, w
    h_in = (rules.max_drop_height_in if drop_height_in is None
            else max(0.0, min(drop_height_in, rules.max_drop_height_in)))
    pad = crush_stroke_m(h_in, rules.sensor_survivable_g) / 0.70
    return PayloadGeometry(l, w, h, pad,
                           rules.container_length_m,
                           rules.container_width_m,
                           rules.container_height_m,
                           rules.container_tare_kg)


def best_bay_arrangement(n_boxes: int, box_l: float, box_w: float,
                         box_h: float) -> Tuple[int, int, int, float, float, float]:
    """Pack n identical boxes into the smallest-wetted-area cuboid bay.

    Returns (rows_along, cols_across, layers_high, bay_l, bay_w, bay_h).

    One long row gives a long thin fuselage, side-by-side gives a short fat
    one; which is cheaper depends on the box proportions, so the arrangement
    is chosen rather than assumed.
    """
    if n_boxes <= 0:
        return 0, 0, 0, 0.0, 0.0, 0.0
    best = None
    for cols in range(1, min(n_boxes, 4) + 1):
        for layers in range(1, min(n_boxes, 3) + 1):
            rows = int(math.ceil(n_boxes / (cols * layers)))
            if rows * cols * layers < n_boxes:
                continue
            bl = rows * box_l
            bw = cols * box_w
            bh = layers * box_h
            # proxy for fuselage cost: wetted area of the bay itself
            wetted = 2.0 * (bl * bw + bl * bh + bw * bh)
            # a bay wider or taller than it is long makes an unflyable fuselage
            if bw > 0.55 * bl or bh > 0.55 * bl:
                wetted *= 1.6
            if best is None or wetted < best[0]:
                best = (wetted, rows, cols, layers, bl, bw, bh)
    _w, rows, cols, layers, bl, bw, bh = best
    return rows, cols, layers, bl, bw, bh
