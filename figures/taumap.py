# -*- coding: utf-8 -*-
"""The diverging colour ramp of the tau* heat maps (olive = tau* > 0, indigo = tau* < 0, neutral grey at 0), as a small module.

 The poles are interpolated in OKLab from the neutral grey: lightness linear, chroma t**gamma. fig1.py builds the same ramp inline.
 wcag() is the WCAG 2 contrast ratio used by the figure scripts to check that numbers and dots remain legible on the cells.
"""
import math

import numpy as np
from matplotlib.colors import ListedColormap, to_rgb

POS_HEX, NEG_HEX, MID_HEX = "#757b10", "#797ec7", "#efefef"   # olive (OKLCH .56/.12/113), indigo (.62/.11/280), neutral grey


def _s2l(x): return x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4
def _l2s(x): x = min(max(x, 0.0), 1.0); return 12.92 * x if x <= 0.0031308 else 1.055 * x ** (1 / 2.4) - 0.055
def oklab(hexc):
    r, g, b = (_s2l(v) for v in to_rgb(hexc))
    l, m, s = (np.cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b), np.cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b),
               np.cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b))
    return np.array([0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s, 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
                     0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s])
def oklab_to_srgb(L, a, b):
    l_, m_, s_ = L + 0.3963377774 * a + 0.2158037573 * b, L - 0.1055613458 * a - 0.0638541728 * b, L - 0.0894841775 * a - 1.2914855480 * b
    l, m, s = l_ ** 3, m_ ** 3, s_ ** 3
    lin = [4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s, -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
           -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s]
    return [_l2s(x) for x in lin], max(max(0 - x, x - 1, 0) for x in lin)
def arm(pole, n=128, gamma=1.0):
    """Neutral grey -> pole at the pole's hue: lightness linear (monotone), chroma t**gamma (gamma < 1 keeps light tints chromatic)."""
    L0 = oklab(MID_HEX)[0]; Lp, ap, bp = oklab(pole); Cp, h = math.hypot(ap, bp), math.atan2(bp, ap); out = []
    for t in np.linspace(0, 1, n):
        rgb, oog = oklab_to_srgb(L0 + (Lp - L0) * t, Cp * t ** gamma * math.cos(h), Cp * t ** gamma * math.sin(h))
        assert oog < 2e-3, (pole, t, oog); out.append(rgb)
    return out
TAUMAP = ListedColormap(arm(NEG_HEX, gamma=1.2)[::-1] + arm(POS_HEX, gamma=0.6), name="tau_olive_indigo")


def wcag(c1, c2):
    """WCAG 2 contrast ratio of two colours."""
    def Y(cc): r, g, b = (_s2l(v) for v in to_rgb(cc)); return 0.2126 * r + 0.7152 * g + 0.0722 * b
    a, b = sorted([Y(c1), Y(c2)], reverse=True); return (a + 0.05) / (b + 0.05)
