"""Faststars specific constant variables.
"""
from astropy import constants as const
from astropy import units as un
import numpy as np


CLIGHT = const.c.cgs.value
KM = (1.0 * un.km).cgs.value

MAX_VISUAL_BANDS = [
    ['B', 'b', 'g'],  # B-like bands first
    ['V', 'G'],       # if not, V-like bands
    ['R', 'r']        # if not, R-like bands
]

RSOLAR = 8.29

DTR = np.pi/180.
