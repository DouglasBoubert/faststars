"""
"""
from .constants import RSOLAR, DTR
import numpy as np

__all__ = ['rgc_to_dhel']


def rgc_to_dhel(RGC,GLON,GLAT):
    if RGC<RSOLAR:
        raise ValueError('RGC<RSOLAR and so there are two possible distances.')
    DHEL = RSOLAR*np.cos(DTR*GLON)*np.cos(DTR*GLAT)+np.sqrt(RGC**2.-(1.-(np.cos(DTR*GLON)*np.cos(DTR*GLAT))**2.)*RSOLAR**2.)
    return DHEL
