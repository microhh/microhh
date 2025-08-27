import numpy as np
from constants import *


def calc_esat(T):
    """
    Calculate saturation vapor pressure (Pa).
    """
    a_ab = 611.21; b_ab = 18.678; c_ab = 234.5; d_ab = 257.14
    return a_ab * np.exp((b_ab - ((T-T0) / c_ab)) * ((T-T0) / (d_ab + (T-T0))))


def calc_qsat(T, p):
    """
    Calculate saturation specific humidity (kg kg-1)
    """
    esat = calc_esat(T)
    return ep*esat / (p-(1.-ep)*esat)


class Emissions:
    """
    Help class to gather emission settings.
    """
    def __init__(self):

        self.source_list = []
        self.strength = []
        self.sw_vmr = []

        self.x0 = []
        self.y0 = []
        self.z0 = []

        self.sigma_x = []
        self.sigma_y = []
        self.sigma_z = []

        self.line_x = []
        self.line_y = []
        self.line_z = []

    def add(
            self, name, strength, sw_vmr,
            x0, y0, z0,
            sigma_x, sigma_y, sigma_z,
            line_x=0, line_y=0, line_z=0):

        self.source_list.append(name)
        self.strength.append(strength)
        self.sw_vmr.append(sw_vmr)

        self.x0.append(x0)
        self.y0.append(y0)
        self.z0.append(z0)

        self.sigma_x.append(sigma_x)
        self.sigma_y.append(sigma_y)
        self.sigma_z.append(sigma_z)

        self.line_x.append(line_x)
        self.line_y.append(line_y)
        self.line_z.append(line_z)
