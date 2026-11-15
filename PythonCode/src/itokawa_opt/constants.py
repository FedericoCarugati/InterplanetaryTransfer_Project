"""Physical constants and unit conventions.

Units throughout the package: km, s, km/s, rad (same as the MATLAB code).

The values marked "thesis" are the ones used in the BSc thesis MATLAB scripts and
are kept identical so that the Python validation reproduces the MATLAB numbers.
"""

# Sun gravitational parameter used in the thesis MATLAB code [km^3/s^2]
MU_SUN_THESIS: float = 1.327e11
# Reference (IAU / DE440) value, to be used from Level 1 onwards [km^3/s^2]
MU_SUN: float = 1.32712440018e11

# Astronomical unit [km]
AU: float = 149_597_870.7

# Universal gravitational constant used in the thesis [m^3 kg^-1 s^-2]
G_THESIS: float = 6.67e-11
# Solar mass used in the thesis [kg]
M_SUN_THESIS: float = 1.989e30
# Itokawa mass used in the thesis [kg]
M_ITOKAWA: float = 3.1464e10

DAY: float = 86_400.0
