"""Level 1 - SPICE ephemerides (DE440 + Itokawa SPK from JPL Horizons) with the same
``rv(t)`` interface as :class:`itokawa_opt.ephemeris.thesis_kepler.KeplerEphemeris`.

TODO: implement with ``spiceypy``; kernels are downloaded by ``scripts/fetch_kernels.py``
into ``data/kernels/`` (not versioned).
"""
