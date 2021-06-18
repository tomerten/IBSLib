#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tests for C++ module twiss.
"""

import ibslib_pb as ibslib
import pandas as pd
import pytest

constants = [
    (ibslib.clight, 299792458.0),
    (ibslib.hbarGeV, 6.582119569e-25),
    (ibslib.electron_mass, 0.51099895000e-3),
    (ibslib.proton_mass, 0.93827208816),
    (ibslib.neutron_mass, 0.93956542052),
    (ibslib.mu_mass, 0.1056583755),
    (ibslib.atomic_mass_unit, 0.93149410242),
    (ibslib.pi, 3.141592653589793),
    (ibslib.electric_charge, 1.602176634e-19),
    (ibslib.euler, 0.577215664901533),
    (ibslib.electron_radius, 2.8179403262e-15),
    (ibslib.proton_radius, 1.5346982671888944e-18),
]


@pytest.mark.parametrize("name, value", constants)
def test_constants(name, value):
    assert name == value


def test_cpp_sigefromsigs():
    print(ibslib.sige_from_sigs(ibslib.pi * 2 * 1.25e6, 0.005, 5e-4, 3326.0, 37.0))
    assert (ibslib.sige_from_sigs(ibslib.pi * 2 * 1.25e6, 0.005, 5e-4, 3326.0, 37.0)) < 1e-2


def test_cpp_sigsfromsige():
    val = ibslib.sigs_from_sige(8.96628617341675e-05, 3326.0, 37.0, 5e-4 * ibslib.pi * 2 * 1.25e6)
    print(val)
    assert (val < 0.0051) & (val > 0.004999)
