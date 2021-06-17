===
IBS
===

Installation
============

Depends on:

- cmake
- make
- g++

.. code-block:: bash

    $ git clone https://github.com/tomerten/IBSLib.git
    $ cd IBSLib
    $ cmake .
    $ make
    $ ./IBSLib

Current Supported Models
========================

- Piwinski smooth lattice approximation
- Piwinski Lattice element by element weighted
- Piwinski Lattice Modified taking some vertical effects into account
- Nagaitsev's high-energy approximation 
- Bjorken-Mtingwa
- Conte-Martini
- Madx (CERN note AB-2006-002) using `TWINT` and `SIMPSONDECADE` methods to perform the integration.

Coublomb Log methods
====================

- twclog - uses element by element twiss data
- twclogtail - uses element by element twiss data
- CoublombLog - uses ring averages 
- TailCutCoulombLog - uses ring averages


Integration methods
===================

- Simpson (standard implementation)
- SimpsonDecade - Simspon per decade for covering large spread in integration ranges (ususally 50 orders of magnitude difference between low and high)
- TWINT, based on GSL QUAD method

Radiation Damping
=================

- Radiation Damping using smooth lattice approximation 
- Radiation Damping element by element
- Equilibrium from pure radiation damping and exitation (taux, tauy, taus, exinf, eyinf, sigeoe2, jx, jy)
- Radiation losses per turn
- Critical omega, theta, photon energy

Numeric Functions
=================

- sigefromsigs
- eta 
- fmohl
- particle_radius
- BetaRelativisticFromGamma
- rds 
- VeffRFeV
- VeffRFeVPrime
- synchronuousphase
- VeffRFeVPotentialWellDistortion
- VeffRFeVPotentialWellDistortionPrime
- synchronuousphasewithPWD
- synchrotronTune
- synchrotronTunewithPWD
- csige (calculates sige from RF settings, radiation losses and sigs)

ODE 
===

The tables below show examples of the output of the ODE method, with `maxsteps` 
set to two in order to limit the output,
written to a csv file.

.. csv-table:: ODE Model using Piwinski Smooth
    :file: ODE_test_output_piwinski_smooth.csv 
    :header-rows: 1

|

.. csv-table:: ODE Model using Piwinski Lattice
    :file: ODE_test_output_piwinski_lattice.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Piwinski Lattice Modified   
    :file: ODE_test_output_piwinski_latticemodified.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Nagaitsev   
    :file: ODE_test_output_nagaitsev.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Nagaitsev Tailcut
    :file: ODE_test_output_nagaitsevtailcut.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using MADX (Zimmerman)
    :file: ODE_test_output_madx.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using MADX (Zimmerman) with Tailcut
    :file: ODE_test_output_madxtailcut.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Bjorken-Mtingwa with standard Simpson integration (Fails for ey)
    :file: ODE_test_output_bjorken_mtingwa2.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Bjorken-Mtingwa with Simpson Decade Integration 
    :file: ODE_test_output_bjorken_mtingwa.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Bjorken-Mtingwa with Simpson Decade Integration and Tailcut
    :file: ODE_test_output_bjorken_mtingwatailcut.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Conte-Martini using Simspon Decade Integration
    :file: ODE_test_output_conte_martini.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Conte-Martini using Simspon Decade Integration and Tailcut
    :file: ODE_test_output_conte_martini_tailcut.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using MADX (Zimmerman) using Simpson Decade Integration 
    :file: ODE_test_output_madxibs.csv 
    :header-rows: 1