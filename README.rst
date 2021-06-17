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

.. csv-table:: ODE Model using Piwinski Smooth 
    :file: ODE_Output_test.csv 
    :header-rows: 1