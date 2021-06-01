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

Integration methods
===================

- Simpson (standard) implementation)
- Simspon per decade for covering large spread in integration ranges (ususally 50 orders of magnitude difference between low and high)
- TWINT, based on GSL QUAD method

Radiation Damping
=================

- Radiation Damping using smooth lattice approximation 
- Radiation Damping element by element
- Equilibrium from pure radiation damping and exitation (taux, tauy, taus, exinf, eyinf, sigeoe2, jx, jy)
- Radiation losses per turn