![Language](https://raw.githubusercontent.com/imronuke/KOMODO/master/docs/images/fortran.png) [![Build Status](https://travis-ci.com/imronuke/KOMODO.svg?branch=master)](https://travis-ci.com/imronuke/KOMODO) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/imronuke/KOMODO/blob/master/LICENSE)  [![codecov](https://codecov.io/gh/imronuke/KOMODO/branch/master/graph/badge.svg)](https://codecov.io/gh/imronuke/KOMODO)

# KOMODO
## An Open Nuclear Reactor Simulator

**Documentation available at**: https://imronuke.github.io/KOMODO/

**Features:**
* GPU accelelator (using OpenACC) is partially supported. Useful for a large problem with many nodes.
* Input is straightforward, modular and in a free-format form
* Solves both static and transient core problems **with or without TH feedback**
* Performs forward, adjoint and fixed-source calculations
* Performs calculations using branched cross sections data. An example of the library format can be seen [here](https://github.com/imronuke/KOMODO/blob/master/smpl/xsec/SERPENT_CMM/m40.tab)
* Critical boron concentration search
* Rod ejection simulation or Reactivity Initiated Accident (RIA)
* Solves multi-group of neutron energy
* Solves problems with Assembly Discontinuity Factors (ADFs)
* CMFD accelerated using two-node problem non-linear iteration
* CMFD matrix is solved with the latest linear system solver: BiCGSTAB
* Thermal-hydraulics solutions are obtained by solving mass and energy conservation equations in an enclosed channel
* Three nodal kernels are available:
  * Traditional Finite Difference Method
  * Polynomial Nodal Method (PNM) which is equivalent to Nodal Expansion Method (NEM)
  * Semi-Analytic Nodal Method (SANM)

# KOMODO
KOMODO is an open nuclear reactor simulator that solves both static and transient neutron diffusion equation for one, two or three dimensional reactor problems in Cartesian geometry. Currently, by default, KOMODO uses Semi-Analytic Nodal Method (SANM) to spatially discretise the neutron diffusion equation. While theta method is used for the time discretisation.

KOMODO development was mainly motivated by the cumbersome process to obtain computer codes for most nuclear engineering students. And even so, some of them are not completely free. KOMODO is a great learning tool for reactor theory classes, and we have been striving hard to make the input is easy to create. It is open and free, so everyone has access to the source codes and play with them.

KOMODO is continuation of [ADPRES](https://github.com/imronuke/ADPRES). Since KOMODO name is considered more neutral and institutional-independent, it is expected that more contributors would join in this project.

# User Guides

Here you can find quick and complete guides on how to use KOMODO. Given you have background in nuclear engineering, **we believe you can create your own KOMODO input within minutes!**
## [Theory and Background](https://imronuke.github.io/KOMODO/method)
## [Installation (Building from Source Codes)](https://imronuke.github.io/KOMODO/install)
## [Quick guides](https://imronuke.github.io/KOMODO/quick-guides)
## [Complete guides](https://imronuke.github.io/KOMODO/card-desc)


# How to give feedbacks
You may raise an issue or contact me at
* makrus.imron[at]gmail.com

# How to cite
If you find this work helpful and use this work for a publication, you may cite it as

**Imron, M. (2019). [Development and verification of open reactor simulator ADPRES](https://doi.org/10.1016/j.anucene.2019.06.049). Annals of Nuclear Energy, 133, 580–588.**


> **"The best of people are those who bring most benefit to the rest of mankind." (THE PROPHET)**
