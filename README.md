# IC-FERSR-REACT

## Overview

This repository contains the principle source API code for coupling [IC-FERST](http://multifluids.github.io/) (Imperial College Finite
Element Reservoir SimulaTor), a next generation three-dimensional (3D) reservoir simulator
based on the Double-Control-Volume Finite Element methods and dynamic unstructured mesh
optimization, with [PHREEQC](https://www.usgs.gov/software/phreeqc-version-3), a state-of-the-art geochemical package.

`IC-FERST` is developed by Imperial College London and `PHREEQC` is developed by the United Stated Geological Survey (USGS).

The Imperial College Finite Element Reservoir Simulator (IC-FERST) is available under the terms of the Lesser General-Purpose License (LGPL). IC-FERST has been parallelized using MPI and also has been tested on the U.K. national supercomputer Archer. Moreover, IC-FERST can solve for Darcy as well as Navier-Stokes formulations in two and three dimensions. 

IC-FERST presents a novel high-order, control-volume-finite-element method for the simulation of compositional multiphase flow in porous media with discontinuous first-order representation for pressure and discontinuous second-order representation for velocity. The general advantages of IC-FERST are (i) allowing for both continuous and discontinuous description of pressure and saturation between elements; (ii) the use of arbitrary high-order polynomial representation for pressure and velocity, and (iii) the use of high-order flux limiter in space and time to avoid introducing non-physical oscillations while achieving high-order accuracy where and when possible. The method is implemented using unstructured tetrahedral meshes to discretize space and can easily handle complex geometries [Salinas et al.2017b](https://onlinelibrary.wiley.com/doi/full/10.1002/fld.4381). The implemented numerical method in IC-FERST was presented by [Gomes et al., 2017](https://abdn.pure.elsevier.com/en/publications/a-force-balanced-control-volume-finite-element-method-for-multi-p).

PHREEQC is a general-purpose geochemical reaction model developed by the US Geological Survey. Written in C and C++, PHREEQC is designed to perform a wide range of geochemical calculations and simulations of geochemical reactions and transport processes in both natural or contaminated water and in laboratory and field-scale experiments. This program is based on the equilibrium chemistry of aqueous solutions interacting with minerals, gases, solid solutions, exchangers, and sorption surfaces, which accounts for the original acronym—PH-REdox-Equilibrium (in C language). The program has evolved to include the capability of modeling kinetic reactions and one-dimensional (1D) advection and dispersion transport. Rate equations
are completely user-specifiable in the form of Basic statements. Kinetic and equilibrium reactants can be interconnected. Extensive chemical databases allow the application of reaction, transport, and inverse-modeling capabilities to almost any chemical reaction that is recognized to influence rainwater, soil-water, groundwater, and surface-water quality. PHREEQC can be used to calculate speciation, saturation indices (SI) of minerals, mixing of solutions, reaction path modelling,
inverse modelling and more. Calculations are made using designated thermodynamic databases which include a wide range of data for mineral phases and compounds. PHREEQC allows the concentration of an element to be adjusted to obtain equilibrium (or a specified saturation index or gas partial pressure) with a specified phase, or to obtain charge balance. Solution compositions
can be specified with a variety of concentration units (Parkhurst and Appelo, 1999, 2013).

## Using the API code

The API code could be run with MATLAB after instaling the PHREEQCRM and IC-FERST:

Instaling the PHREEQCRM: https://www.usgs.gov/software/phreeqc-version-3

Instaling the IC-FERST: http://multifluids.github.io/license/
