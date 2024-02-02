---
title: The IMAGINE Model library
tags:
- C++
- Python
- astronomy
- Galactic magnetic field
- Milky Way
- non-thermal radiation 
- CMB foreground
- Cosmic rays
authors:
- name: Sebastian Hutschenreuter
  orcid: 0000-0002-6952-9688
  affiliation: "1, 2"
- name: Teresa Bister
  orcid: 0000-0003-4005-0857
  affiliation: "3, 4"
- name: Ashley Stock 
  orcid: 0000-0001-5351-824X
  affiliation: "5, 6, 7"
- name: Vasundhara Shaw
  orcid: 0000-0002-5824-7191
  affiliation: "8"


affiliations:
- name: Institute for Astronomy (IfA), University of Vienna, T\"urkenschanzstrasse 17, A-1180 Vienna
  index: 1
- name: Department of Astrophysics/IMAPP, Radboud University, P.O. Box 9010,6500 GL Nijmegen, The Netherlands
  index: 2
- name: Department of High Energy Physics/IMAPP, Radboud University, P.O. Box 9010,6500 GL Nijmegen, The Netherlands
  index: 3
- name: Nationaal Instituut voor Kernfysica en Hoge Energie Fysica (NIKHEF), SciencePark, Amsterdam, The Netherlands
  index: 4
- name: David A. Dunlap Department of Astronomy and Astrophysics, University of Toronto, 50 St. George Street, Toronto, ON M5S 3H4, Canada  
  index: 5
- name: Dunlap Institute for Astronomy and Astrophysics, University of Toronto, 50 St. George Street, Toronto, ON M5S 3H4, Canada
  index: 6
- name: Canadian Institute for Theoretical Astrophysics, University of Toronto, 60 St. George Street, Toronto, ON M5S 3H8, Canada
  index: 7
- name: Jodrell Bank Centre for Astrophysics, Alan Turing Building, Department of Physics and Astronomy, School of Natural Sciences, Oxford Road, Manchester, M13 9PL, U.K.
  index: 8
date: x x 2023
bibliography: paper.bib
---

# Statement of need

The Galactic magnetic field (GMF) is an integral component of the Milky Way, contributing a significant share to the energy budget of the interstellar medium (ISM) and influencing the dynamics of almost all of its components such as gas, dust, and cosmic rays. [@2013Beck]
The interactions of particles with the GMF within the ISM also produce a multitude of observables, ranging from the radio regime to highly energetic cosmic rays. 
Additionally, possible GMF interactions with e.g. axions make the GMF an ideal laboratory for the search for beyond-Standard-Model particles.
While all of these signals provide us with a wealth of information on the Galactic ISM, they also obscure the universe beyond. This makes it necessary to remove the Galactic contribution in the study of most extragalactic objects, which requires a detailed understanding of Galactic components such as the GMF.

Constraining the GMF, however, is a delicate task as all available tracers are convolved with other components of the Milky Way, e.g. with the thermal electrons in the case of Faraday rotation, relativistic electrons in case of synchrotron radiation or dust properties in the case of (optical and infrared) stellar polarization.

Over the decades, the community has produced several models of the large scale GMF constrained by different tracers. 
Additionally, more local measurements have allowed to constrain the GMF close to the sun's location, specifically in the local bubble. 
These models vary considerably between each other, highlighting the need of a holistic effort to constrain the Galactic magnetic fields. 
[@Boulanger2018] have proposed the IMAGINE framework, which aims to connect this multitude of models and data with modern statistical inference engines.

The aim of this model library is to provide a unified infrastructure to access most available GMF models, not only within the context of the IMAGINE framework, but also as a standalone package with a simple, user-friendly interface that allows the user to include this package easily within their framework. 
Explicitly, we have implemented bridges to the simulation engines Hammurabi [@Wang2020] and CRPRopa [@2022Batista], as well as the Bayesian inference software Nifty8 [@Nifty2019].  *This sentence might have to be revised depending on what is actually implemented when submitting*

While the focus of the library is on parametric GMF models, it also provides infrastructure for scalar fields and specific implementations of thermal electron density models. 
A further extension to other Galactic ISM components such as e.g. relativistic electrons is envisioned.


# Features of the library

The library provides a simple interface to access the magnetic field vector for a particular model at each position in Galacto-centric Cartesian coordinates. 
It furthermore provides convenience functionality to evaluate the models on regular and irregular grids, as well as point clouds.  
Optionally (and requiring the installation of additional packages), the library can also compute turbulent magnetic fields based on the FFTW algorithm [@FFTW2005] and calculate the derivative with respect to model parameters using the [@autodiff2018] package.
The full list of models is available [here](https://github.com/IMAGINE-Consortium/imagine-models#list-of-models). 
At the time of writing this document, 18 magnetic field models and three scalar field models were implemented. 
These models encompass a wide range of complexity from trivial models such as a uniform magnetic field to models depending on dozens of parameters such as the widely used Jansson & Farrar model [@JF2012].  

The library is written in C++ and interfaced to Python using [@pybind112017]. 
A small subset of models is directly implemented on the python side, and hence not available in C++.  
While we aim for almost complete implementation on the C++ in the basic functionality of the library, the fact the models can be implemented in pure python is aimed at users who want to implement their own models within this framework, without having to reside to the somewhat opaque issue of binding C++ code to python.


# Acknowledgements

We acknowledge ... 

# References
