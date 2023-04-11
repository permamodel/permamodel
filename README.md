[![Test](https://github.com/permamodel/permamodel/actions/workflows/test.yml/badge.svg)](https://github.com/permamodel/permamodel/actions/workflows/test.yml)
[![Coverage Status](https://coveralls.io/repos/github/permamodel/permamodel/badge.svg?branch=main)](https://coveralls.io/github/permamodel/permamodel?branch=main)

Permamodel
==========

Permamodel is a collection of numerical permafrost models with a range of capability and complexity.
Permamodel includes multiple sets of sample inputs representing a variety of conditions and locations.
The Permamodel project is intended to facilitate the broader use of permafrost models.
We hope that the simple Python interfaces and open source licensing can make permafrost models accessible to a broad audience well beyond the permafrost research community, such as educators, students, and policy makers. 

Frost number
-----------
Calculates the *frost number*, an indication of the probability
of finding permafrost, at a site.  There are different versions of the
frost number depending on what data area available at the site.  

* The air frost number requires some indication of the annual temperature cycle.
* The surface frost number also requires indication of the year's precipitation.
* The Stefan frost number additionally incorporates soil information.

Ku
--

Implements a semi-empirical/analytical solution to soil conditions. Please cite:

> Wang, K., Jafarov, E., & Overeem, I. (2020). Sensitivity evaluation of the Kudryavtsev permafrost model. Science of The Total Environment, 137538. https://doi.org/10.1016/j.scitotenv.2020.137538

> Overeem, I., E. Jafarov, K. Wang, K. Schaefer, S. Stewart, G. Clow, M. Piper, and Y. Elshorbany (2018). A modeling toolbox for permafrost landscapes. Eos. https://doi.org/10.1029/2018EO105155.

GIPL
----

GIPL is a numerical model that solves for the temperature profile
of a soil column given its material properties and the temperature and
precipitation conditions it experiences. For more information please see https://github.com/Elchin/GIPL and https://github.com/permamodel/GIPL-BMI-Fortran.

Installation
------------

Permamodel can be installed with `pip`:
```
$ pip install permamodel
```
or with `conda`:
```
$ conda install -c conda-forge permamodel
```
We recommend installing permaodel into Python virtual environment.
