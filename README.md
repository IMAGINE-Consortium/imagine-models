# IMAGINE Model Library

A model library usable for Galactic inference engines. 
The library is (mostly) written C++, but can be accessed via both C++ and Python.
Some (non-essential) elements are written in pure Python, which presently cannot be accessed in the C++ version of the library.
A full list of implemented models can be found [here](#list-of-models).

<span style="color:red;font-weight:700;font-size:18px">
    Please note that this library is work in progress and in early development. 
    By far not everything has been tested, and interfaces may still change.
    Don't hesitate to use this software, but beware of the risks.
</span>

## Requirements

Requires:

- [Python 3](https://www.python.org/) (>=3.6)
- [C++](https://www.python.org/) (>=17)

Python Libraries:

- [NumPy](https://numpy.org/) (>=1.22)

Optional:

- [FFTW3](http://fftw.org/) (>3.3, necessary for random models)
- [autodiff](https://autodiff.github.io/) (>1.0, necessary for building gradients w.r.t. to model parameters)

Optional (Developers):

- [pybind11](https://pybind11.readthedocs.io/en/stable/installing.html)

Note that pybind is included in the current repository, so in the current state you should not need to install it.
This may change in the future.  


## Installation 

The installation procedures are preliminary only 


### Installation (Python)

Easiest via

    pip3 install --user git+https://github.com/IMAGINE-Consortium/imagine-models.git

If you install the library in a virtual environment, remove the `--user` tag.

If you want to specify a branch, you can do so by adding @branch-name to the above command.

Alternatively, (e.g. if you want to add your own model) you need to clone the repository via

    git clone --recursive https://github.com/IMAGINE-Consortium/imagine-models.git

The `--recursive` flag makes sure that also the pybind module is cloned (to the `/extern` folder).

The package can then be installed with

    python3 -m pip install -e folder/where/setup/py/is/

where `-e` flag makes it possible to edit the source files directly, which is convenient for developing. 

The installer will automatically figure out which of the optional dependencies you have installed. 
If you want to disable those, you can do so by defining the `USE_AUTODIFF=OFF` and `USE_FFTW=OFF` environment variables BEFORE you run pip. 


### Installation (C++)

First, one needs download the source files, e.g. via cloning the repository with

    git clone --recursive https://github.com/IMAGINE-Consortium/imagine-models.git

The `--recursive` flag makes sure that also the pybind11 module is cloned (to the `/extern` folder), in case you also want to build the Python package. 
One can then install the C++ library via:
 
```
mkdir build
cd build 
cmake ..
sudo make install 
```

The installer will automatically figure out which of the optional dependencies you have installed. 

If you want to disable those, you can do so by defining the `USE_AUTODIFF=OFF` and `USE_FFTW=OFF` environment variables BEFORE you run cmake. 


## Examples

### Including the library

Example scripts demonstrating how to include both the python and C++ version are located in the ./demo folder. 
In the python case, we also include Jupyter notebooks. In these notebooks you can find an overview with information and plots for each of the models listed below.

### Adding new models **TBD** 

Defining your own models is easy in both Python and C++. 
For that you can start with the model templates provided in the `\templates` folder. 
Binding a C++ model to Python is a bit more involved. For the simplest case, a template exists as well. 


## List of Models

### Magnetic/Vector Fields

| MODEL NAME          | PYTHON      | C++         | reference    | notes               | original implementation |  notebook and plots |
| -----------         | ----------- | ----------- | -----------  | -----------         | -----------             | -----------  |
| **Regular models**  |             |             |              |                     |                         |  |
| Uniform             | &#x2714;    | &#x2714;    |              | used for unit tests |                         |  |
| Helix               | &#x2714;    | &#x2714;    |              |                     |                         |  |
| Axissymetric spiral | &#x2714;    | &#x2718;    | Pelgrims, V. |                     |  Pelgrims, V.           |  |
| Archimedean spiral | &#x2714;    | &#x2714;    |  | simple demonstrative ASS model |   [CRPropa](https://github.com/CRPropa/CRPropa3/tree/master)        |   [ipynb](https://github.com/IMAGINE-Consortium/imagine-models/blob/main/demos/python/model_examples/archimedes_demo.ipynb) |
| Local Bubble | &#x2714;    | &#x2718;    | [Pelgrims et al.](https://www.aanda.org/articles/aa/full_html/2020/04/aa37157-19/aa37157-19.html) |   Only defined on the shell                   |  Pelgrims, V.           |
| Jaffe               | &#x2714;    | &#x2714;    | [Jaffe et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010MNRAS.401.1013J/abstract)        | based on ASS-A spiral with modifications, parameter values taken from hammurabi, not from any publication | [Hammurabi X](https://github.com/hammurabi-dev/hammurabiX)  | [ipynb](https://github.com/IMAGINE-Consortium/imagine-models/blob/main/demos/python/model_examples/jaffe_demo.ipynb) |
| Sun2007              | &#x2714;    | &#x2714;    | [Sun et al. (2007)](https://www.aanda.org/articles/aa/abs/2008/02/aa8671-07/aa8671-07.html) | ASS+Ring as disk field implemented, toroidal asymmetric halo with updated halo parameter from [Sun et al. (2010)](https://iopscience.iop.org/article/10.1088/1674-4527/10/12/009), central part of disk field is constant in z-direction (unphysical) | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/)  | [ipynb](https://github.com/IMAGINE-Consortium/imagine-models/blob/main/demos/python/model_examples/sun_demo.ipynb) |
| Han2018               | &#x2714;    | &#x2714;    | [Han et al. (2018)](https://iopscience.iop.org/article/10.3847/1538-4365/aa9c45)        | BSS-S disk field | / | [ipynb](https://github.com/IMAGINE-Consortium/imagine-models/blob/main/demos/python/model_examples/han_demo.ipynb) |
| Pshirkov | &#x2714;    | &#x2714;    | [Pshirkov et al.](https://iopscience.iop.org/article/10.1088/0004-637X/738/2/192)  | possible to change between ASS-S and BSS-S, and to switch halo field on/off which is asymmetrical with respect to the plane (A)|   [CRPropa](https://github.com/CRPropa/CRPropa3/tree/master)        | [ipynb](https://github.com/IMAGINE-Consortium/imagine-models/blob/main/demos/python/model_examples/pshirkov_demo.ipynb) |
| HMR              | &#x2714;    | &#x2714;    | [Harari et al. (1999)](https://arxiv.org/abs/astro-ph/9906309)        | BSS-S model | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/) and [Kachelrieß (2007)](https://arxiv.org/pdf/astro-ph/0510444.pdf )  | [ipynb](https://github.com/IMAGINE-Consortium/imagine-models/blob/main/demos/python/model_examples/hmr_demo.ipynb) |
| TT              | &#x2714;    | &#x2714;    | [Tinyakov and Tkachev (2017)](https://arxiv.org/abs/astro-ph/0111305)        | BSS-A model implemented (eq. 5 in ref.) | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/) and [Kachelrieß (2007)](https://arxiv.org/pdf/astro-ph/0510444.pdf )  | [ipynb](https://github.com/IMAGINE-Consortium/imagine-models/blob/main/demos/python/model_examples/tt_demo.ipynb) |
| TF              | &#x2714;    | &#x2714;    | [Terral and Ferriere (2017)](https://arxiv.org/abs/1611.10222)        | different halo and disk models available.  Note that only the halo has been fitted using data, leading to very strong magnetic field strengths and unexpected features in the disk fields. Also, the halo fields can converge to infinite field strengths at infinite r/z. Thus, the field should be viewed more as a mathematical exercise than an actual model and should not be used for e.g. cosmic ray propagation. | [CRPropa](https://github.com/CRPropa/CRPropa3/tree/master)        | [ipynb](https://github.com/IMAGINE-Consortium/imagine-models/blob/main/demos/python/model_examples/tf17_demo.ipynb) |
| Fauvet              | &#x2714;    | &#x2714;    | [Fauvet et al. (2012)](https://arxiv.org/abs/1201.5742)        | modified logarithmic spiral with two arms (BSS-S) and z-component, fitted to simulated data (no real data used!) | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/)  | [ipynb](https://github.com/IMAGINE-Consortium/imagine-models/blob/main/demos/python/model_examples/fauvet_demo.ipynb) |
| Stanev            | &#x2714;    | &#x2714;    | [Stanev (1996)](https://arxiv.org/abs/astro-ph/9607086)        | BSS-S model implemented, change in halo field at abs(z)=0.5 was not in hammurabi implementation | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/)  | [ipynb](https://github.com/IMAGINE-Consortium/imagine-models/blob/main/demos/python/model_examples/stanev_demo.ipynb) |
| WMAP              | &#x2714;    | &#x2714;    | [Page et al. (2007)](https://iopscience.iop.org/article/10.1086/513699)  | logarithmic spiral with constant amplitude(B) and z-component, parameters are taken from original publication, not from update mentioned in [Ruiz-Granados et al 2010](https://www.aanda.org/articles/aa/full_html/2010/14/aa12733-09/aa12733-09.html). | / | [ipynb](https://github.com/IMAGINE-Consortium/imagine-models/blob/main/demos/python/model_examples/wmap_demo.ipynb) |
| Jansson Farrar      | &#x2714;    | &#x2714;    | [Jansson & Farrar (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...757...14J/abstract)        | regular JF12 field (disk + symmetric toroidal halo + X-field in z-direction)  | [Hammurabi X](https://github.com/hammurabi-dev/hammurabiX)   | [ipynb](https://github.com/IMAGINE-Consortium/imagine-models/blob/main/demos/python/model_examples/jf12_regular_demo.ipynb) |
| **Random models**   |             |             |              |                     |                         |
| Jansson Farrar      | &#x2714;    |&#x2714;     | [Jansson & Farrar (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...761L..11J/abstract)        | depends on the JF12 regular model   |[Hammurabi X](https://github.com/hammurabi-dev/hammurabiX)    | [ipynb](https://github.com/IMAGINE-Consortium/imagine-models/blob/main/demos/python/model_examples/jf12_random_demo.ipynb) |
| Ensslin Steininger  | &#x2714;    |&#x2714;     |              |                      | [Hammurabi X](https://github.com/hammurabi-dev/hammurabiX)             |  |

### Thermal electron/Scalar Fields

| MODEL NAME  | PYTHON      | C++         | reference   | notes       |  original implementation | notebook and plots |
| ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ---------- |
| **Regular models**        |             |             |             |             |            |  |
| YMW16       | &#x2714;    | &#x2714;    | [Yao et al. (2016)](https://ui.adsabs.harvard.edu/abs/2017ApJ...835...29Y/abstract)     |             |       [Hammurabi X](https://github.com/hammurabi-dev/hammurabiX)         |  |
| **Random models** |             |             |             |             |             |
| GaussianScalar | &#x2714;|&#x2714; |             |      used for random number unit testing (TBD!)       |             | |
| LogNormalScalar | &#x2714;|&#x2714; |             |             |             |
