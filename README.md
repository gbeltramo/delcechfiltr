# Delaunay-&Ccaron;ech filtrations

This package provides functions for computing the smallest enclosing ball
radius (miniball radius) of sets of three points in two and three dimensions and
four points in three dimensions.

These are used to produce Delaunay-&Ccaron;ech filtrations, which can be used to
compute &Ccaron;ech persistence diagrams of points in two and three dimensions.

Wrappers are provided for computing &Ccaron;ech persistence diagrams with
<a href="http://gudhi.gforge.inria.fr/python/latest/">`gudhi`</a>
for points in two and three dimensional space.

## Installation

The `setuptools`, `numpy`, `scipy`, and `gudhi` Python packages are
prerequisites for using this package.

The last one can be installed with `conda` running

```
>>> conda install -c conda-forge gudhi
```

You'll also need `pybind11`, `cmake` and a `C++` compiler in order to install this
package. It is recommended to install both `pybind11` and `cmake` via `conda`.

```
>>> conda install -c anaconda cmake
>>> conda install -c conda-forge pybind11
```

Finally run the following command to build and install this Python package

```
>>> pip install delcechfiltr
```

If the command above fails to build this package, then you might want to clone
this repository to a directory `delcechfiltr/` on your computer and run

```
>>> cd /<path>/<to>/<cloned>/<repo>/delcechfiltr/
>>> python setup.py install
```

which outputs more information on the errors causing the build to fail.

**Windows.** After installing `conda`, run the above commands within an
`Anaconda prompt`. For the C++ compiler install
<a href="https://visualstudio.microsoft.com/vs/">Visual Studio community</a>.

## Basic usage

```python
import numpy as np
import delcechfiltr.dgms

n, d = 20, 3
np.random.seed(0)
pts = np.random.rand(n, d)

dgm_H0_cech1, dgm_H1_cech1                        = delcechfiltr.dgms.cech(pts, persistence_dim_max=False)
dgm_H0_cech2, dgm_H1_cech2, dgm_H2_cech2          = delcechfiltr.dgms.cech(pts, persistence_dim_max=True)
dgm_H0_delcech1, dgm_H1_delcech1                  = delcechfiltr.dgms.delcech_3D(pts, persistence_dim_max=False)
dgm_H0_delcech2, dgm_H1_delcech2, dgm_H2_delcech2 = delcechfiltr.dgms.delcech_3D(pts, persistence_dim_max=True)
```
