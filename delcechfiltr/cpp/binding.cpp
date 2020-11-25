#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "inc/triangle.hpp"

using namespace std;

PYBIND11_MODULE(binding, m) {
    m.doc() = "";
    m.def("cross_3D", &delcechfiltr_tri::cross_3D);
    m.def("triangle_circumradius_2D", &delcechfiltr_tri::circumradius_2D);
    m.def("triangle_circumradius_3D", &delcechfiltr_tri::circumradius_3D);
}
