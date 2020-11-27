#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "inc/triangle.hpp"
#include "inc/tetrahedron.hpp"

using namespace std;

PYBIND11_MODULE(binding, m) {
    m.doc() = "";

    m.def("triangle_circumradius_2D", &delcechfiltr_tri::circumradius_2D);
    m.def("triangle_circumradius_3D", &delcechfiltr_tri::circumradius_3D);
    m.def("triangle_cech_parameter_2D", &delcechfiltr_tri::cech_parameter_2D);
    m.def("triangle_cech_parameter_3D", &delcechfiltr_tri::cech_parameter_3D);
    m.def("triangle_cech_param_list_2D", &delcechfiltr_tri::cech_param_list_triangles_2D);
    m.def("triangle_cech_param_list_3D", &delcechfiltr_tri::cech_param_list_triangles_3D);
    m.def("triangle_circumcenter_2D", &delcechfiltr_tri::circumcenter_2D);
    m.def("triangle_circumcenter_3D", &delcechfiltr_tri::circumcenter_3D);
    m.def("triangle_miniball_center_2D", &delcechfiltr_tri::miniball_center_2D);
    m.def("triangle_miniball_center_3D", &delcechfiltr_tri::miniball_center_3D);

    m.def("tetrahedron_circumradius", &delcechfiltr_tetra::circumradius);
    m.def("tetrahedron_det_2x2", &delcechfiltr_tetra::determinant_2x2);
    m.def("tetrahedron_det_3x3", &delcechfiltr_tetra::determinant_3x3);
    m.def("tetrahedron_inv_3x3", &delcechfiltr_tetra::invert_3x3);
    m.def("tetrahedron_is_on_correct_side", &delcechfiltr_tetra::is_on_correct_side);
    m.def("tetrahedron_circumcenter", &delcechfiltr_tetra::circumcenter);
    m.def("tetrahedron_cech_parameter", &delcechfiltr_tetra::cech_parameter);
    m.def("tetrahedron_cech_param_list", &delcechfiltr_tetra::cech_param_list_tetrahedra);
}
