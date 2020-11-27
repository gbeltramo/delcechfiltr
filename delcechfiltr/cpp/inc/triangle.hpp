#ifndef DELCECHFILTR_TRIANGLE_HPP
#define DELCECHFILTR_TRIANGLE_HPP

using namespace std;

namespace delcechfiltr_tri {

    double circumradius_2D(vector<vector<double>> points);
    double circumradius_3D(vector<vector<double>> points);
    double cech_parameter_2D(vector<vector<double>> points);
    double cech_parameter_3D(vector<vector<double>> points);
    vector<double> circumcenter_2D(vector<vector<double>> points);
    vector<double> circumcenter_3D(vector<vector<double>> points);
    vector<double>  miniball_center_2D(vector<vector<double>> points);
    vector<double>  miniball_center_3D(vector<vector<double>> points);
    vector<double> cech_param_list_triangles_2D(vector<vector<double>> points,
                                                vector<vector<size_t>> triangles);
    vector<double> cech_param_list_triangles_3D(vector<vector<double>> points,
                                                vector<vector<size_t>> triangles);
}

#endif
