#ifndef DELCECHFILTR_TETRAHEDRON_HPP
#define DELCECHFILTR_TETRAHEDRON_HPP

using namespace std;

namespace delcechfiltr_tetra {

    double circumradius(vector<vector<double>> points);
    double determinant_2x2(vector<vector<double>> matrix);
    double determinant_3x3(vector<vector<double>> matrix);
    vector<vector<double>> invert_3x3(vector<vector<double>> matrix);
    bool is_on_correct_side(vector<double> A, vector<double> B, vector<double> C,
                            vector<double> test_point, vector<double> circum);
    vector<double> circumcenter(vector<vector<double>> points);

    double cech_parameter(vector<vector<double>> points);
    vector<double> cech_param_list_tetrahedra(vector<vector<double>> points,
                                              vector<vector<size_t>> tetra);
}

#endif
