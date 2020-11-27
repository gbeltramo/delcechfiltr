#include <vector>
#include <cmath>       // std::sqrt
#include <algorithm>   // std::sort
#include "../inc/triangle.hpp"   // delcechfiltr_tri::cech_parameter_3D

using namespace std;

namespace delcechfiltr_tetra {

    double dot(vector<double> A, vector<double> B) {
        double tot = 0.0;
        for (size_t ind = 0; ind < A.size(); ++ind) {
            tot = tot + A[ind] * B[ind];
        }
        return tot;
    }

    double magnitude(vector<double> A) {
        double tot = dot(A, A);
        return sqrt(tot);
    }

    vector<double> cross_3D(vector<double> A, vector<double> B) {
        vector<double> cross_out = {0.0, 0.0, 0.0};
        cross_out[0] = A[1] * B[2] - A[2] * B[1];
        cross_out[1] = A[2] * B[0] - A[0] * B[2];
        cross_out[2] = A[0] * B[1] - A[1] * B[0];
        return cross_out;
    }

    double volume_tetrahedron(vector<vector<double>> points) {
        vector<double> A = {points[0][0] - points[3][0],
                            points[0][1] - points[3][1],
                            points[0][2] - points[3][2]};
        vector<double> B = {points[1][0] - points[3][0],
                            points[1][1] - points[3][1],
                            points[1][2] - points[3][2]};
        vector<double> C = {points[2][0] - points[3][0],
                            points[2][1] - points[3][1],
                            points[2][2] - points[3][2]};
        vector<double> cross_BC = cross_3D(B, C);
        double vol = fabs(dot(A, cross_BC)) / 6.0;
        return vol;
    }

    double circumradius(vector<vector<double>> points) {
        vector<double> A_B = {points[0][0] - points[1][0],
                              points[0][1] - points[1][1],
                              points[0][2] - points[1][2]};
        vector<double> A_C = {points[0][0] - points[2][0],
                              points[0][1] - points[2][1],
                              points[0][2] - points[2][2]};
        vector<double> A_D = {points[0][0] - points[3][0],
                              points[0][1] - points[3][1],
                              points[0][2] - points[3][2]};
        vector<double> C_D = {points[2][0] - points[3][0],
                              points[2][1] - points[3][1],
                              points[2][2] - points[3][2]};
        vector<double> B_D = {points[1][0] - points[3][0],
                              points[1][1] - points[3][1],
                              points[1][2] - points[3][2]};
        vector<double> B_C = {points[1][0] - points[2][0],
                              points[1][1] - points[2][1],
                              points[1][2] - points[2][2]};

        double len_edge_AB = magnitude(A_B);
        double len_opposite_AB = magnitude(C_D);
        double prod_AB = len_edge_AB * len_opposite_AB;
        double len_edge_AC = magnitude(A_C);
        double len_opposite_AC = magnitude(B_D);
        double prod_AC = len_edge_AC * len_opposite_AC;
        double len_edge_AD = magnitude(A_D);
        double len_opposite_AD = magnitude(B_C);
        double prod_AD = len_edge_AD * len_opposite_AD;

        double numer = sqrt((prod_AB + prod_AC + prod_AD) * (prod_AB + prod_AC - prod_AD) * (prod_AB - prod_AC + prod_AD) * (-prod_AB + prod_AC + prod_AD));
        double denom = 24 * volume_tetrahedron(points);

        return numer / denom;
    }

    double determinant_2x2(vector<vector<double>> matrix) {
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
    }

    double determinant_3x3(vector<vector<double>> matrix) {
        double coeff1 = matrix[0][0];
        vector<vector<double>> matrix1 = { {matrix[1][1], matrix[1][2]},
                                           {matrix[2][1], matrix[2][2]} };
        double coeff2 = matrix[0][1];
        vector<vector<double>> matrix2 = { {matrix[1][0], matrix[1][2]},
                                           {matrix[2][0], matrix[2][2]} };
        double coeff3 = matrix[0][2];
        vector<vector<double>> matrix3 = { {matrix[1][0], matrix[1][1]},
                                           {matrix[2][0], matrix[2][1]} };

        double det = coeff1 * determinant_2x2(matrix1) - coeff2 * determinant_2x2(matrix2) + coeff3 * determinant_2x2(matrix3);
        return det;
    }

    vector<vector<double>> invert_3x3(vector<vector<double>> matrix) {
        double det = determinant_3x3(matrix);
        vector<vector<double>> inverse = { {0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0} };

        vector<vector<double>> matrix1 = { {matrix[1][1], matrix[1][2]},
                                           {matrix[2][1], matrix[2][2]} };
        vector<vector<double>> matrix2 = { {matrix[1][0], matrix[1][2]},
                                           {matrix[2][0], matrix[2][2]} };
        vector<vector<double>> matrix3 = { {matrix[1][0], matrix[1][1]},
                                           {matrix[2][0], matrix[2][1]} };

        vector<vector<double>> matrix4 = { {matrix[0][1], matrix[0][2]},
                                          {matrix[2][1], matrix[2][2]} };
        vector<vector<double>> matrix5 = { {matrix[0][0], matrix[0][2]},
                                          {matrix[2][0], matrix[2][2]} };
        vector<vector<double>> matrix6 = { {matrix[0][0], matrix[0][1]},
                                          {matrix[2][0], matrix[2][1]} };

        vector<vector<double>> matrix7 = { {matrix[0][1], matrix[0][2]},
                                           {matrix[1][1], matrix[1][2]} };
        vector<vector<double>> matrix8 = { {matrix[0][0], matrix[0][2]},
                                           {matrix[1][0], matrix[1][2]} };
        vector<vector<double>> matrix9 = { {matrix[0][0], matrix[0][1]},
                                           {matrix[1][0], matrix[1][1]} };

        inverse[0][0] = determinant_2x2(matrix1) / det;
        inverse[0][1] = - determinant_2x2(matrix4) / det;
        inverse[0][2] = determinant_2x2(matrix7) / det;
        inverse[1][0] = - determinant_2x2(matrix2) / det;
        inverse[1][1] = determinant_2x2(matrix5) / det;
        inverse[1][2] = - determinant_2x2(matrix8) / det;
        inverse[2][0] = determinant_2x2(matrix3) / det;
        inverse[2][1] = - determinant_2x2(matrix6) / det;
        inverse[2][2] = determinant_2x2(matrix9) / det;

        return inverse;
    }

    vector<double> circumcenter(vector<vector<double>> points) {
        vector<double> A_D = {points[0][0] - points[3][0],
                              points[0][1] - points[3][1],
                              points[0][2] - points[3][2]};
        vector<double> B_D = {points[1][0] - points[3][0],
                              points[1][1] - points[3][1],
                              points[1][2] - points[3][2]};
        vector<double> C_D = {points[2][0] - points[3][0],
                              points[2][1] - points[3][1],
                              points[2][2] - points[3][2]};
        vector<vector<double>> matrix = { A_D, B_D, C_D };
        vector<vector<double>> inverse = invert_3x3(matrix);
        vector<double> arr = { (dot(points[0], points[0]) - dot(points[3], points[3])) / 2.0,
                               (dot(points[1], points[1]) - dot(points[3], points[3])) / 2.0,
                               (dot(points[2], points[2]) - dot(points[3], points[3])) / 2.0 };
        vector<double> circumcenter = { dot(arr, inverse[0]),
                                        dot(arr, inverse[1]),
                                        dot(arr, inverse[2]) };
        return circumcenter;
    }

    bool is_on_correct_side(vector<double> A, vector<double> B, vector<double> C,
                            vector<double> test_point, vector<double> circum) {
        vector<double> B_A = { B[0] - A[0], B[1] - A[1], B[2] - A[2] };
        vector<double> C_A = { C[0] - A[0], C[1] - A[1], C[2] - A[2] };
        vector<double> test_A = { test_point[0] - A[0], test_point[1] - A[1], test_point[2] - A[2] };
        vector<double> circum_A = { circum[0] - A[0], circum[1] - A[1], circum[2] - A[2] };
        vector<double> cross_BA_CA = cross_3D(B_A, C_A);
        double dot_test = dot(cross_BA_CA, test_A);
        double dot_circum = dot(cross_BA_CA, circum_A);

        if (dot_test >= 0 && dot_circum >= 0) {
            return true;
        } else if (dot_test <= 0 && dot_circum <= 0) {
            return true;
        } else {
            return false;
        }
    }

    bool circumcenter_is_inside(vector<vector<double>> points) {
        vector<double> A = points[0];
        vector<double> B = points[1];
        vector<double> C = points[2];
        vector<double> D = points[3];
        vector<double> circum = circumcenter(points);

        if ( !is_on_correct_side(A, B, C, D, circum) ) {
            return false;
        } else if ( !is_on_correct_side(A, B, D, C, circum) ) {
            return false;
        } else if ( !is_on_correct_side(A, C, D, B, circum) ) {
            return false;
        } else if ( !is_on_correct_side(C, B, D, A, circum) ) {
            return false;
        } else {
            return true;
        }
    }

    double cech_parameter(vector<vector<double>> points) {
        vector<double> A = points[0];
        vector<double> B = points[1];
        vector<double> C = points[2];
        vector<double> D = points[3];
        vector<double> circum = circumcenter(points);

        if ( circumcenter_is_inside(points) ) {
            return circumradius(points);
        } else {
            vector<vector<double>> ABC = { A, B, C };
            vector<vector<double>> ABD = { A, B, D };
            vector<vector<double>> ACD = { A, C, D };
            vector<vector<double>> BCD = { B, C, D };

            double param_ABC = delcechfiltr_tri::cech_parameter_3D(ABC);
            double param_ABD = delcechfiltr_tri::cech_parameter_3D(ABD);
            double param_ACD = delcechfiltr_tri::cech_parameter_3D(ACD);
            double param_BCD = delcechfiltr_tri::cech_parameter_3D(BCD);
            vector<double> params_tri = { param_ABC, param_ABD, param_ACD, param_BCD};
            sort(params_tri.begin(), params_tri.end());
            return params_tri[3];
        }
    }

    vector<double> cech_param_list_tetrahedra(vector<vector<double>> points,
                                              vector<vector<size_t>> tetra) {
        vector<double> parameters;
        for (size_t ind = 0; ind < tetra.size(); ++ind) {
            vector<size_t> tmp_vert_i = tetra[ind];
            vector<vector<double>> vertices = { points[tmp_vert_i[0]],
                                                points[tmp_vert_i[1]],
                                                points[tmp_vert_i[2]],
                                                points[tmp_vert_i[3]] };
            double tmp_param = cech_parameter(vertices);
            parameters.push_back(tmp_param);
        }
        return parameters;
    }

}
