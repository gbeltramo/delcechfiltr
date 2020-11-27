#include <vector>
#include <cmath>       // std::sqrt
#include <algorithm>   // std::sort

using namespace std;

namespace delcechfiltr_tri {

    double dot(vector<double> A, vector<double> B) {
        double tot = 0.0;
        for (size_t ind = 0; ind < A.size(); ++ind) {
            tot = tot + A[ind] * B[ind];
        }
        return tot;
    }

    double euclidean(vector<double> A, vector<double> B) {
        double tot = 0.0;
        for (size_t ind = 0; ind < A.size(); ++ind) {
            tot = tot + (A[ind] - B[ind]) * (A[ind] - B[ind]);
        }
        return sqrt(tot);
    }

    double magnitude(vector<double> A) {
        double tot = dot(A, A);
        return sqrt(tot);
    }

    double cross_2D(vector<double> A, vector<double> B) {
        return A[0] * B[1] - A[1] * B[0];
    }

    vector<double> cross_3D(vector<double> A, vector<double> B) {
        vector<double> cross_out = {0.0, 0.0, 0.0};
        cross_out[0] = A[1] * B[2] - A[2] * B[1];
        cross_out[1] = A[2] * B[0] - A[0] * B[2];
        cross_out[2] = A[0] * B[1] - A[1] * B[0];
        return cross_out;
    }

    double circumradius_2D(vector<vector<double>> points) {
        vector<double> A = {points[0][0] - points[2][0],
                            points[0][1] - points[2][1]};
        vector<double> B = {points[1][0] - points[2][0],
                            points[1][1] - points[2][1]};
        vector<double> A_minus_B = {A[0] - B[0], A[1] - B[1]};
        double mag_A = magnitude(A);
        double mag_B = magnitude(B);
        double mag_A_minus_B = magnitude(A_minus_B);
        double mag_cross_AB = fabs(cross_2D(A, B));
        return mag_A * mag_B * mag_A_minus_B / mag_cross_AB / 2.0;
    }

    double circumradius_3D(vector<vector<double>> points) {
        vector<double> A = {points[0][0] - points[2][0],
                            points[0][1] - points[2][1],
                            points[0][2] - points[2][2]};
        vector<double> B = {points[1][0] - points[2][0],
                            points[1][1] - points[2][1],
                            points[1][2] - points[2][2]};
        vector<double> A_minus_B = {A[0] - B[0], A[1] - B[1], A[2] - B[2]};
        double mag_A = magnitude(A);
        double mag_B = magnitude(B);
        double mag_A_minus_B = magnitude(A_minus_B);
        vector<double> cross_AB = cross_3D(A, B);
        double mag_cross_AB = magnitude(cross_AB);
        return mag_A * mag_B * mag_A_minus_B / mag_cross_AB / 2.0;
    }

    vector<double> circumcenter_2D(vector<vector<double>> points) {
        vector<double> A = {points[0][0] - points[2][0],
                            points[0][1] - points[2][1]};
        vector<double> B = {points[1][0] - points[2][0],
                            points[1][1] - points[2][1]};
        double mag_cross_AB = fabs(cross_2D(A, B));
        double mag_A2 = magnitude(A) * magnitude(A);
        double mag_B2 = magnitude(B) * magnitude(B);
        vector<double> numer1 = { mag_A2 * B[0] - mag_B2 * A[0],
                                  mag_A2 * B[1] - mag_B2 * A[1] };
        vector<double> numer2 = { dot(numer1, B) * A[0] - dot(numer1, A) * B[0],
                                  dot(numer1, B) * A[1] - dot(numer1, A) * B[1] };
        double coeff = 2.0 * mag_cross_AB * mag_cross_AB;
        vector<double> circ = { numer2[0] / coeff + points[2][0],
                                numer2[1] / coeff + points[2][1] };
        return circ;
    }

    vector<double> circumcenter_3D(vector<vector<double>> points) {
        vector<double> A = {points[0][0] - points[2][0],
                            points[0][1] - points[2][1],
                            points[0][2] - points[2][2]};
        vector<double> B = {points[1][0] - points[2][0],
                            points[1][1] - points[2][1],
                            points[1][2] - points[2][2]};
        vector<double> cross_AB = cross_3D(A, B);
        double mag_cross_AB = magnitude(cross_AB);
        double mag_A2 = magnitude(A) * magnitude(A);
        double mag_B2 = magnitude(B) * magnitude(B);
        vector<double> numer1 = { mag_A2 * B[0] - mag_B2 * A[0],
                                  mag_A2 * B[1] - mag_B2 * A[1],
                                  mag_A2 * B[2] - mag_B2 * A[2] };
        vector<double> numer2 = { dot(numer1, B) * A[0] - dot(numer1, A) * B[0],
                                  dot(numer1, B) * A[1] - dot(numer1, A) * B[1],
                                  dot(numer1, B) * A[2] - dot(numer1, A) * B[2] };
        double coeff = 2.0 * mag_cross_AB * mag_cross_AB;
        vector<double> circ = { numer2[0] / coeff + points[2][0],
                                numer2[1] / coeff + points[2][1],
                                numer2[2] / coeff + points[2][2] };
        return circ;
    }

    bool triangle_is_not_obtuse(vector<double> lengths) {
        double a2 = lengths[0] * lengths[0];
        double b2 = lengths[1] * lengths[1];
        double c2 = lengths[2] * lengths[2];

        bool ret = a2 + b2 > c2;
        return ret;
    }

    double cech_parameter_2D(vector<vector<double>> points) {
        double len_AB = euclidean(points[0], points[1]);
        double len_AC = euclidean(points[0], points[2]);
        double len_BC = euclidean(points[1], points[2]);
        vector<double> lengths = {len_AB, len_AC, len_BC};
        sort(lengths.begin(), lengths.end());

        if (triangle_is_not_obtuse(lengths)) {
            return circumradius_2D(points);
        } else {
            return lengths[2] / 2.0;
        }
    }

    double cech_parameter_3D(vector<vector<double>> points) {
        double len_AB = euclidean(points[0], points[1]);
        double len_AC = euclidean(points[0], points[2]);
        double len_BC = euclidean(points[1], points[2]);
        vector<double> lengths = {len_AB, len_AC, len_BC};
        sort(lengths.begin(), lengths.end());

        if (triangle_is_not_obtuse(lengths)) {
            return circumradius_3D(points);
        } else {
            return lengths[2] / 2.0;
        }
    }

    vector<double> cech_param_list_triangles_2D(vector<vector<double>> points,
                                                vector<vector<size_t>> triangles) {
        vector<double> parameters;
        for (size_t ind = 0; ind < triangles.size(); ++ind) {
            vector<size_t> tmp_vert_i = triangles[ind];
            vector<vector<double>> vertices = { points[tmp_vert_i[0]],
                                                points[tmp_vert_i[1]],
                                                points[tmp_vert_i[2]] };
            double tmp_param = cech_parameter_2D(vertices);
            parameters.push_back(tmp_param);
        }
        return parameters;
    }

    vector<double> cech_param_list_triangles_3D(vector<vector<double>> points,
                                                vector<vector<size_t>> triangles) {
        vector<double> parameters;
        for (size_t ind = 0; ind < triangles.size(); ++ind) {
            vector<size_t> tmp_vert_i = triangles[ind];
            vector<vector<double>> vertices = { points[tmp_vert_i[0]],
                                                points[tmp_vert_i[1]],
                                                points[tmp_vert_i[2]] };
            double tmp_param = cech_parameter_3D(vertices);
            parameters.push_back(tmp_param);
        }
        return parameters;
    }

    vector<double>  miniball_center_2D(vector<vector<double>> points) {
        vector<double> A = points[0];
        vector<double> B = points[1];
        vector<double> C = points[2];

        vector<double> lengths = { euclidean(A, B),
                                   euclidean(A, C),
                                   euclidean(B, C) };
        sort(lengths.begin(), lengths.end());

        if ( triangle_is_not_obtuse(lengths) ) {
            return circumcenter_2D(points);
        } else {
            vector<double> lengths = { euclidean(A, B),
                                       euclidean(A, C),
                                       euclidean(B, C) };
            if (lengths[0] >= lengths[1] && lengths[0] >= lengths[2]) {
                return { (A[0] + B[0]) / 2.0,
                         (A[1] + B[1]) / 2.0 };
            } else if (lengths[1] >= lengths[0] && lengths[1] >= lengths[2]) {
                return { (A[0] + C[0]) / 2.0,
                         (A[1] + C[1]) / 2.0 };
            } else {
                return { (B[0] + C[0]) / 2.0,
                         (B[1] + C[1]) / 2.0 };
            }
        }
    }

    vector<double>  miniball_center_3D(vector<vector<double>> points) {
        vector<double> A = points[0];
        vector<double> B = points[1];
        vector<double> C = points[2];

        vector<double> lengths = { euclidean(A, B),
                                   euclidean(A, C),
                                   euclidean(B, C) };
        sort(lengths.begin(), lengths.end());

        if ( triangle_is_not_obtuse(lengths) ) {
            return circumcenter_3D(points);
        } else {
            vector<double> lengths = { euclidean(A, B),
                                       euclidean(A, C),
                                       euclidean(B, C) };
            if (lengths[0] >= lengths[1] && lengths[0] >= lengths[2]) {
                return { (A[0] + B[0]) / 2.0,
                         (A[1] + B[1]) / 2.0,
                         (A[2] + B[2]) / 2.0 };
            } else if (lengths[1] >= lengths[0] && lengths[1] >= lengths[2]) {
                return { (A[0] + C[0]) / 2.0,
                         (A[1] + C[1]) / 2.0,
                         (A[2] + C[2]) / 2.0 };
            } else {
                return { (B[0] + C[0]) / 2.0,
                         (B[1] + C[1]) / 2.0,
                         (B[2] + C[2]) / 2.0 };
            }
        }
    }

}
