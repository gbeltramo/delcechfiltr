#include <vector>
#include <cmath>   // std::sqrt
// #include <tuple>
// #include <map>
// #include <algorithm>
// #include <limits>


using namespace std;

namespace delcechfiltr_tri {

    double magnitude(vector<double> v) {
        double tot = 0.0;
        for (size_t ind = 0; ind < v.size(); ++ind) {
            tot = tot + v[ind] * v[ind];
        }
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
        vector<double> A = {0.0, 0.0};
        vector<double> B = {0.0, 0.0};
        A[0] = points[0][0] - points[2][0];
        A[1] = points[0][1] - points[2][1];
        B[0] = points[1][0] - points[2][0];
        B[1] = points[1][1] - points[2][1];
        vector<double> A_minus_B = {A[0] - B[0], A[1] - B[1]};
        double mag_A = magnitude(A);
        double mag_B = magnitude(B);
        double mag_A_minus_B = magnitude(A_minus_B);
        double mag_cross_AB = fabs(cross_2D(A, B));
        return mag_A * mag_B * mag_A_minus_B / mag_cross_AB / 2.0;
    }

    double circumradius_3D(vector<vector<double>> points) {
        vector<double> A = {0.0, 0.0, 0.0};
        vector<double> B = {0.0, 0.0, 0.0};
        A[0] = points[0][0] - points[2][0];
        A[1] = points[0][1] - points[2][1];
        A[2] = points[0][2] - points[2][2];
        B[0] = points[1][0] - points[2][0];
        B[1] = points[1][1] - points[2][1];
        B[2] = points[1][2] - points[2][2];
        vector<double> A_minus_B = {A[0] - B[0], A[1] - B[1], A[2] - B[2]};
        double mag_A = magnitude(A);
        double mag_B = magnitude(B);
        double mag_A_minus_B = magnitude(A_minus_B);
        vector<double> cross_AB = cross_3D(A, B);
        double mag_cross_AB = magnitude(cross_AB);
        return mag_A * mag_B * mag_A_minus_B / mag_cross_AB / 2.0;
    }

}
