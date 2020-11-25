#ifndef DELCECHFILTR_TRIANGLE_HPP
#define DELCECHFILTR_TRIANGLE_HPP

using namespace std;

namespace delcechfiltr_tri {
    vector<double> cross_3D(vector<double> A, vector<double> B);
    double circumradius_2D(vector<vector<double>> points);
    double circumradius_3D(vector<vector<double>> points);
}

#endif
