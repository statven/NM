#ifndef NEWTON_METHOD_H
#define NEWTON_METHOD_H

#include <iostream>
#include <cmath>

using namespace std;

void solve(double x1, double x2, double M, bool analytical_jacobian,
    double (*f1)(double x1, double x2), double (*f2)(double x1, double x2),
    double (*df1_dx1)(double x1, double x2), double (*df1_dx2)(double x1, double x2),
    double (*df2_dx1)(double x1, double x2), double (*df2_dx2)(double x1, double x2),
    const double eps1 = 1e-9,
    const double eps2 = 10,
    const int NIT = 100);

#endif // NEWTON_METHOD_H