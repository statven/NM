#include <iostream>
#include <cmath>
#include "newton_method.h"

using namespace std;

const double eps1 = 1e-9; // заданная точность
const double eps2 = 10; // заданная точность
const int NIT = 100; // максимальное количество итераций

double f1(double x1, double x2) {
    return 1.5 * pow(x1, 3) - pow(x2, 2) - 1;
}

double f2(double x1, double x2) {
    return x1 * pow(x2, 3) - x2 - 4;
}

double df1_dx1(double x1, double x2) {
    return 4.5 * pow(x1, 2);
}

double df1_dx2(double x1, double x2) {
    return -2 * x2;
}

double df2_dx1(double x1, double x2) {
    return pow(x2, 3);
}

double df2_dx2(double x1, double x2) {
    return x1 * 3 * pow(x2, 2) - 1;
}

int main() { 
    double x1 = 1, x2 = 1;

    solve(x1, x2, 0.01, true, f1, f2, df1_dx1, df1_dx2, df2_dx1, df2_dx2, eps1, eps2, NIT);
    solve(x1, x2, 0.01, false, f1, f2, df1_dx1, df1_dx2, df2_dx1, df2_dx2, eps1, eps2, NIT);
    solve(x1, x2, 0.05, true, f1, f2, df1_dx1, df1_dx2, df2_dx1, df2_dx2, eps1, eps2, NIT);
    solve(x1, x2, 0.05, false, f1, f2, df1_dx1, df1_dx2, df2_dx1, df2_dx2, eps1, eps2, NIT);
    solve(x1, x2, 0.1, true, f1, f2, df1_dx1, df1_dx2, df2_dx1, df2_dx2, eps1, eps2, NIT);
    solve(x1, x2, 0.1, false, f1, f2, df1_dx1, df1_dx2, df2_dx1, df2_dx2, eps1, eps2, NIT);

    return 0;
}