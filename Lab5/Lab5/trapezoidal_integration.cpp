#include "trapezoidal_integration.h"
#include <cmath>

double f2(double x) {
    return ((1 + x + x * x) / pow((x * x * x - 1), 1 / 2));
}

double trapezoidalIntegration(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f2(a) + f2(b));
    for (int i = 1; i < n; i++) {
        sum += f2(a + i * h);
    }
    return h * sum;
}
