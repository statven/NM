#include "simpson_integration.h"
#include <cmath>

double f(double x) {
    return ((1+x+x*x)/pow((x*x*x-1),1/2)); 
}

double simpsonIntegration(double a, double b, int n) {
    double h = (b - a) / (2 * n);
    double sum = f(a) + f(b);

    for (int i = 1; i < (2 * n); i += 2) {
        sum += 4 * f(a + i * h);
    }

    for (int i = 2; i < (2 * n) - 1; i += 2) {
        sum += 2 * f(a + i * h);
    }

    return (h / 3) * sum;
}
