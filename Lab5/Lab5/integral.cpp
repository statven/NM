#include "integral.h"
#include <cmath>

// Подынтегральная функция
double f(double x, double y) {
    return x * x + 2 * y;
}

// Метод Симпсона для вычисления двойного интеграла
double doubleIntegralSimpson(double a, double b, double c, double d, int n) {
    double h = (b - a) / n;
    double k = (d - c) / n;
    double sum = 0;

    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            double x = a + i * h;
            double y = c + j * k;
            double coeff = 1.0;

            if (i == 0 || i == n) {
                coeff *= 1.0;
            }
            else if (i % 2 == 1) {
                coeff *= 4.0;
            }
            else {
                coeff *= 2.0;
            }

            if (j == 0 || j == n) {
                coeff *= 1.0;
            }
            else if (j % 2 == 1) {
                coeff *= 4.0;
            }
            else {
                coeff *= 2.0;
            }

            sum += coeff * f(x, y);
        }
    }

    return (h * k / 9) * sum;
}


