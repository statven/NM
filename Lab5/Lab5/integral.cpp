#include "integral.h"
#include <cmath>

// ��������������� �������
double f(double x, double y) {
    return x * x + 2 * y;
}

// ����� �������� ��� ���������� �������� ���������
double doubleIntegralSimpson(double a, double b, double c, double d, int n,int m) {
    double h = (b - a) / (2*n);
    double k = (d - c) / m;
    double sum = 0;

    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
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

            if (j == 0 || j == m) {
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


