#include <iostream>
#include <cmath>
#include "trapezoidal_integration.h"
#include "simpson_integration.h"
#include "integral.h"
using namespace std;
int main() {
    setlocale(LC_ALL, "RUS");
    double a = 1.0; // Нижний предел интегрирования
    double b = 2.631; // Верхний предел интегрирования
    int n = 100,m=100; // Количество разбиений отрезка
    double eps1 = pow(10, -4); // Критерий завершения для eps = 10^(-4)
    double eps2 = pow(10, -5); // Критерий завершения для eps = 10^(-5)

    double result_simpson_1, result_simpson_2, result_trapezoidal_1, result_trapezoidal_2;

    do {
        result_simpson_1 = simpsonIntegration(a, b, n);
        result_simpson_2 = simpsonIntegration(a, b, n * 2);
        n *= 2;
    } while (abs(result_simpson_2 - result_simpson_1) > 3 * eps1);
    cout << "Результат интегрирования методом Симпсона: " << result_simpson_2 << endl;
    n = 100;
    do {
        result_trapezoidal_1 = trapezoidalIntegration(a, b, n);
        result_trapezoidal_2 = trapezoidalIntegration(a, b, n * 2);
        n *= 2;
    } while (abs(result_trapezoidal_2 - result_trapezoidal_1) > 15 * eps2);

    cout << "Результат интегрирования методом трапеции: " << result_trapezoidal_1 << endl;
    n = 10;
    double result = doubleIntegralSimpson(0, 2, 0, 1, n,m);
    cout << "Результат вычисления двойного интеграла методом Симпсона: " << result << endl;
    
    return 0;
}
