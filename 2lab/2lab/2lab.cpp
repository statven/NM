#include <iostream>
#include "newton_method.h"

using namespace std;

int main() {
    double x1 = 1, x2 = 1; // Начальное приближение
    double epsilon1 = 1e-9, epsilon2 = 1e-9; // Точность
    int NIT = 15; // Предельное число итераций

    newtonMethod(x1, x2, epsilon1, epsilon2, NIT);

    return 0;
}
