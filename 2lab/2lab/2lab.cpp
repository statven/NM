#include <iostream>
#include <cmath>

using namespace std;

const double eps1 = 1e-9;
const double eps2 = 10.0;
const int NIT = 1000;

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
    return x2 * x2 * 3 * x1;
}

double df2_dx2(double x1, double x2) {
    return x1 * 3 * pow(x2, 2) - 1;
}

int main() {
    double x1 = 1, x2 = 1;
    int k = 0;
    cout << "k\t\tx1\t\tx2\t\tf1\t\tf2\n";
    while (true) {
        double f_1 = f1(x1, x2);
        double f_2 = f2(x1, x2);
        double df[2][2] = { {df1_dx1(x1, x2), df1_dx2(x1, x2)}, {df2_dx1(x1, x2), df2_dx2(x1, x2)} };
        double det = df[0][0] * df[1][1] - df[0][1] * df[1][0];
        double dx1 = (df[1][1] * f_1 - df[0][1] * f_2) / det;
        double dx2 = (-df[1][0] * f_1 + df[0][0] * f_2) / det;
        x1 += dx1;
        x2 += dx2;
        k++;
        cout << k << "\t\t" << x1 << "\t\t" << x2 << "\t\t" << f_1 << "\t\t" << f_2 << endl;
        if (abs(f_1) < eps1 && abs(f_2) < eps1) {
            cout << "Solution found: x1 = " << x1 << ", x2 = " << x2 << endl;
            break;
        }
        if (k >= NIT) {// проверяем условие максимального числа итераций
            cout << "Maximum number of iterations " << endl;
            break;
        }
        if (abs(dx1) > eps2 || abs(dx2) > eps2) {// проверяем условие сходимости
            cout << "Divergence detected." << endl;// выводим сообщение о расходимости
            break;
        }
    }
    return 0;
}