#include "newton_method.h"

void solve(double x1, double x2, double M, bool analytical_jacobian,
    double (*f1)(double x1, double x2), double (*f2)(double x1, double x2),
    double (*df1_dx1)(double x1, double x2), double (*df1_dx2)(double x1, double x2),
    double (*df2_dx1)(double x1, double x2), double (*df2_dx2)(double x1, double x2),
    const double eps1,
    const double eps2,
    const int NIT) {

    cout << "M = " << M << endl;
    cout << "Analytical Jacobian: " << (analytical_jacobian ? "yes" : "no") << endl;
    cout << "Initial approximation: x1 = " << x1 << ", x2 = " << x2 << endl;
    cout << "eps1 = " << eps1 << ", eps2 = " << eps2 << ", NIT = " << NIT << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "k\t\tx1\t\tx2\t\t||F||\t\t||dx||\t\t||x||" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;

    double x1_prev = x1, x2_prev = x2;
    double F1, F2, J11, J12, J21, J22, detJ, dx1, dx2, x_norm, F_norm;
    int k = 0;

    do {
        F1 = f1(x1_prev, x2_prev);
        F2 = f2(x1_prev, x2_prev);

        if (analytical_jacobian) {
            J11 = df1_dx1(x1_prev, x2_prev);
            J12 = df1_dx2(x1_prev, x2_prev);
            J21 = df2_dx1(x1_prev, x2_prev);
            J22 = df2_dx2(x1_prev, x2_prev);
        }
        else {
            double h = 1e-6;
            J11 = (f1(x1_prev + h, x2_prev) - F1) / h;
            J12 = (f1(x1_prev, x2_prev + h) - F1) / h;
            J21 = (f2(x1_prev + h, x2_prev) - F2) / h;
            J22 = (f2(x1_prev, x2_prev + h) - F2) / h;
        }

        detJ = J11 * J22 - J12 * J21;
        dx1 = (-F1 * J22 + F2 * J12) / detJ;
        dx2 = (F1 * J21 - F2 * J11) / detJ;

        double x1_next = x1_prev + dx1;
        double x2_next = x2_prev + dx2;

        x_norm = sqrt(x1_next * x1_next + x2_next * x2_next);
        F_norm = sqrt(F1 * F1 + F2 * F2);

        cout << k << "\t\t" << x1_prev << "\t\t" << x2_prev << "\t\t" << F_norm << "\t\t" << sqrt(dx1 * dx1 + dx2 * dx2) << "\t\t" << x_norm << endl;

        if (F_norm < eps1 && x_norm < eps2) {
            cout << "--------------------------------------------------------------------------------" << endl;
            cout << "Solution: x1 = " << x1_next << ", x2 = " << x2_next << endl;
            return;
        }

        x1_prev = x1_next;
        x2_prev = x2_next;

        k++;
    } while (k < NIT);

    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Solution not found within " << NIT << " iterations." << endl;
}