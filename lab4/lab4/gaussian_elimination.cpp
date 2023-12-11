#include "gaussian_elimination.h"

#include <iostream>
#include <cmath>

void solveGauss(double A[][7], double B[], double X[]) {
   int  N = 6;
    for (int k = 0; k < N - 1; k++) {
        for (int i = k + 1; i < N; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < N; j++) {
                A[i][j] -= factor * A[k][j];
            }
            B[i] -= factor * B[k];
        }
    }

    // Пример реализации метода Гаусса (обратный ход)
    for (int i = N - 1; i >= 0; i--) {
        X[i] = B[i];
        for (int j = i + 1; j < N; j++) {
            X[i] -= A[i][j] * X[j];
        }
        X[i] /= A[i][i];
    }
}