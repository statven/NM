#include "least_squares.h"
#include "gaussian_elimination.h"

#include <iostream>
#include <cmath>

void performMNKApproximation(double x[6], double y[6], int m) {
    const int M = 6;
    double POWERX[2 * M] = { 0 }; // ������ ��� �������� ���� (x_i)^k
    double SUMX[M + 1][7] = { 0 }; // ������� �������������
    double PRAW[M + 1] = { 0 }; // ������ ����� �������
    double a[M + 1] = { 0 }; // ������������ �������������

    // ���������� ���� (x_i)^k
    for (int i = 0; i < N; i++) {
        for (int k = 1; k <= m; k++) {
            POWERX[k - 1] += std::pow(x[i], k);
        }
    }

    // ������������ ������� �������������
    for (int l = 1; l <= m + 1; l++) {
        for (int j = 1; j <= m + 1; j++) {
            SUMX[l - 1][j - 1] = POWERX[l + j - 2];
        }
    }

    // ������������ ������ ����� �������
    for (int l = 1; l <= m + 1; l++) {
        for (int i = 0; i < N; i++) {
            PRAW[l - 1] += std::pow(x[i], l - 1) * y[i];
        }
    }

    // ������� ������� �������� ��������� ������� ������
    solveGauss(SUMX, PRAW, a);

    // ���������� ���������� ���������
    double S = 0.0;
    for (int i = 0; i < N; i++) {
        double fi_x = 0.0;
        for (int k = 0; k <= m; k++) {
            fi_x += a[k] * std::pow(x[i], k);
        }
        S += std::pow(y[i] - fi_x, 2);
    }
    S /= (N - m - 1);
    double RMSE = std::sqrt(S);

    // ����� ������������� � ��������������������� ����������
    std::cout << "������������ ������������� (a_0, a_1, ..., a_m):" << std::endl;
    for (int k = 0; k <= m; k++) {
        std::cout << "a_" << k << " = " << a[k] << std::endl;
    }
    std::cout << "�������������������� ����������: " << RMSE << std::endl;

   
}