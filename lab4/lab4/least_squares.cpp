#include "least_squares.h"
#include "gaussian_elimination.h"
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
void performMNKApproximation(double x[6], double y[6], int m) {
    const int M = 1;
    double POWERX[2 * M] = { 0 }; // Массив для хранения сумм (x_i)^k
    double SUMX[M + 1][M+1] = { 0 }; // Матрица коэффициентов
    double PRAW[M + 1] = { 0 }; // Правая часть системы
    double S2 = 0;
    cout << "H:" << endl;
    for (int i = 0; i < N; i++) {
        cout << x[i] << setw(5);
        if (i % 5 == 0 && i != 0)
            cout << endl;
    }
    cout << endl;
    cout << "mu:" << endl;
    for (int i = 0; i < N; i++) {
        cout << y[i] << "  ";
        if (i % 5 == 0 && i != 0)
            cout << endl;
    }
    cout << endl;
    cout << "POWERX:" << endl;
    for (int k = 0; k < 2 * m; k++)
    {
        POWERX[k] = 0;
        for (int i = 0; i < N; i++)
        {
            POWERX[k] += pow(x[i], k + 1);
        }
        cout << POWERX[k] << " ";
    }
    cout << endl;
    cout << endl;

    // Вычисление сумм (x_i)^k
    for (int i = 0; i < m + 1; i++)
        for (int j = 0; j < m + 1; j++)
            i + j ? SUMX[i][j] = POWERX[i + j - 1] : SUMX[0][0] = N;

    cout << "SUMX:" << endl;
    for (int i = 0; i < m + 1; i++) {
        for (int j = 0; j < m + 1; j++)
            cout << SUMX[i][j] << " ";
        cout << endl;
    }
    cout << endl;

    for (int i = 0; i < m + 1; i++)
    {
        PRAW[i] = 0;
        for (int j = 0; j < N; j++)
        {
            PRAW[i] += y[j] * pow(x[j], i);
        }
    }

  
    cout << "PRAW:" << endl;
    for (int i = 0; i < m + 1; i++)
    {
        cout << PRAW[i] << " ";
    }
    cout << endl;

    double* a = methodGauss(SUMX, PRAW);

    cout << "a:" << endl;
    for (int i = 0; i < m + 1; i++) {
        cout << a[i] << " ";
    }
    cout << endl;

    for (int i = 0; i < N; i++)
    {
        double sum = y[i];
        for (int j = 0; j < m + 1; j++)
        {
            sum -= a[j] * pow(x[i], j);
        }
        S2 += pow(sum, 2);
    }
    S2 /= N - m - 1;
    double sigma = sqrt(S2); 
}