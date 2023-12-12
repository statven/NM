#include "least_squares.h"

#include <iostream>
using namespace std;
double calculateAverage(const double* values, size_t size) {
    double sum = 0.0;
    for (size_t i = 0; i < size; ++i) {
        sum += values[i];
    }
    return sum / size;
}

double calculateCoefficientA(const double* H, const double* mu, size_t size, double H_avg, double mu_avg) {
    double S_Hmu = 0.0;
    double S_H = 0.0;
    for (size_t i = 0; i < size; ++i) {
        S_Hmu += (H[i] - H_avg) * (mu[i] - mu_avg);
        S_H += (H[i] - H_avg) * (H[i] - H_avg);
    }
    return S_Hmu / S_H;
}

double calculateCoefficientB(double a, double H_avg, double mu_avg) {
    return mu_avg - a * H_avg;
}


int main() {
    setlocale(LC_ALL, "Russian");
    double H[N] = { 0.164, 0.328, 0.656, 0.984, 1.312, 1.640 }; // Массив значений H
    double mu[N] = { 0.448, 0.432, 0.421, 0.417, 0.414, 0.412 }; // Массив значений mu

    performMNKApproximation(H, mu,1);
    size_t size = sizeof(H) / sizeof(H[0]);

    double H_avg = calculateAverage(H, size);
    double mu_avg = calculateAverage(mu, size);

    double a = calculateCoefficientA(H, mu, size, H_avg, mu_avg);
    double b = mu_avg - a * H_avg;
    cout << "The coefficients : a = " << a << ", b = " << b << endl;

    return 0;
}