#include "least_squares.h"

#include <iostream>

int main() {
    double H[N] = { 0.164, 0.328, 0.656, 0.984, 1.312, 1.640 }; // Массив значений H
    double mu[N] = { 0.448, 0.432, 0.421, 0.417, 0.414, 0.412 }; // Массив значений mu

    performMNKApproximation(H, mu,6);

  
    return 0;
}