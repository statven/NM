#pragma once
#ifndef NEWTON_METHOD_H
#define NEWTON_METHOD_H
const int n = 2;
double f1(double x1, double x2);
double f2(double x1, double x2);
double df1_dx1(double x1, double x2);
double df1_dx2(double x1, double x2);
double df2_dx1(double x1, double x2);
double df2_dx2(double x1, double x2);
void newtonMethod(double& x1, double& x2, double epsilon1, double epsilon2, int NIT);
double* methodGauss(double mA[n][n], double cB[n]);

#endif // NEWTON_METHOD_H

