#include "newton_method.h"
#include <iostream>
#include <iomanip>
#include<algorithm>
#include <math.h>


using namespace std;
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
	return pow(x2, 3);
}

double df2_dx2(double x1, double x2) {
	return 3 * pow(x1, 2) * pow(x2, 2) - 1;
}


void newtonMethod(double& x1, double& x2, double epsilon1, double epsilon2, int NIT) {

	int k = 1;

	cout << "k" << setw(12) << "l1" << setw(14) << "l2" << setw(16) << "x1" << setw(14) << "x2" << endl;
	cout << "------------------------------------------------------------------" << endl;
	
	double* dXk = new double[n];
	double dXk1 = 0, dXk2 = 0;
	double l1 = 1, l2 = 1, l1_2 = 0, l2_2 = 0;
	cout << k << setw(10) << l1 << setw(10) << l2 << setw(10) << x1 << setw(10) << x2 << endl;
	try {
		while (l1 > epsilon1 || l2 > epsilon2) {
			double Fx[] = { -f2(x1,x2), -f2(x1,x2) };
			double Jacobian[n][n] = {
					{df1_dx1(x1,x2), df1_dx2(x1,x2)},
					{df2_dx1(x1,x2), df2_dx2(x1,x2)}
			};
			
		
			dXk = methodGauss(Jacobian, Fx);
			dXk1 = dXk[0] + x1;
			dXk2 = dXk[1] + x2;
			l1 = abs(f1(x1, x2));
			l1_2 = abs(f2(x1, x2));
			l1_2 > l1 ? l1 = l1_2 : l1_2;
			dXk1 >= 1 ?
				l2 = abs(dXk1 - x1) / dXk1 :
				l2 = abs(dXk1 - x1);
			dXk2 >= 1 ?
				l2_2 = abs(dXk2 - x2) / dXk2 :
				l2_2 = abs(dXk2 - x2);
			if (l2_2 > l2)
				l2 = l2_2;
			x1 = dXk1; x2 = dXk2;
			dXk[0] = x1; dXk[1] = x2;
			cout << k << "\t " << l1 << " \t" << l2 << " \t" << x1 << " \t" << x2 << endl;
			k++;
			if (k > NIT)
			{

				throw 1;

			}
			if (l1 <= epsilon1 && l2 <= epsilon2) {
				cout << "solution!!!" << endl;
				break;
			}
		}
	}
	catch (int err) {
		cout << endl;
		cout << "IER=2" << endl;

	}
	k = 0;

}
double* methodGauss(double mA[n][n], double cB[n]) {
	double X[n] = { 0 };
	for (int i = 0; i < n; i++) {
		int maxIndex = i;
		double max = mA[i][i];
		for (int j = i + 1; j < n; j++) {
			if (abs(max) < abs(mA[j][i])) {
				maxIndex = j;
				max = mA[j][i];
			}
		}
		if (i != maxIndex) {
			double root = cB[i];
			cB[i] = cB[maxIndex];
			cB[maxIndex] = root;
			for (int j = 0; j < n; j++) {
				double x = mA[i][j];
				mA[i][j] = mA[maxIndex][j];
				mA[maxIndex][j] = x;
			}
		}
		double a = mA[i][i];
		for (int j = i; j < n; j++)
		{
			mA[i][j] /= a;
		}
		cB[i] /= a;
		for (int j = i + 1; j < n; j++)
		{
			double s = mA[j][i];
			for (int k = i; k < n; k++)
			{
				mA[j][k] -= s * mA[i][k];
			}
			cB[j] -= s * cB[i];
		}
	}

	for (int k = n - 1; k >= 0; k--)
	{
		X[k] = cB[k];
		for (int i = n - 1; i > k; i--)
		{
			X[k] -= mA[k][i] * X[i];
		}
	}
	cout << endl;
	return X;
}
void solutionM(double x1, double x2, double M,double e1, double e2, int NIT) {
	int k = 1;
	cout << "k" << setw(12) << "l1" << setw(14) << "l2" << setw(16) << "x1" << setw(14) << "x2" << setw(14) << "M" << endl;
	cout << "------------------------------------------------------------------" << endl;
	double* dXk = new double[n];
	double dXk1 = 0, dXk2 = 0;
	double l1 = 1, l2 = 1, l1_2 = 0, l2_2 = 0;
	cout << k << setw(10) << l1 << setw(10) << l2 << setw(10) << x1 << setw(10) << x2 << endl;

	try {
		while (l1 > e1 || l2 > e2) {
			double Fx[] = { -f1(x1,x2), -f2(x1,x2) };
			double Jacobian[n][n] = {
				{(f1(x1 + M * x1, x2) - f1(x1, x2)) / (M * x1), (f1(x1, x2 + M * x2) - f1(x1, x2)) / (M * x2)},
				{(f2(x1 + M * x1, x2) - f2(x1, x2)) / (M * x1), (f2(x1, x2 + M * x2) - f2(x1, x2)) / (M * x2)}
			};
			dXk = methodGauss(Jacobian, Fx);
			dXk1 = dXk[0] + x1;
			dXk2 = dXk[1] + x2;
			l1 = abs(f1(x1, x2));
			l1_2 = abs(f2(x1, x2));
			if (l1_2 > l1) {
				l1 = l1_2;
			}
			dXk1 >= 1 ?
				l2 = abs(dXk1 - x1) / dXk1 :
				l2 = abs(dXk1 - x1);
			dXk2 >= 1 ?
				l2_2 = abs(dXk2 - x2) / dXk2 :
				l2_2 = abs(dXk2 - x2);
			if (l2_2 > l2)
				l2 = l2_2;
			x1 = dXk1; x2 = dXk2;
			dXk[0] = x1; dXk[1] = x2;
			cout << k << "\t " << l1 << " \t" << l2 << " \t" << x1 << " \t" << x2 << " \t" << M << endl;
			k++;
			if (k >= NIT)
			{
				throw 1;
			}
			if (l1 <= e1 && l2 <= e2) {
				break;
			}
		}
	}
	catch (int e) {
		cout << endl;
		cout << "IER=2" << endl;
	}
	k = 0;

}