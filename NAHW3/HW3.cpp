#define _USE_MATH_DEFINES // for C++  

#include<iomanip>
#include<math.h>
#include <iostream>
using namespace std;

double f(double x) {
	return x * (x - 1.0) * (x - 2.0) * (x - 3.0) * (x - 4.0) + 7.0;
}

double integral(double a0, double an) {
	double result_0 = (2 * pow(a0, 6) - 24 * pow(a0, 5) + 105 * pow(a0, 4) - 200 * pow(a0, 3) + 144 * pow(a0, 2) + 84 * a0) / 12;
	double result_1 = (2 * pow(an, 6) - 24 * pow(an, 5) + 105 * pow(an, 4) - 200 * pow(an, 3) + 144 * pow(an, 2) + 84 * an) / 12;
	return result_1 - result_0;
}

double Simpson(double t0, double tn, int n) {
	double integration = f(t0) + f(tn);
	double init = t0;
	double h = (tn - t0) / n;
	for (int i = 1; i < n; i++) {
		init += h;
		if (i % 2 == 0)
			integration += 2 * f(init);
		else
			integration += 4 * f(init);
	}
	integration = integration * h / 3.0;
	return integration;
}

double Trapezoid(double a0, double an, int n) {
	double h = (an - a0) / n;
	double integration = 0;
	double init = a0;
	int count = 0;
	while (init < an) {
		integration += (f(init) + f(init + h)) * h / 2;
		init += h;
	}
	return integration;
}

double Lagrange(double *x, double *y, double input, int n) {
	double sum = 0.0;
	double temp;
	for (int i = 0; i < n; i++) {
		temp = y[i];
		for (int j = 0; j < n; j++) {
			if (i != j) temp = temp * ((input - x[j]) / (x[i] - x[j]));
		}
		sum += temp;
	}
	return sum;
}

double conform_map(double t, double a, double b)
{
	double   x;
	x = (b - a)*t / 2.0 + (a + b) / 2.0;
	return x;
}

double GaussQuadrature(double a0, double an, int num) {
	//Sample points in [-1, 1], 1-, 2-, 3-, and 4-points. Ignore extra points (the 0s).
	double  p[4][4] = { { 0.0, 0.0, 0.0, 0.0 },{ 0.5773502691896257, -0.5773502691896257, 0.0, 0.0 },
	{ 0.0, 0.7745966692414834, -0.7745966692414834, 0.0 },
	{ 0.3399810435848563, -0.3399810435848563, 0.8611363115940526, -0.8611363115940526 } };

	//weights, for 1-, 2-, 3-, and 4-points. Extra weights = 0s.
	double  wgt[4][4] = { { 2.0, 0.0, 0.0, 0.0 },{ 1.0, 1.0, 0.0, 0.0 },
	{ 0.888888888888888888888, 0.555555555555555556, 0.555555555555555556, 0.0 },
	{ 0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538 } };
	double   integration = 0, range = (an - a0) / 8, init = a0;	
	for (int j = 0; j < 8; j++) {
		for (int i = 0; i < num; i++) {
			double point = conform_map(p[num-1][i], init, init + range);
			integration += wgt[num-1][i] * f(point);
		}
		init += range;
	}
	integration *= range / 2;
	return integration;
}

void main() {
	double result, exact_result, err;
	const double up = -1, low = 4;
	int n;
	exact_result = integral(up, low);
	n = 8;
	//problem_3
	cout << "Problem3" << endl << endl;

	//n = 8, Simpson's Rule
	result = Simpson(up, low, n);
	cout << "N = 8, Simpson's Rule : " << endl;
	cout << "result:" << fixed << setprecision(15) << result << ", exact solution:" << fixed << setprecision(15) << exact_result << ", error:" << fixed << setprecision(15) << exact_result - result << endl << endl;

	//n = 8, Trapezoid's Rule
	result = Trapezoid(up, low, n);
	cout << "N = 8, Trapezoid's Rule : " << endl;
	cout << "result:" << fixed << setprecision(15) << result << ", exact solution:" << fixed << setprecision(15) << exact_result << ", error:" << fixed << setprecision(15) << exact_result - result << endl << endl;

	//problem_4
	cout << "Problem4" << endl << endl;
	//Simpson
	while (true) {
		result = Simpson(up, low, n);
		err = exact_result - result;
		 
		if (err < 0.00001){
			if (err < 0)
				err *= (-1);
			break;
		}
		cout << "Simpson's rule : N = " << n << endl;
		cout << "result:" << fixed << setprecision(15) << result << ", exact solution:" << fixed << setprecision(15) << exact_result << ", error:" << fixed << setprecision(15) << err << endl << endl;
		n *= 2;
	}
	

	//Trapezoid
	n = 8;
	while (true) {
		result = Trapezoid(up, low, n);
		err = exact_result - result;
		if (err < 0.00001) {
			if (err < 0)
				err *= (-1);
			break;
		}
		cout << "Trapezoid's rule : N = " << n << endl;
		cout << "result:" << fixed << setprecision(15) << result << ", exact solution:" << fixed << setprecision(15) << exact_result << ", error:" << fixed << setprecision(15) << err << endl << endl;
		n *= 2;
	}


	//problem_5
	cout << "Problem5" << endl << endl;
	for (n = 2; n <= 4; n++) {
		result = GaussQuadrature(up, low, n);
		err = exact_result - result;
		cout << "Gaussian quadrature : sample_point = " << n << endl;
		cout << "result:" << fixed << setprecision(15) << result << ", exact solution:" << fixed << setprecision(15) << exact_result << ", error:" << fixed << setprecision(15) << err << endl << endl;

	}
	system("pause");
}