#define _USE_MATH_DEFINES // for C++  

#include<iomanip>
#include<math.h>
#include <iostream>
using namespace std;

double r = 10.0;
bool LagrangeFlag = true; //true : normal Lagrange, false : cubic Lagrange

void Norms_2(int n, double *qx, double *qy) {
	double delta = 2 * M_PI / 360;
	double px[360], py[360], t[360];
	double result = 0.0;
	for (int i = 0; i < 360; i++) {
		t[i] = i * delta;
		px[i] = r * cos(t[i]);
		py[i] = r * sin(t[i]);
		if (i % 30 == 0) {
		/*	cout << fixed << setprecision(3) << i << "  " << fixed << setprecision(6) << px[i] << "  " << fixed << setprecision(6) << qx[i] << "  " << fixed << setprecision(6) << py[i] << "  " << fixed << setprecision(6) << qy[i] << endl;
			//左逼近
			cout << "Left  :" << fixed << setprecision(6) << (qy[i - 1] - py[i]) / (qx[i - 1] - px[i]) << endl;
			//右逼近
			cout << "Right  :" << fixed << setprecision(6) << (qy[i - 1] - py[i]) / (qx[i - 1] - px[i]) << endl;*/
		}
	}
	for (int j = 0; j < 360; j++)
		result += (px[j] - qx[j])*(px[j] - qx[j]) + (py[j] - qy[j])*(py[j] - qy[j]);	
	result = sqrt(result);//Square root
	/*cout << "--------------------------------------------------" << endl;
	cout << "Norm_2 :" << fixed << setprecision(9) << result << endl;*/
}

void NormsInfi(int n, double *qx, double *qy) {
	double delta = 2 * M_PI / 360;
	double px[360], py[360], t[360];
	double max = 0.0, angle, x, y;
	int index;
	for (int i = 0; i < 360; i++) {
		t[i] = i * delta;
		px[i] = r * cos(t[i]);
		py[i] = r * sin(t[i]);
		
	}
	for (int j = 0; j < 360; j++) {
		if (sqrt(pow((px[j] - qx[j]), 2) + pow((py[j] - qy[j]), 2) > max)) {
			max = sqrt(pow((px[j] - qx[j]), 2) + pow((py[j] - qy[j]), 2));
			x = qx[j];
			y = qy[j];
			angle = t[j];
			index = j;
		}
	}
	/*cout << "--------------------------------------------------" << endl;
	cout << "infinite norms:" << endl;
	cout << "max:" << fixed << setprecision(6) << max << "\tindex:" << fixed << setprecision(6) << index << "\tangle:" << fixed << setprecision(6) << angle << "\t(x,y) = (" << fixed << setprecision(6) << x << "," << fixed << setprecision(6) << y << ")" << endl;*/
}

double Lagrange(double *x, double *y, double in, int n) {
	double sum = 0.0;
	double temp;
	for (int i = 0; i < n; i++) {
		temp = y[i];
		for (int j = 0; j < n; j++) {
			if (i != j)
				temp *= ((in - x[j]) / (x[i] - x[j]));
		}
		sum += temp;
	}
	return sum;
}

double CubicLagrange(double *x, double *y, double in, int index, int n) {
	double sum = 0.0, temp;
	int k = index == 0 ? index + n : index;
	k += 2;
	int m;
	for (int i = (index - 1 + n) % n; i < k + 1; i++) {
		if (i == n) {
			i = i % n;
			k = k % n;
		}
		temp = y[i];
		m = index == 0 ? index + 2 + n : index + 2;
		for (int j = (index - 1 + n) % n; j < m + 1; j++) {
			if (j == n) {
				j %= n;
				m %= n;
			}
			if (i != j) {
				temp *= ((in - x[j]) / (x[i] - x[j]));
			}
		}
		sum += temp;
	}
	return sum;
}
void ProInterpolation(int n) {
	/*cout << "index\tpx[i]\tqx[i]\tpy[i]\tqy[i]" << endl;
	cout << "--------------------------------------------------" << endl;*/
	double delta = 2 * M_PI / 360.0;
	double *qx, *qy, *sample_x, *sample_y, *sample_t;
	qx = (double*)malloc(sizeof(double) * 360);
	qy = (double*)malloc(sizeof(double) * 360);
	sample_x = (double*)malloc(sizeof(double) * n);
	sample_y = (double*)malloc(sizeof(double) * n);
	sample_t = (double*)malloc(sizeof(double) * n);
	double delta_2 = 2 * M_PI / n;
	for (int i = 0; i < n; i++) {
		sample_t[i] = i * delta_2;
		sample_x[i] = r * cos(sample_t[i]);
		sample_y[i] = r * sin(sample_t[i]);
		
	}
	for (int i = 0; i < 361; i++) {
		if (LagrangeFlag) {//normal Lagrange
			qx[i] = Lagrange(sample_t, sample_x, i * delta, n);
			qy[i] = Lagrange(sample_t, sample_y, i * delta, n);
		}
		else {//cubic Lagrange
			qx[i] = CubicLagrange(sample_t, sample_x, i * delta, (i * n / 360), n);
			qy[i] = CubicLagrange(sample_t, sample_y, i * delta, (i * n / 360), n);
		}
		cout << fixed << setprecision(6) << qx[i] << " " << fixed << setprecision(6) << qy[i] << endl;
	}
	Norms_2(n, qx, qy);
	NormsInfi(n, qx, qy);
	cout << endl << endl;
}

void main() {
	LagrangeFlag = false;//true : normal Lagrange, false : cubic Lagrange	
	//ProInterpolation(12);
	ProInterpolation(24);
	//ProInterpolation(36);
	system("pause");
}
