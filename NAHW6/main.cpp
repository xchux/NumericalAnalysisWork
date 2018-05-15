#define _CRT_SECURE_NO_WARNINGS
#include<iomanip>
#include<math.h>
#include <iostream>

#include "function.h"

using namespace std;

void promblem1() {
	int n = 4;
	double **A, **P, **originA;
	double *v;					//eigenvalue array
	A = alloc_mtx(n);
	originA = alloc_mtx(n);
	P = alloc_mtx(n);
	v = alloc_vec(n);

	gen_mtx(A, n);				//initialize A, symmetric
	copy_mtx(A, originA, n);	//store original A in originA

	jacobian_method(A, P, 0, n);
	extract_eigen(A, P, v, n);

	cout << "eigenvalue, eigenvector, result" << endl;
	cout << "--------------------------"<< endl;
	for (int i = 0; i < n; i++) {
		cout << "(", v[i];
		for (int j = 0; j < n - 1; j++) {
			cout << P[i][j] << " ";
		}
		double norm2 = cal_norm(originA, P[i], v[i], n);
		cout << fixed << setprecision(3) << P[i][n - 1] << " ) " << fixed << setprecision(3) << norm2 << endl;
	}
	cout << endl;
	eigen_vector_inner_product(P, n);

	return;
}

void promblem2() {
	int n = 20;
	double **A, **P, **originA;
	double *v;					//eigenvalue array
	A = alloc_mtx(n);
	originA = alloc_mtx(n);
	P = alloc_mtx(n);
	v = alloc_vec(n);

	gen_mtx(A, n);				//initialize A, symmetric
	copy_mtx(A, originA, n);	//store original A in originA

	jacobian_method(A, P, 1, n);
	extract_eigen(A, P, v, n);
	cout << "--------------------------" << endl;
	cout << "eigenvalue : []= "<<endl;
	for (int i = 0; i < n - 1; i++){
		cout << fixed << setprecision(6) << v[i] << " ,";
		if(i%6==5)
			cout <<endl;
	}
	cout << fixed << setprecision(6) << v[n - 1]  << endl;
	return;
}

void promblem3() {
	double iteration[20] = { 0 };
	char str[20] = "IterationsFile";
	char item1[10] = "N";
	char item2[15] = "K";
	for (int n = 2; n < 20; n++) {
		double **A, **P, **originA;
		double *v;					//eigenvalue array
		A = alloc_mtx(n + 1);
		originA = alloc_mtx(n + 1);
		P = alloc_mtx(n + 1);
		v = alloc_vec(n + 1);

		gen_mtx(A, n + 1);				//initialize A, symmetric
		copy_mtx(A, originA, n + 1);	//store original A in originA

		iteration[n] = jacobian_method(A, P, 2, n + 1);
		for (int i = 0; i < n; i++) {
			cout << iteration[i] << " ";
		}
		cout << endl;
		extract_eigen(A, P, v, n + 1);
	}

	create_marks_csv(str, item1, item2, iteration, 20);

	return;
}

void main() {
	promblem1();
	promblem2();
	promblem3();

	system("pause");
}