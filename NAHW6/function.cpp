#define _CRT_SECURE_NO_WARNINGS
#include<iomanip>
#include<math.h>
#include <iostream>
#include "function.h"

using namespace std;

void gen_mtx(double **A, int n){
	srand(0);
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++) {
			A[i][j] = (double)(rand() % 50);;
			A[j][i] = A[i][j];
		}
	}
}

void print_mtx(double **A, int n){
	for (int i = 0; i < n; i++) {
		
		for (int j = 0; j < n; j++){
			cout << fixed << setprecision(8) << A[i][j]<<" ";
			
		}
		cout << endl;
	}
	cout   << "--------------------------" << endl<< endl;
}

void copy_mtx(double **A, double **cA, int n){
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			cA[i][j] = A[i][j];
}

void make_identity_mtx(double **A, int n){
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = (i == j) ? 1 : 0;
}

void make_rotate_mtx(double **A, int p, int q, double c, double s, int n){
	make_identity_mtx(A, n);
	A[p][p] = c;	A[p][q] = -s;
	A[q][p] = s;	A[q][q] = c;
}

double **alloc_mtx(int n){
	double **A;
	A = (double **)malloc(sizeof(double*) * n);
	for (int i = 0; i < n; i++)
		A[i] = (double*)malloc(sizeof(double) * n);
	return A;
}

//A = R * A
void pre_mtx_mult(double **R, double **A, int n){
	double **T;//temp mtx
	T = alloc_mtx(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			T[i][j] = 0.0;
			for (int k = 0; k < n; k++)
				T[i][j] += R[i][k] * A[k][j];
		}
	}
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = T[i][j];
}

//A = A * R
void post_mtx_mult(double **R, double **A, int n){
	double **T;//temp mtx
	T = alloc_mtx(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			T[i][j] = 0.0;
			for (int k = 0; k < n; k++)
				T[i][j] += A[i][k] * R[k][j];
		}
	}
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = T[i][j];
}

void transpose_mtx(double **A, int n){
	double temp;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			temp = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = temp;
		}
	}
}

double max_off_diag_entry(double **A, int *p, int *q, double* offdiag_sum, int n){
	//initial
	double max = A[0][1];
	*p = 0;
	*q = 1;
	*offdiag_sum = 0;
	//因為symmetric，找一半就好
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			*offdiag_sum += fabs(A[i][j]);//sigma | Aij |
			if (fabs(A[i][j]) > fabs(max)) {
				max = A[i][j];
				*p = i;
				*q = j;
			}
		}
	}
	return max;
}

double inner_product(double *a, double *b, int n){
	double sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += a[i] * b[i];
	return sum;
}

double *alloc_vec(int n){
	double *vector;
	vector = (double*)malloc(sizeof(double)*n);
	return vector;
}

//a(nx1) = A(nxn) * b(nx1)
void mtx_vec_mult(double *a, double **A, double *b, int n){
	for (int i = 0; i < n; i++) {
		a[i] = 0.0;
		for (int j = 0; j < n; j++)
			a[i] += A[i][j] * b[j];
	}
}

double vec_norm2(double *a, int n){
	double sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += a[i] * a[i];
	return sqrt(sum);
}

int jacobian_method(double **A, double **P, int flag, int n){
	double **R = alloc_mtx(n);				//rotation mtx
	int stp_cnt = 0;						//step count
	int p, q;								//index
	double A_pq, c, s, theta, offdiag_sum;
	double *offdiag_array;
	offdiag_array = (flag == 0) ? alloc_vec(50) : alloc_vec(1000);
	char str[20] = "OffDiagFile";
	char item1[11] = "iterations";
	char item2[10] = "offDiag";
	make_identity_mtx(P, n);				//set P to identity mtx	
	A_pq = max_off_diag_entry(A, &p, &q, &offdiag_sum, n);
	/* B = R.inverse * A * R = R.transpose * A * R
	 * after this transformation,
	 * B_qp = B_pq
	 * B_pq = (A_qq - A_pp) * sin(theta) * cos(theta) + A_pq( cos^2(theta) - sin^2(theta) )

	 * To elimate the entry B_pq -> find proper theta
	 * theta = (1/2) * atan(2 * A_pq / (A_pp - A_qq))
	 */
	while (fabs(A_pq) > EPSILON && stp_cnt < MAXSTEP){
		if ((A[p][p] - A[q][q]) == 0.0)
			theta = (A_pq > 0.0) ? PI / 2.0 : 3.0 * PI / 2.0;
		else
			theta = atan(2 * A_pq / (A[p][p] - A[q][q]));
		theta *= 0.5;
		c = cos(theta);
		s = sin(theta);
		make_rotate_mtx(R, p, q, c, s, n);
		post_mtx_mult(R, P, n);						//P = P * R : eigenvector
		post_mtx_mult(R, A, n);						//A = A * R
		transpose_mtx(R, n);						//R = R^t
		pre_mtx_mult(R, A, n);						//A = R^t * A * R;
		A[p][q] = A[q][p] = 0;						//enforce A[p][q] = A[q][p] = 0
		A_pq = max_off_diag_entry(A, &p, &q, &offdiag_sum, n);	//find next A_pq
		offdiag_array[stp_cnt] = offdiag_sum;
		stp_cnt++;
		if (flag == 0) {
			cout << "iteration= " << stp_cnt << ", offDiag[" << stp_cnt << "] = " << fixed << setprecision(3) << offdiag_sum << ", A[][]=" << endl;
			print_mtx(A, n);
		}
		else if (flag == 1) {
			if (stp_cnt % 20 == 0)
				cout << "iteration= " << stp_cnt << ", offDiag[" << stp_cnt << "] = " << fixed << setprecision(3) << offdiag_sum << endl;
		}

	}
	if (flag == 1)
		create_marks_csv(str, item1, item2, offdiag_array, stp_cnt - 1);
	return (stp_cnt - 1);
}

void extract_eigen(double **A, double **P, double *v, int n) {
	transpose_mtx(P, n);		//eigenvector
	for (int i = 0; i < n; i++)
		v[i] = A[i][i];			//eigenvalue
}

double cal_norm(double **A, double *vector, double value, int n) {
	double *t = alloc_vec(n);
	double norm2;
	//set t = A*v
	for (int i = 0; i < n; i++) {
		t[i] = 0.0;
		for (int j = 0; j < n; j++)
			t[i] += A[i][j] * vector[j];
	}
	//set t = A*v - lamda*v
	for (int i = 0; i < n; i++)
		t[i] = t[i] - value * vector[i];
	norm2 = vec_norm2(t, n);
	return norm2;
}

void eigen_vector_inner_product(double **A, int n) {
	double * sum = alloc_vec(n - 1);
	cout << "--------------------------" << endl;
	cout << "eigen vector inner product: " << endl;
	for (int i = 0; i < n - 1; i++) {
		sum[i] = 0.0;
		for (int j = 0; j < n; j++) {
			sum[i] += A[i][j] * A[i + 1][j];
		}
	}
	for (int i = 0; i < n - 1; i++) {
		cout << "( ";
		for (int j = 0; j < n - 1; j++) {
			cout << fixed << setprecision(3) << A[i][j]<<" ";
		}
		cout << fixed << setprecision(3) << A[i][n - 1] << ")。( ";
		for (int j = 0; j < n - 1; j++) {
			cout << fixed << setprecision(3) << A[i + 1][j] << " , ";
		}
		cout << fixed << setprecision(3) << A[i + 1][n - 1] << " ) = " << fixed << setprecision(3) << sum[i] << endl;
	}
	cout << "--------------------------" << endl;
	return;
}

void create_marks_csv(char *filename, char *item1, char *item2, double *offDiag, int n) {
	FILE *fp;
	filename = strcat(filename, ".csv");
	fp = fopen(filename, "w+");
	fprintf(fp, item1);
	fprintf(fp, " ,");
	fprintf(fp, item2);
	for (int i = 0; i < n; i++) {
		fprintf(fp, "\n%d,%lf", i + 1, offDiag[i]);
	}
	fclose(fp);
}
