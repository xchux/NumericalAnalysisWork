#define _CRT_SECURE_NO_WARNIGS

#include<iomanip>
#include<math.h>
#include <iostream>
using namespace std;

//allocate
/*--------------------------------------------------------
* Procedure to allocate an n X n matrix and return the
* pointer to the matrix.
*/
double **alloc_mtx(int n)
{
	double  **B;
	int     i;

	B = (double **)malloc(sizeof(double *)*n);
	for (i = 0; i < n; i++)
		B[i] = (double *)malloc(sizeof(double)*n);
	return (B);
}
/*-----------------------------------------------------
* Procedure to allocate space for an n-dimensional
* vector. Return the pointer to the vector.
*/
double *alloc_vec(int n)
{
	double *t;
	t = (double *)malloc(sizeof(double)*n);
	return(t);
}

//initialize
/*------------------------------------------------------
* Create a Hilbert linear system.
*   A: the coef. mtx,
*   b: the right hand side.
*   n: dimension of matrix.
*/
void Hilbert_linear_system(double **A, double *b, int n)
{
	int  i, j;
	for (i = 0; i < n; i++) {
		b[i] = 0.0;
		for (j = 0; j < n; j++) {
			A[i][j] = 1.0 / (i + j + 1.0);
			b[i] += A[i][j];
		}
	}
	return;
}
/*------------------------------------------------------
* Create a Householder vector for a matrix. This vector
* is used to create a Householder matrix H such that
* H*A makes the entries below A[i][j] of A[][j] be 0.
*/
void create_H_vec(double **A, int i, int j, double *v, int n)
{
	int      k;
	double   norm2;

	for (k = 0; k < i; k++)
		v[k] = 0.0;
	norm2 = 0.0;
	for (k = i; k < n; k++)
		norm2 += A[k][j] * A[k][j];//calaulate ||A.j||
	norm2 = sqrt(norm2);

	if (A[i][j] >= 0.0)
		v[i] = A[i][j] + norm2;
	else
		v[i] = A[i][j] - norm2;

	for (k = i + 1; k < n; k++)
		v[k] = A[k][j];
}
/*void testData(double** A, double* b, int n) {

	A[0][0] = 0;	A[0][1] = 3;	A[0][2] = 3;	A[0][3] = 2;
	A[1][0] = 2;	A[1][1] = 2;	A[1][2] = 4;	A[1][3] = 4;
	A[2][0] = 8;	A[2][1] = 4;	A[2][2] = 4;	A[2][3] = 0;
	A[3][0] = 1;	A[3][1] = 2;	A[3][2] = 4;	A[3][3] = 3;
	b[0] = 219;
	b[1] = 326;
	b[2] = 372;
	b[3] = 276;
}*/

//print
/*------------------------------------------------------
* Procedure to print out a matrix.
/*void print_mtx(double **A, int n)
{
	int i, j;

	for (i = 0; i < n; i++) {
		cout << endl;
		for (j = 0; j < n; j++)
			cout <<  A[i][j]<<" ";
	}
	cout << endl;
	cout << "--------------------------------" << endl;

}*/
/*--------------------------------------------------------
* Procedure to print out a vector.
*/
void print_vec(double *x, int n)
{
	for (int i = 0; i < n; i++) {
		fprintf(stderr, " %lf", x[i]);
		if (i % 4 == 3)
			cout << endl;
	}
}

//matrix & vector calculate
/*---------------------------------------------------------
* Compute a = A*b, where a and b are vectors and A an n x n
* matrix.
*/
void mtx_vec_mult(double *a, double **A, double *b, int n)
{
	for (int i = 0; i < n; i++) {
		a[i] = 0.0;
		for (int j = 0; j < n; j++)
			a[i] += A[i][j] * b[j];
	}
}
/*--------------------------------------------------
* Compute inner product <a, b>, where a and b are
* n-dimensional vectors.
*/
double  inner_product(double *a, double *b, int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += a[i] * b[i];
	return (sum);
}
/*----------------------------------------------------------------
* Procedure to solve a lower triangular system
*    L*x = b.
*/
void back_substitute(double **U, double *x, double *b, int n)
{
	int i, j;

	for (i = n - 1; i >= 0; i--) {
		x[i] = b[i] / U[i][i];
		for (j = i - 1; j >= 0; j--)
			b[j] -= U[j][i] * x[i];
	}
}
//calculate 2 norm
double norm_2(double a, double b) {
	return sqrt((a - b)*(a - b));
}

//algoritm
/*
* Gaussian_elmination.
*   A: the coef. mtx.
*   b: the right hand side.
*   x: the unknown vector.
*   n: dimension of matrix.
*/
void gauss_elm(double** A, double* b, double* x, int n) {
	int k;
	double ratio, temp;
	for (int i = 0; i < (n - 1); i++) {
		//Partial pivoting
		k = i; //k record the index of pivot in which row
		temp = fabs(A[i][i]);
		for (int j = i + 1; j < n; j++) {
			if (fabs(A[j][i]) > temp) {
				k = j;
				temp = fabs(A[j][i]);
			}
		}

		//swapping if necessary, make pivot to the top
		if (i != k) {
			for (int j = i; j < n; j++) {
				temp = A[i][j];
				A[i][j] = A[k][j];
				A[k][j] = temp;
			}

			temp = b[i];
			b[i] = b[k];
			b[k] = temp;
		}

		//Elimination of rows below
		for (int j = i + 1; j < n; j++) {//each row
			if (A[j][i] == 0.0) continue;
			ratio = A[j][i] / A[i][i];
			for (k = i; k < n; k++)//each col
				A[j][k] -= ratio * A[i][k];
			b[j] -= ratio * b[i];
		}//end_for(j)
	}//end_for(i)
	cout << endl;

	
}
/*------------------------------------------------------
* Procedure to reflect matrix A into an upper triangular
* matrix by using a sequence of Householder matrices.
*    H(n-2)H(n-3)...H(0)Ax = H(n-2)H(n-3)...H(0)b.
* ---> Rx = b'.
*/
void QR_reflect(double **A, double *b, int n)
{
	int     i, j, k;
	double  *v, vv, vx;
	v = alloc_vec(n);
	//Eliminate each column to make A[][] into U[][].
	for (j = 0; j < n - 1; j++) {
		// Create vector v[] = A.j + ||A.j|| * e1.
		create_H_vec(A, j, j, v, n);
		vv = inner_product(v, v, n);
		// Update each column A.k, j<=k<=n-1.
		for (k = j; k < n; k++) {
			// Compute vx=<v, A.k>
			vx = 0.0;
			for (i = j; i < n; i++)
				vx += A[i][k] * v[i];
			// A.k = A.k - 2(<v, A.k>/<v,v>)v ;
			for (i = j; i < n; i++)
				A[i][k] -= 2.0*(vx / vv)*v[i];
		}
		// Update b. b = b - 2(<v,b>/<v,v>)v ;
		vx = inner_product(v, b, n);
		for (i = j; i < n; i++)
			b[i] -= 2.0*(vx / vv)*v[i];
	}
	//print_mtx(A, n);
}

//main branch
/*---------------------------------------------------
* Procedure to solve a linear system by using QR-
* decompsotion.
* Algm:
*    1. convert Ax=b into an upper triangular system.
*    2. Solve the upper trianglular system by using
*       backward substitution.
*/
void QR_solver(int start, int end)
{
	double **A;
	double *b, *x, *r;
	for (int N = start; N <= end; N += 2) {//loop 8~16 by 2
		A = alloc_mtx(N);
		b = alloc_vec(N);
		x = alloc_vec(N);
		r = alloc_vec(N);
		//testData(A, b, N);
		Hilbert_linear_system(A, b, N);//recover A & b
		QR_reflect(A, b, N);// Reflect A into an upper triangular matrix.
		back_substitute(A, x, b, N);// Solve the upper triangular system by using backward substitution.
		Hilbert_linear_system(A, b, N);//recover A & b

		cout << "N = " << N << " ,x[] = " << endl;
		print_vec(x, N);
		cout << endl << endl;
		cout << "e[] = " << endl;
		for (int i = 0; i < N; i++) {
			cout << fixed << setprecision(6) << norm_2((1.0 - x[i]), x[i]) << " ";
			if (i % 8 == 7)
				cout << endl;
		}
		cout << endl << endl;
		cout << "r[] = " << endl;
		mtx_vec_mult(r, A, b, N);
		for (int i = 0; i < N; i++) {
			cout << fixed << setprecision(6) << norm_2(r[i], x[i]) << " ";
			if (i % 8 == 7)
				cout << endl;
		}
		cout << endl;
		cout << "--------------------------------" << endl;
		free(A); free(b); free(x); free(r);
	}
}
void main_gauss(int start, int end) {
	double **A;
	double *b, *x, *r;
	for (int N = start; N <= end; N += 2) {
		A = alloc_mtx(N);
		b = alloc_vec(N);
		x = alloc_vec(N);
		r = alloc_vec(N);
		Hilbert_linear_system(A, b, N);//recover A & b
									   //testData(A, b, N);
		gauss_elm(A, b, x, N);
		back_substitute(A, x, b, n);
		cout << "N = " << N << " ,x[] = " << endl;
		print_vec(x, N);
		cout << endl << endl;

		cout << "e[] = " << endl;
		for (int i = 0; i < N; i++) {
			cout << fixed << setprecision(6) << norm_2((1.0 - x[i]), x[i]) << " ";
			if (i % 8 == 7)
				cout << endl;
		}
		Hilbert_linear_system(A, b, N);//recover A & b
		cout << endl<<endl;
		cout << "r[] = "<<endl;
		mtx_vec_mult(r, A, b, N);
		for (int i = 0; i < N; i++) {
			cout << fixed << setprecision(6) << norm_2(r[i], x[i]) << " ";
			if (i % 8 == 7)
				cout << endl;
		}
		cout << endl;
		cout << "--------------------------------" << endl;
		free(A); free(b); free(x); free(r);
	}
}

void main() {
	int start = 8, end = 16;
	cout << "Gaussian_Elimination from " << start << " to " << end << endl;
	main_gauss(start, end);
	cout << endl;
	cout << "QR_Method from " << start << " to " << end << endl;
	QR_solver(start, end);
	system("pause");
}