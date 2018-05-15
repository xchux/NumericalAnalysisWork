#include<iomanip>
#include<math.h>
#include <iostream>
using namespace std;
#define start 2000
#define end 0


double T[21][21] = { 0 }, x[21] = { 0 }, y[21] = { 0 }, b[21] = {0};
double h = 0.05;
double sum = 1;

void SetBoard() {
	for (int i = 0; i < 21; i++)
	{
		T[i][0] = 20;
		T[i][20] = 20;
		T[20][i] = 20;
	}

}
double SOR_iteration(double w)
{
	double   error, oldT;
	error = 0;
	for (int i = 1; i < 20; i++)
	{
		for (int j = 1; j < 20; j++)
		{
			oldT = T[i][j];
			T[i][j] += w*(T[i + 1][j] + T[i][j + 1] + T[i - 1][j] + T[i][j - 1] - 4.0 * T[i][j]) / 4.0;
			//T[i][j] = T[i][j] + (w / 4.0)*(-s*h*h + (T[i - 1][j] + T[i][j - 1] + T[i][j + 1] + T[i + 1][j]) - 4.0*T[i][j]);
			if (i == j&&i == 10) {
				T[i][j] += start*0.0025;
			}
			if (error < fabs(oldT - T[i][j]))
				error = fabs(oldT - T[i][j]);
		}
	}
	for (int i = 1; i < 20; i++)
		T[0][i] += w*(T[0 + 1][i] + T[0][i + 1] + T[0][i - 1] + end - 4.0 * T[0][i]) / 4.0;
	return(error);
}
/*double gauss_seidel(double w)
{
	double  error =0, sum,oldT;
	for (int i = 0; i < 20; i++) {
		sum = b[i];
		for (int j = 0; j < 20; j++) {
		
			if (j != i) 
				sum-= T[i][j] * x[j];
			T[i][j] += w*(T[i + 1][j] + T[i][j + 1] + T[i - 1][j] + T[i][j - 1] - 4.0 * T[i][j]) / 4.0;
		}
		oldT = x[i];
		x[i] += w*sum / T[i][i];	
		if (error < fabs(oldT - x[i]))
			error = fabs(oldT - x[i]);
	}	
	for (int i = 1; i < 20; i++)
		T[0][i] += w*(T[0 + 1][i] + T[0][i + 1] + T[0][i - 1] + end - 4.0 * T[0][i]) / 4.0;
	
	return(error);
}*/
int main()
{
	int counter = 0;
	double w = 1.0, error = 1;

	T[10][10] = start;
	cout << "gauss_seidel" << endl;
	SetBoard();

	for (int i = 0; i < 21; i++) {
		b[i] = 0.0;
		for (int j = 0; j < 21; j++) {
			b[i] += T[i][j];
		}
		//cout << b[i] << " ";
	}
	cout << endl;
	while (error > 0.000001)
	{
		error = SOR_iteration(1.5);

	}
	for (int i = 0; i < 21; i++)
	{
		for (int j = 0; j < 21; j++)
		{
			cout << fixed << setprecision(3) << T[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	////////////	////////////	////////////	////////////	////////////
	cout << "SOR_method" << endl;
	w = 1.5;
	error = 1;
	SetBoard();
	//Guess an initial solution.
	for (int i = 0; i < 20; i++)
		y[i] = x[i] = 0.0;

	while (error > 0.000001)
	{
		error = SOR_iteration(w);

	}
	//cout << "20 20 "<<	1.0/20 << endl;
	for (int i = 0; i < 21; i++)
	{
		for (int j = 0; j < 21; j++)
		{
			cout << fixed << setprecision(3) << T[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;

	system("pause");
}