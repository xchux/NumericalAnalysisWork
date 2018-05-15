#include <stdio.h>
#include <math.h>
#include<iomanip>
#include <iostream>
using namespace std;

double func(double x)
{
	//double x2 = x*x;
	//return (x2/2 +x+ 1- exp(x));//->1
	double x2 = x*x,x3= x*x*x,x4= x*x*x*x;
	return (x4 -2*x3+2*x2-2*x+1);//->2


}

// funcd(x) = func'(x)
double funcd(double x)
{
	//double x2 = x*x;
	//return ( x+1 - exp(x));//->1
	double x2 = x*x, x3 = x*x*x;
	return (4*x3-6*x2+4*x-2);//->2

}

// -------------------------------------------------------
double NewtonRoot(double x0,              /*   初點  */
	double(*fx)(double),   /* 適應函式*/
	double(*fd)(double),   /* 微分函式*/
	double eps            /* 容許誤差*/)        
{
	double x = x0,ex;
	int i = 1;
	do {
		x0 = x;
		//x = x0 - fx(x0) / fd(x0);//->1 & 2(a)
		//x = x0 - 2*fx(x0) / fd(x0);//->2(b)
		x = x0 - sqrt(2) * fx(x0) / fd(x0);//->2(c)
		
		ex = fabs(x - x0);
		//cout << i << "  "<<fixed<<setprecision(9)<< x << "  "<< fixed << setprecision(9)<< fx(x) <<"  "<< fixed << setprecision(9) << ex << endl;//->1 & 2
		cout << fixed << setprecision(3) <<i<<"  "<<fixed<<setprecision(9)<< x << "  " <<fixed<<setprecision(9)<< 1-x <<"  "<<fixed<<setprecision(9)<< x-x0 << endl;//->2(d)
		i++;
	} while (ex>eps);
	//cout << x << " "<< ex << endl;
	return x;

}

int main()
{
	const double eps = 1E-9;
	const int max_iterator = 100;
	double x0, x;

	//x0 = 1;//->1
	x0 = 11;//->2
	x = NewtonRoot(x0, func, funcd, eps);
	//cout << "1" << " " << x << endl;
	

	system("pause");
	return 0;
}