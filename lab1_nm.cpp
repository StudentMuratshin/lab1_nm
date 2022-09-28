#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "gnuplot.h"
#include "windows.h"


using namespace std;

double f(double x)
{
	return log(x + 1) * sin(x);
}

double Lagrange(double x0, int n, double* x, double* y)
{
	double S = 0;
	for (int i = 0; i <= n; i++)
	{
		double P = 1;
		for (int j = 0; j <= n; j++)
		{
			if (i == j) continue;
			P *= (x0 - x[j]) / (x[i] - x[j]);
		}
		S += P * y[i];
	}
	return S;
}

double Chebyshev(int i, int n, double a, double b)
{
	return (b + a) / 2 + (b - a) / 2 * cos(3.14 * (2 * i + 1) / (2 * n));
}

int main()
{
	Gnuplot plot;
	string command_out = "set xrange [-1:10];"
		"set yrange [-5:5];"
		"plot 'output.csv' using 1:2 w l lt 1 lw 5 , log(x + 1)* sin(x)";
	string command_err = "set xrange [0:20];"
		"set yrange [0:77];"
		"plot 'errors.csv' using 1:2 w l";
	int n = 15, div = 2;
	double a = 0, b = 3 * 3.14;
	double h = abs(a - b) / n;
	double L, L2, fx;
	double x[15], y[15];
	ofstream out("output.csv");
	ofstream err("errors.csv");

	for (int i = 0; i < n; i++)
	{
		x[i] = i * h;
		y[i] = f(x[i]);
	}

	h = abs(a - b) / (div * n);
	
	//задание 2
	for (int i = 0; i < div * n; i++)
	{
		L = Lagrange(i * h, n, x, y);
		fx = f(i * h);
		cout << std::fixed;
		cout << std::setprecision(8);
		cout << i+1 << ": " << fx << "  -  " << L << "  =  " << abs(fx - L) << endl << endl;
		out << i * h << " " << L << endl;
		if (i != 0) {
			//plot(command_out);
			//system("pause");
			//Sleep(500);
		}
	}

	//задание 3
	double prev_error;
	int n0, key = 0;
	for (n = 1; n <= 15; n++)
	{
		double error = 0;
		cout << "n:__" << n << "_______________________________________" << endl << endl;
		h = abs(a - b) / n;
		for (int i = 0; i < n; i++)
		{
			L = Lagrange(i * h, n, x, y);
			fx = f(i * h);
			cout << std::fixed;
			cout << std::setprecision(8);
			cout << "i: " << i << "    " << fx << "  -  " << L << "  =  "  << abs(fx - L) << endl << endl;
			
			if (abs(fx - L) > error)
			{
				error = abs(fx - L);
			}
		}
		err << n << " " << error << endl;
		if (key == 0 && error != 0) { prev_error = error; key = 1; }

		if (error != 0 && prev_error > error)
		{
			prev_error = error;
			n0 = n;
			
		}
	}
	cout << "Error: " << prev_error << "   Best n0: " << n0;
	plot(command_err);
	system("pause");
	//plot("exit");
	out.close();
	err.close();

	//Задача 2
	//Задание 1
	n = 15;
	double x2[15];
	double y2[15];
	for (int i = 0; i < n; i++)
	{
		x2[i] = Chebyshev(i, n, a, b);
		y2[i] = f(x2[i]);
	}
	//Задание 2
	for (int i=0; i < n; i++)
	{
		L = Lagrange(i * h, n, x, y);
		L2 = Lagrange(i * h, n, x2, y2);
		fx = f(i * h);
		cout << abs(L2 - L) << endl;
	}
	return 0;
}