#include <iostream> 
#include <iomanip> 
#include <math.h> 
#define M_PI 3.1415926535
using namespace std;

double f(double x)
{
	return sin(sqrt(x) + 1)* exp(sqrt(x)) / (sqrt(x));
}

double integral(double a, double b) {
	double n = 1e5;
	double step = (b - a) / n;
	double ans = 0.0;  
	for (int i = 0; i < n; i++) {
		ans += f(a + (i + 0.5) * step) * step; 
	}
	return ans;
}

int main()
{
	cout.precision(7);
	double x;
	cout << integral(M_PI, 3* M_PI);
	return 0;
}