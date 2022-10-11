#pragma once


double ff(double& a, double& b, double& c, double& d, double x, double xi)
{
	return a + b * (x - xi) + c * (x - xi) * (x - xi) + d * (x - xi) * (x - xi) * (x - xi);
}

void inicilization_h(const int& n, const vector<pair<double, double>>& T, vector<double>& h)
{
	for (int i = 1; i <= n; i++)
	{
		h[i] = T[i].first - T[i - 1].first;
	}
}

void inicilization_a(const int& n, const vector<pair<double, double>>& T, vector<double>& a)
{
	for (int i = 0; i < n; i++)
	{
		a[i] = T[i].second;
	}
}

void inicilization_b(const int& n, const vector<pair<double, double>>& T, vector<double>& b, vector<double>& h, vector<double>& c)
{
	for (int i = 1; i <= n - 1; i++)
	{
		b[i - 1] = (T[i].second - T[i - 1].second) / h[i] - h[i] / 3 * (c[i] + 2 * c[i - 1]);
	}
	b[n - 1] = (T[n].second - T[n - 1].second) / h[n] - (2 * h[n] * c[n - 1]) / 3;
}

void inicilization_d(const int& n, const vector<pair<double, double>>& T, vector<double>& d, vector<double>& h, vector<double>& c)
{
	for (int i = 1; i < n; i++)
	{
		d[i - 1] = (c[i] - c[i - 1]) / (3 * h[i]);
	}
	d[n - 1] = (-c[n - 1] / (3 * h[n]));
}


void inicilization_ABGO(const int& n, vector<double>& h, const vector<pair<double, double>>& T, vector<double> &Alpha, vector<double>& Betta, vector<double>& Gamma, vector<double>& Omega)
{
	for (int i = 2; i <= n; i++)
	{
		Alpha[i - 1] = h[i - 1];
		Betta[i - 1] = 2 * (h[i - 1] + h[i]);
		Gamma[i - 1] = h[i];
		Omega[i - 1] = 3 * ((T[i].second - T[i - 1].second) / h[i] - (T[i - 1].second - T[i - 2].second) / h[i - 1]);
	}
}

void matrix_solution(const int& n, vector<double>& c, vector<double>& Alpha, vector<double>& Betta, vector<double>& Gamma, vector<double>& Omega)
{
	vector <double> P(n + 2, { 0 });
	vector <double> Q(n + 2, { 0 });

	for (int i = 1; i <= n; i++)
	{
		/*P[i] = -Gamma[i - 1] / (Alpha[i] * P[i - 1] - Betta[i - 1]);
		Q[i] = (Omega[i - 1] - Alpha[i - 1] * Q[i - 1]) / (Alpha[i - 1] * P[i - 1] - Betta[i - 1]);*/

		P[i + 1] = Gamma[i] / (Betta[i] - Alpha[i] * P[i]);
		Q[i + 1] = -(Omega[i] - Alpha[i] * Q[i]) / (Alpha[i] * P[i] - Betta[i]);
	}
	for (int i = n; i > 0; i--)
	{
		c[i - 1] = -P[i] * c[i] + Q[i];
	}
}