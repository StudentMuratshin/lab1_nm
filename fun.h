#pragma once
#define M_PI 3.14159265358979323846
double f(double x)
{
    return log(x + 1) * sin(x);
}

double ft(double t)
{
    return log((3 * t + 2) / 2.0) * sin(3 * t / 2.0);
}


double ak(const int n, int k, const vector<pair<double, double>>& T)
{
    double sum = 0;
    if (k == 0)
    {
        for (int j = 0; j <= 2 * n; j++)
        {
            sum += ft(T[j].first);
        }
        return sum / (2 * n + 1);
    }

    for (int j = 0; j <= 2 * n; j++)
    {
        sum += ft(T[j].first) * cos((2 * M_PI * k * j) / (2 * n + 1));
    }
    return 2 * sum / (2 * n + 1);
}

double bk(const int n, int k, const vector<pair<double, double>>& T)
{
    if (k == 0) return 0;
    double sum = 0;
    for (int j = 0; j <= 2 * n; j++)
    {
        sum += ft(T[j].first) * sin((2 * M_PI * k * j) / (2 * n + 1));
    }
    return 2 * sum / (2 * n + 1);
}

double Fourier(int n, double t, const vector<pair<double, double>>& T)
{
    double sum = 0;
    for (int k = 0; k <= n; k++)
    {
        /*cout << std::setprecision(20);
        cout << ak(n, k, t) << " " << bk(n, k, t) << endl;
        cout << cos(k * t) << " " << sin(k * t) << endl;
        cout << ak(n, k, t) * cos(k * t) << " " << bk(n, k, t) * sin(k * t) << endl;*/
        sum += ak(n, k, T) * cos(k * t) + bk(n, k, T) * sin(k * t);
        //cout << sum << endl << endl;
    }
    return sum;
}

struct opt_dbl {
    double val = 0.0;
    bool has_value = false;
};

vector<opt_dbl> u_map;
double u(const vector<pair<double, double>>& T, int start_idx, int end_idx) {
    if (u_map[start_idx * T.size() + end_idx].has_value)
        return u_map[start_idx * T.size() + end_idx].val;
    if (start_idx == end_idx) {
        u_map[start_idx * T.size() + end_idx].val = T[start_idx].second;
        u_map[start_idx * T.size() + end_idx].has_value = true;
        return T[start_idx].second;
    }
    if (end_idx - start_idx == 1) {
        auto v = (T[end_idx].second - T[start_idx].second) / (T[end_idx].first - T[start_idx].first);
        u_map[start_idx * T.size() + end_idx].val = v;
        u_map[start_idx * T.size() + end_idx].has_value = true;
        return v;
    }
    auto v = (u(T, start_idx + 1, end_idx) - u(T, start_idx, end_idx - 1)) / (T[end_idx].first - T[start_idx].first);

    u_map[start_idx * T.size() + end_idx].val = v;
    u_map[start_idx * T.size() + end_idx].has_value = true;
    return v;
}

double Nn(const vector<pair<double, double>>& T, double x) {
    u_map = vector<opt_dbl>(T.size() * T.size(), { 0,false });
    double sum = u(T, 0, 0);
    for (int i = 1; i < T.size(); i++) {
        double Ui = u(T, 0, i);
        for (int j = 0; j < i; j++) {
            Ui *= (x - T[j].first);
        }
        sum += Ui;
    }
    return sum;
}

double Lagrange(double x0, int n, const vector<pair<double, double>>& T)
{
    double S = 0;
    for (int i = 0; i < n; i++)
    {
        double P = 1;
        for (int j = 0; j < n; j++)
        {
            if (i == j) continue;
            P *= (x0 - T[j].first) / (T[i].first - T[j].first);
        }
        S += P * T[i].second;
    }
    return S;
}

double Chebyshev(int i, int n, double a, double b)
{
    return (b + a) / 2 + (b - a) / 2 * cos(3.14 * (2 * i + 1) / (2 * n));
}

double Max_Error(int n, double a, double b, const vector<pair<double, double>>& T, string Name)
{
    const int k = 1e5;
    double h = abs(a - b) / k;
    double error = 0;
    double fT, fx, L, Newton, four;
    if (Name == "Lagrange")
    {
        for (int i = 0; i <= k; i++)
        {
            fx = f(i * h);
            L = Lagrange(i * h, T.size(), T);
            double e = abs(fx - L);
            if (e > error) error = e;
        }
    }
    else if (Name == "Newton")
    {
        for (int i = 0; i <= k; i++)
        {
            fx = f(i * h);
            Newton = Nn(T, i * h);
            double e = abs(fx - Newton);
            if (e > error) error = e;
        }
    }
    else
    {
        for (int i = 0; i <= k; i++)
        {
            fT = ft(i * h);
            four = Fourier(n, i * h, T);
            double e = abs(fT - four);
            if (e > error) error = e;
        }
    }
    return error;
}