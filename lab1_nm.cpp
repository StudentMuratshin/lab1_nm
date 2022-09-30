#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "gnuplot.h"
#include "windows.h"


using namespace std;

double f(double x)
{
    return log(x + 1) * sin(x);
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
    vector<pair<double, double>> table(n);
    ofstream out("output.csv");
    ofstream err("errors.csv");

    for (int i = 0; i < n; i++)
    {
        table[i] = { i * h,f(i * h) };
    }

    //h = abs(a - b) / (div * n);

    //задание 2
    //for (int i = 0; i < div * n; i++)
    //{
    //    L = Lagrange(i * h, n, table);
    //    fx = f(i * h);
    //    cout << std::fixed;
    //    cout << std::setprecision(8);
    //    cout << i + 1 << ": " << fx << "  -  " << L << "  =  " << abs(fx - L) << endl << endl;
    //    out << i * h << " " << L << endl;
    //    if (i != 0) {
    //        plot(command_out);
    //        //system("pause");
    //        Sleep(500);
    //    }
    //}

    //задание 3
    for (n = 1; n <= 15; n++)
    {
        vector<pair<double, double>> table3(n + 1);
        const int k = 1e5;
        for (int i = 0; i <= n; i++)
        {
            h = abs(a - b) / n;
            for (int i = 0; i <= n; i++) {
                table3[i] = { i * h,f(i * h) };
            }
        }
        h = abs(a - b) / k;
        double error = 0;
        for (int i = 0; i <= k; i++)
        {
            fx = f(i * h);
            L = Lagrange(i * h, table3.size(), table3);
            double e = abs(fx - L);
            if (e > error) error = e;
        }
        cout << std::fixed;
        cout << std::setprecision(8);
        cout << "n = " << n << ", " << error << endl;
        err << n << " " << error << endl;
        
    }
    plot(command_err);
    system("pause");
    //plot("exit");
    out.close();
    err.close();

    //Задача 2
    //Задание 1
    //n = 15;
    //double x2[15];
    //double y2[15];
    //for (int i = 0; i < n; i++)
    //{
    //    x2[i] = Chebyshev(i, n, a, b);
    //    y2[i] = f(x2[i]);
    //}
    ////Задание 2,4
    //for (int i = 0; i < n; i++)
    //{
    //    L = Lagrange(i * h, n, x, y);
    //    L2 = Lagrange(i * h, n, x2, y2);
    //    fx = f(i * h);
    //    cout << abs(L2 - L) << endl;
    //}

    //// Задача 4
    //// Задание 1
    //int n0 = 3;
    //h = abs(a - b) / n0;
    //vector<pair<double, double>> table{ {1, 6}, {3, 24},{4, 45} };
    //for (int i = 0; i < n0; i++) {
    //    table[i] = { i * h,f(i * h) };
    //}
    //for (int i = 0; i < n0; i++) {
    //    double e = (f(i * h) - Nn(table, i * h));
    //}
    //cout << Nn(table, 5);
    return 0;
}