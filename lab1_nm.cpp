#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "gnuplot.h"
#include "windows.h"
#include <chrono>
#define M_PI 3.14159265358979323846


using namespace std;

double f(double x)
{
    return log(x + 1) * sin(x);
}

double ft(double t)
{
    return log((3 * t + 2) / 2.0) * sin(3 * t / 2.0);
}

double a0(const int n, const vector<pair<double, double>>& T)
{
    double sum = 0;
    for (int j = 0; j <= 2 * n; j++)
    {
        sum += ft(T[j].first);
    }
    return sum / (2 * n + 1);
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


int main()
{
    Gnuplot plot;
    string command_out = "set xrange [0:10];"
        "set yrange [-5:5];"
        "plot 'output.csv' using 1:2 w l lt 1 lw 5 , log(x + 1)* sin(x)";
    string command_err = "set xrange [0:20];"
        "set yrange [0:5];"
        "plot 'errors.csv' using 1:2 w l";
    string command_Fou = "set xrange [-1:10];"
        "set yrange [-2:3];"
        "plot 'outF.csv' using 1:2 w l lt 1 lw 3 , log((3 * x + 2) / 2.0) * sin(3 * x / 2.0)";
    const int n = 15;
    double a = 0, b = 3 * M_PI;
    double h = abs(a - b) / n;
    double L, L2, fx;
    double x[15], y[15];
    vector<pair<double, double>> table(n);
    ofstream out("output.csv");
    ofstream err("errors.csv");
    ofstream Fou("outF.csv");

    for (int i = 0; i < n; i++)
    {
        table[i] = { i * h,f(i * h) };
    }

    

    //задание 2
    //h = abs(a - b) / (div * n);
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

    //задание 2, 3
    for (int n = 1; n <= 15; n++)
    {
        vector<pair<double, double>> table2(n + 1);
        for (int i = 0; i <= n; i++)
        {
            h = abs(a - b) / n;
            table2[i] = { i * h,f(i * h) };
            //out << table2[i].first << " " << table2[i].second << endl;
        }
        //plot(command_out);
        //Sleep(3000);
        //ofstream file("output.csv");
        double error = Max_Error(n, a, b, table2, "Lagrange");
        cout << std::fixed;
        cout << std::setprecision(8);
        int m = 9;
        if (n < 10) m = 11;
        cout << "n = " << n <<  ", " << setw(m) << error << endl;
        err << n << " " << error << endl;
        
    }
    //plot(command_err);
    //system("pause");
    //plot("exit");

    //Задача 2
    //Задание 1
    cout << endl << "Chebyshev - Lagrange:" << endl << endl;
    vector<pair<double, double>> table_Cheb(n);
    for (int i = 0; i < n; i++)
    {
        table_Cheb[i].first = Chebyshev(i, n, a, b);
        table_Cheb[i].second = f(table_Cheb[i].first);
    }
    //задание 2,4
    //cout << "ERROR Lagrange:  " << Max_Error(a, b, table, "Lagrange") << endl; 
    cout << "ERROR Chebyshev: " << Max_Error(n, a, b, table_Cheb, "Lagrange") << endl;

    auto start = chrono::high_resolution_clock::now();
    cout << "ERROR Lagrange:  " << Max_Error(n, a, b, table, "Lagrange") << endl;
    auto end = chrono::high_resolution_clock::now();
    double seconds = chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1e3;
    printf("Lagrange worked for %.3f seconds\n", seconds);

    start = chrono::high_resolution_clock::now();
    Max_Error(n, a, b, table, "Newton");
    end = chrono::high_resolution_clock::now();
    seconds = chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1e3;
    printf("Newton worked for %.3f seconds\n", seconds);

    //Задача 4
    //Заданаие 1
    const int N = 100;
    vector<pair<double, double>> table_F;
    for (int i = 0; i <= 2 * N; i++)
    {
        table_F.push_back( { 2 * M_PI * i / (2 * N + 1), 0});
    }
    cout << endl;
    for (int i = 0; i <= 2 * N; i++)
    {
        table_F[i].second = Fourier(N, table_F[i].first, table_F); 
        cout << table_F[i].first << " " << table_F[i].second << endl << endl;
        Fou << table_F[i].first << " " << table_F[i].second << endl;
    }
    cout << Max_Error(N, 0, 2 * M_PI, table_F, "Fourier") << endl;
    plot(command_Fou);
    system("pause");
    out.close();
    err.close();
    return 0;
}