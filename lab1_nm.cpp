#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "gnuplot.h"
#include "fun.h"
#include "spline.h"
#include "windows.h"
#include <chrono>
#define M_PI 3.14159265358979323846


using namespace std;
int main()
{
    Gnuplot plot;
    string command_out = "set xrange [0:10];"
        "set yrange [-5:5];"
        "plot 'output.csv' using 1:2 w l lt 1 lw 5 , log(x + 1)* sin(x)";
    string command_spline = "set xrange [0:10];"
        "set yrange [-5:5];"
        "plot 'spline.csv' using 1:2 w l lt 1 lw 3, log(x + 1)* sin(x)";
    string command_err = "set xrange [0:20];"
        "set yrange [0:5];"
        "plot 'errors.csv' using 1:2 w l";
    string command_Fou = "set xrange [-1:10];"
        "set yrange [-2:3];"
        "plot 'outF.csv' using 1:2 w l lt 1 lw 3 , log((3 * x + 2) / 2.0) * sin(3 * x / 2.0)";
    int n = 15;
    double A = 0, B = 3 * M_PI;
    double H = abs(A - B) / n;
    vector<pair<double, double>> table(n);
    ofstream out("output.csv");
    ofstream err("errors.csv");
    ofstream Fou("outF.csv");
    ofstream spline("spline.csv");

    for (int i = 0; i < n; i++)
    {
        table[i] = { i * H,f(i * H) };
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
            H = abs(A - B) / n;
            table2[i] = { i * H,f(i * H) };
            out << table2[i].first << " " << table2[i].second << endl;
        }
        //plot(command_out);
        //Sleep(3000);
        ofstream file("output.csv");
        double error = Max_Error(n, A, B, table2, "Lagrange");
        cout << std::fixed;
        cout << std::setprecision(8);
        int m = 9;
        if (n < 10) m = 11;
        cout << "n = " << n <<  ", " << setw(m) << error << endl;
        err << n << " " << error << endl;
        
    }
    //plot(command_err);
    //system("pause");

    //Задача 2
    //Задание 1
    cout << endl << "Chebyshev - Lagrange:" << endl << endl;
    vector<pair<double, double>> table_Cheb(n);
    for (int i = 0; i < n; i++)
    {
        table_Cheb[i].first = Chebyshev(i, n, A, B);
        table_Cheb[i].second = f(table_Cheb[i].first);
    }
    //задание 2,4
    cout << "ERROR Chebyshev: " << Max_Error(n, A, B, table_Cheb, "Lagrange") << endl;

    auto start = chrono::high_resolution_clock::now();
    cout << "ERROR Lagrange:  " << Max_Error(n, A, B, table, "Lagrange") << endl;
    auto end = chrono::high_resolution_clock::now();
    double seconds = chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1e3;
    printf("Lagrange worked for %.3f seconds\n", seconds);

    start = chrono::high_resolution_clock::now();
    Max_Error(n, A, B, table, "Newton");
    end = chrono::high_resolution_clock::now();
    seconds = chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1e3;
    printf("Newton worked for %.3f seconds\n", seconds);

    //Задача 4
    //Заданаие 1
    const int N = 15;
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

    //задание 6
    //vector<pair<double, double>> table{ { 1,2 }, { 2,3 }, { 4,1 }, { 7,4 } };
    n = table.size() - 1;
    vector<double> h(n + 1);
    vector<double> a(n);
    vector<double> b(n);
    vector<double> c(n + 1, { 0 });
    vector<double> d(n);

    vector<double> Alpha(n + 1, { 0 });
    vector<double> Betta(n + 1, { 1 });
    vector<double> Gamma(n + 1, { 0 });
    vector<double> Omega(n + 1, { 0 });

    inicilization_h(n, table, h);
    inicilization_ABGO(n, h, table, Alpha, Betta, Gamma, Omega);

    matrix_solution(n, c, Alpha, Betta, Gamma, Omega);

    inicilization_a(n, table, a);
    inicilization_b(n, table, b, h, c);
    inicilization_d(n, table, d, h, c);

    for (int i = 0; i < table.size() - 1; i++)
    {
        for (int j = 0; j < 1e3; j++)
        {
            H = table[i].first + j * (table[i + 1].first - table[i].first) / 1e3;
            spline << H <<  " " << ff(a[i], b[i], c[i], d[i], H, table[i].first) << endl;
        }
    }
    plot(command_spline);
    system("pause");
    /*for (int i = 0; i < c.size() - 1; i++)
    {
        cout << i << " a " << a[i] << endl;
        cout << i << " b " << b[i] << endl;
        cout << i << " c " << c[i] << endl;
        cout << i << " d " << d[i] << endl;
        cout << "__________________________" << endl;
    }*/
    spline.close();
    out.close();
    err.close();
    return 0;
}