#include <iostream>
#include <fstream>
#include "math.h"
#include <cmath>
#include<array>

using namespace std;

//Variables
double t_0 = 0.;
double T = 5.;
double x_0 = 0.1;
double eps = 0.000001;
int p = 4;

const int J = 1000;
int M = 10;

double tau = (T - t_0) / M;
double t[J + 1];
double x[J + 1];

//array<double, J> x {};
//array<double, J> t {};

double w_1, w_2, w_3, w_4;
double x_main, x_new, x_curr;
double tau_opt, tau_old;

double a(double t) {
    //return 1;
    return -1;
}

double b(double t) {
    //return -2 * exp(t);
    return 2 * t;
}

double c(double t) {
    //return exp(2 * t) + exp(t);
    return 5 - t * t;
}

double f(double x, double t, double a(double), double b(double), double c(double)) {
    return a(t) * x * x + b(t) * x + c(t);
}

double ExSol(double t, double x0, double t0) {
    //return exp(t) + 1. / (t0 + 1. / (x0 - exp(t0)) - t);
    return t + 2 + 4./(exp(4*(t - t_0))*(4./(x_0 - t_0 - 2) + 1)-1);
}

int SolverRik(double *x, double *t, double a(double), double b(double), double c(double), double x_0, double t_0) {
    
    x[0] = x_0;
    t[0] = t_0;
    int m = 0;
    while (t[m] < T) {

        w_1 = f(x[m], t[m], a, b, c);
        w_2 = f(x[m] + tau * 1. / 2 * w_1, t[m] + tau * 1. / 2, a, b, c);
        w_3 = f(x[m] + tau * 1. / 2 * w_2, t[m] + tau * 1. / 2, a, b, c);
        w_4 = f(x[m] + tau * w_3, t[m] + tau, a, b, c);
        x_main = x[m] + tau * (1. / 6 * w_1 + 1. / 3 * w_2 + 1. / 3 * w_3 + 1. / 6 * w_4);

        x_curr = x[m];

        for (int k = 0; k <= 1; k++) {
            w_1 = f(x_curr, t[m] + k * tau / 2, a, b, c);
            w_2 = f(x_curr + tau / 2 * 1. / 2 * w_1, t[m] + k * tau / 2 + tau / 2 * 1. / 2, a, b, c);
            w_3 = f(x_curr + tau / 2 * 1. / 2 * w_2, t[m] + k * tau / 2 + tau / 2 * 1. / 2, a, b, c);
            w_4 = f(x_curr + tau / 2 * w_3, t[m] + k * tau / 2 + tau / 2, a, b, c);
            x_curr = x_curr + tau / 2 * (1. / 6 * w_1 + 1. / 3 * w_2 + 1. / 3 * w_3 + 1. / 6 * w_4);
        }

        x_new = x_curr;

        tau_old = tau;
        tau_opt = tau * pow(eps * (pow(2, p) - 1) / (T - t_0) / abs(x_main - x_new), (1. / p));


        tau = tau_opt;
        w_1 = f(x[m], t[m], a, b, c);
        w_2 = f(x[m] + tau * 1. / 2 * w_1, t[m] + tau * 1. / 2, a, b, c);
        w_3 = f(x[m] + tau * 1. / 2 * w_2, t[m] + tau * 1. / 2, a, b, c);
        w_4 = f(x[m] + tau * w_3, t[m] + tau, a, b, c);
        x[m + 1] = x[m] + tau * (1. / 6 * w_1 + 1. / 3 * w_2 + 1. / 3 * w_3 + 1. / 6 * w_4);

        t[m + 1] = t[m] + tau;
        m = m + 1;
    }
    return m;
}

int main() {

    int m = 0;
    m = SolverRik(x, t, a, b, c, x_0, t_0);
    

    std::cout << m << std::endl;
    //Запись в файл и вывод на экран

    ofstream fout;
    fout.open("X_Sol.txt");
    
    for (int k = 0; k <= m; k++) {
        std::cout << t[k] << "  " << x[k] << "  " << abs(ExSol(t[k], x_0, t_0) - x[k]) << std::endl;
        fout << t[k] << "\t" << x[k] << "\t" << abs(ExSol(t[k], x_0, t_0) - x[k]) << "\n";
       }
    std::cout << std::endl;
    fout.close();

    return 0;
}