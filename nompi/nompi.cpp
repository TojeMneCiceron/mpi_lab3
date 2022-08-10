#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
using namespace std;

const double T = 1;

const double a = 0, b = 1, c = 0, d = T;

double f(double x)
{
    if (x < 0.1)
        return 120 * x + 8;
    return -170.0 / 9 * x + 197.0 / 9;
}

void print_u(vector<vector<double>> u, ostream& fout, double h)
{
    fout << 0;

    for (int i = 0; i < u.size(); i++)
        fout << ',' << fixed << setprecision(5) << i * h;
    fout << endl;

    for (int i = 0; i < u[i].size(); i++)
    {
        fout << i * h;
        for (int j = 0; j < u.size(); j++)
            fout << ',' << fixed << setprecision(5) << u[j][i];
        fout << endl;
    }
}

vector<double> pr(int n, double tau, double h, vector<vector<double>> u, vector<double> y1, int j)
{
    vector<double> alpha(n + 1), beta(n + 1);
    alpha[0] = 0;
    beta[0] = u[0][j];

    double gamma = tau / h / h;

    for (int i = 0; i < n; i++)
    {
        alpha[i + 1] = gamma / (1 + 2 * gamma - gamma * alpha[i]);
        beta[i + 1] = (gamma * beta[i] + y1[i]) / (1 + 2 * gamma - gamma * alpha[i]);
    }

    vector<double> y(n + 1);
    y[n] = u[n][j];

    for (int i = n - 1; i >= 0; i--)
        y[i] = alpha[i + 1] * y[i + 1] + beta[i + 1];

    return y;
}

vector<vector<double>> P(int m1, int n, double tau, double h, vector<vector<double>> u, int m_)
{
    double t = 0;

    vector<double> y1(n + 1);

    for (int i = 0; i <= n; i++)
        y1[i] = f(i * h);

    for (int k = 0; k <= m1; k++)
    {
        //cout << k << endl;
        for (int i = 1; i <= n; i++)
        {
            u[i][k] = y1[i];
        }
        for (int j = 0; j < m_; j++)
        {
            y1 = pr(n, tau, h, u, y1, m1);

            t += tau;
        }
    }

    return u;
}

int main()
{
    double time = clock();

    ofstream fout("output.txt");

    int n = 63;
    double h = (b - a) / n;

    int m = n;
    int m_ = 2;
    int m1 = m / m_;
    double tau = (d - c) / m;

    //матрица
    vector<vector<double>> u(n + 1);
    for (int i = 0; i <= n; i++)
        u[i].resize(m1 + 1);

    double u_0_t = f(0), u_1_t = f(1);

    for (int i = 0; i <= n; i++)
        u[i][0] = f(i * h);

    for (int i = 0; i <= m1; i++)
    {
        u[0][i] = u_0_t;
        u[n][i] = u_1_t;
    }

    vector<vector<double>> newU = P(m1, n, tau, h, u, m_);

    print_u(newU, fout, h);

    fout.close();

    cout << "time = " << (clock() - time) / CLOCKS_PER_SEC;
}