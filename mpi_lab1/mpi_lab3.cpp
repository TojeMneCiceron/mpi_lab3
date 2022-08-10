#include <iostream>
#include "mpi.h"
#include <iomanip>
#include <fstream>
#include <vector>
using namespace std;

enum DIRECTIONS { LEFT, RIGHT };
const char* neighbours_names[] = { "left", "right" };

const double T = 1;

const double a = 0, b = 1, c = 0, d = T;

double f(double x)
{
    if (x < 0.1)
        return 120 * x + 8;
    return -170.0 / 9 * x + 197.0 / 9;
}

void print_res_u(vector<vector<double>> u, ostream& fout, double h)
{
    fout << 0;

    for (int i = 0; i < u.size(); i++)
        fout << ',' << fixed << setprecision(2) << i * h;
    fout << endl;

    for (int i = 0; i < u.size(); i++)
    {
        fout << i * h;
        for (int j = 0; j < u[i].size(); j++)
            fout << ',' << fixed << setprecision(2) << u[i][j];
        fout << endl;
    }
}

void print_u(vector<vector<double>> u, int rank)
{
    cout << endl << endl << rank << endl;

    for (int i = 0; i < u.size(); i++)
    {
        for (int j = 0; j < u[i].size(); j++)
            cout << fixed << setprecision(2) << u[i][j] << '\t';
        cout << endl;
    }
}

vector<double> pr(int n, double tau, double h, vector<vector<double>> u, vector<double> y1, int j,
    MPI_Comm comm_1d, int neighbours_ranks[], int size, int rank)
{
    vector<double> alpha(size), beta(size);

    double gamma = tau / h / h;

    //double prev_y1, prev_beta;

    if (neighbours_ranks[LEFT] == MPI_PROC_NULL)
    {
        alpha[0] = 0;
        beta[0] = u[0][j];
    }
    else
    {
        //ждем начальное значение альфы и беты
        MPI_Recv(&alpha[0], 1, MPI_DOUBLE, neighbours_ranks[LEFT], 0, comm_1d, MPI_STATUSES_IGNORE);
        MPI_Recv(&beta[0], 1, MPI_DOUBLE, neighbours_ranks[LEFT], 0, comm_1d, MPI_STATUSES_IGNORE);

        //cout << "\n alpha " << alpha[0];
        //cout << "\n beta " << beta[0];
    }
    for (int i = 0; i < size - 1; i++)
    {
        alpha[i + 1] = gamma / (1 + 2 * gamma - gamma * alpha[i]);
        beta[i + 1] = (gamma * beta[i] + y1[i]) / (1 + 2 * gamma - gamma * alpha[i]);
    }

    if (neighbours_ranks[RIGHT] != MPI_PROC_NULL)
    {
        //отправляем, когда посчитали
        MPI_Send(&alpha[size - 1], 1, MPI_DOUBLE, neighbours_ranks[RIGHT], 0, comm_1d);
        //cout << "\n alpha " << alpha[size - 1];
        MPI_Send(&beta[size - 1], 1, MPI_DOUBLE, neighbours_ranks[RIGHT], 0, comm_1d);
        //cout << "\n beta " << beta[size - 1];
    }

    /*for (int i = 0; i < size - 1; i++)
    {
        alpha[i + 1] = gamma / (1 + 2 * gamma - gamma * alpha[i]);
        beta[i + 1] = (gamma * beta[i] + y1[i]) / (1 + 2 * gamma - gamma * alpha[i]);
    }*/

    vector<double> y(size);

    if (neighbours_ranks[RIGHT] == MPI_PROC_NULL)
    {
        y[size - 1] = u[size - 1][j];
    }
    else
    {
        //ждем начальное значение
        MPI_Recv(&y[size - 1], 1, MPI_DOUBLE, neighbours_ranks[RIGHT], 0, comm_1d, MPI_STATUSES_IGNORE);
        //cout << endl << rank << " " << y[size - 1];
    }
    //if (neighbours_ranks[LEFT] != MPI_PROC_NULL)
    //    MPI_Send(&y[0], 1, MPI_DOUBLE, neighbours_ranks[LEFT], 0, comm_1d);

    for (int i = size - 2; i >= 0; i--)
        y[i] = alpha[i + 1] * y[i + 1] + beta[i + 1];

    //отправляем, когда посчитали
    if (neighbours_ranks[LEFT] != MPI_PROC_NULL)
        MPI_Send(&y[0], 1, MPI_DOUBLE, neighbours_ranks[LEFT], 0, comm_1d);

    return y;
}

vector<vector<double>> P(int m1, int n, double tau, double h, 
    vector<vector<double>> u, MPI_Comm comm_1d, int neighbours_ranks[], int rank, int size, int m_)
{
    double t = 0;

    vector<double> y1(size);

    int offset = size * rank;

    for (int i = offset; i < size + offset; i++)
        y1[i - offset] = f(i * h);

    for (int k = 0; k <= m1; k++)
    {
        //cout << k << endl;
        for (int i = 0; i < size; i++)
        {
            if (u[i][k] == 0)
                u[i][k] = y1[i];
        }
        for (int j = 0; j < m_; j++)
        {
            y1 = pr(n, tau, h, u, y1, m1, comm_1d, neighbours_ranks, size, rank);

            t += tau;
        }
    }

    return u;
}

vector<vector<double>> init_u_part(int rank, int neighbours_ranks[], int size, double h, int n)
{
    vector<vector<double>> u_part(size);

    for (int i = 0; i < size; i++)
        u_part[i].resize(n + 1);

    if (neighbours_ranks[LEFT] == MPI_PROC_NULL)
    {
        for (int i = 0; i <= n; i++)
            u_part[0][i] = f(0);
    }
    
    if (neighbours_ranks[RIGHT] == MPI_PROC_NULL)
    {
        for (int i = 0; i <= n; i++)
            u_part[size - 1][i] = f(1);
    }

    /*if (neighbours_ranks[LEFT] == MPI_PROC_NULL)
    {
        if (neighbours_ranks[UP] == MPI_PROC_NULL)
        {
            for (int i = 0; i < size; i++)
                u_part[i][0] = f(i*h);
        }
        else
        {
            for (int i = size; i < size * 2; i++)
                u_part[i - size][0] = f(i*h);
        }
    }*/

    return u_part;
}

double* v_to_d(vector<vector<double>> u)
{
    double* res = new double[u.size()* u[0].size()];
    for (int i = 0; i < u.size(); i++)
        for (int j = 0; j < u[i].size(); j++)
            res[i* u[i].size() + j] = u[i][j];
    return res;
}

void print_d(double* d, int count, int n, int m)
{
    //cout << count << " " << m;
    cout << endl;
    for (int i = 0; i < count; i++)
    {
        cout << fixed << setprecision(2) << d[i] << "\t";
        if ((i + 1) % m == 0)
            cout << endl;
    }
}

void print_to_f(double* d, int count, int n, int m, ostream& fout, double h)
{
    //cout << count << " " << m;
    fout << 0;

    for (int i = 0; i < n; i++)
        fout << ',' << fixed << setprecision(2) << i * h;
    fout << endl;

    int j = 0;
    fout << j * h;
    for (int i = 0; i < count; i++)
    {
        fout << "," << fixed << setprecision(2) << d[i];
        if ((i + 1) % m == 0)
        {
            fout << endl;
            j++;
            if (i < count - 1)
                fout << j * h;
        }
    }
}

int main(int argc, char** argv)
{
    MPI_Status status;
	MPI_Init(&argc, &argv);
    double time = MPI_Wtime();

    int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int dims[1] = { 0 };
    MPI_Dims_create(size, 1, dims);

    int periods[1] { false };

    MPI_Comm comm_1d;

    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, true, &comm_1d);

    int neighbours_ranks[2];
    MPI_Cart_shift(comm_1d, 0, 1, &neighbours_ranks[LEFT], &neighbours_ranks[RIGHT]);

    MPI_Comm_rank(comm_1d, &rank);

    int n = 63;
    double h = (b - a) / n;

    int m = n;
    int m_ = 2;
    int m1 = m / m_;
    double tau = (d - c) / m;

    int part_size = (n + 1) / size;

    vector<vector<double>> u_part = init_u_part(rank, neighbours_ranks, part_size, h, n);

    vector<vector<double>> new_u_part = P(m1, n, tau, h, u_part, comm_1d, neighbours_ranks, rank, part_size, m_);

    double* part = v_to_d(new_u_part);

    double* u = new double[1];

    if (rank == 0) {

        u = new double[(n+1)*(m1+1)];
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gather(part, part_size*(m1+1), MPI_DOUBLE, u, part_size * (m1 + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        ofstream fout("output.txt");
        print_d(u, (n + 1) * (m1 + 1), n + 1, m1 + 1);
        //print_to_f(u, (n + 1) * (m1 + 1), n + 1, m1 + 1, fout, h);
        fout.close();
    }

    time = MPI_Wtime() - time;

    cout << endl << rank << " time = " << fixed << setprecision(4) << time;

    double max_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0)
        cout << "\ntime = " << fixed << setprecision(4) << max_time;

	MPI_Finalize();
	return 0;
}