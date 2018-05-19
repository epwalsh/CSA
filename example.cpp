#include <cmath>
#include <iostream>
#include <random>

#include <csa.hpp>


// The dimension.
const int n = 10;


// Generate a pseudo random number according to a UNIF(min, max) distribution.
double uniform(double min, double max)
{
    return min + (max - min) * rand() / RAND_MAX;
}


// Schwefel function
double f(void* instance, double* x)
{
    double sum = 0.;
    for (int i = 0; i < n; ++i)
        sum += x[i] * std::sin(std::sqrt(std::fabs(x[i])));
    return 418.9829 * n - sum;
}


void step(void* instance, double* y, const double* x)
{
    int i;
    double tmp;
    for (i = 0; i < n; ++i) {
        tmp = x[i] + uniform(-10, 10);
        if (tmp < -500.)
            tmp = -500.;
        else if (tmp > 500.)
            tmp = 500.;
        y[i] = tmp;
    }
}


void progress(void* instance)
{
    return ;
}


int main()
{
    srand(0);

    double* x = new double[n];
    for (int i = 0; i < n; ++i)
        x[i] = uniform(-500., 500.);

    CSA::Solver<double, double> solver;
    solver.minimize(n, x, f, step, progress, nullptr);

    for (int i = 0; i < n; ++i)
        std::cout << x[i] << " ";
    std::cout << std::endl;

    delete[] x;

    return EXIT_SUCCESS;
}
