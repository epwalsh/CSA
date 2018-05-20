#include <cmath>
#include <iostream>

#include <csa.hpp>


#ifndef PI
    #define PI 3.14159265358979323846264338327
#endif

#define DIM 10


// Generate a pseudo random number according to a UNIF(min, max) distribution.
double uniform(double min, double max)
{
    return min + (max - min) * drand48();
}


// Schwefel function
double f(void* instance, double* x)
{
    double sum = 0.;
    for (int i = 0; i < DIM; ++i)
        sum += 500 * x[i] * std::sin(std::sqrt(std::fabs(500 * x[i])));
    return 418.9829 * DIM - sum;
}


void step(void* instance, double* y, const double* x, float tgen)
{
    int i;
    double tmp;
    for (i = 0; i < DIM; ++i) {
        tmp = std::fmod(x[i] + tgen * std::tan(PI * (drand48() - 0.5)), 1.);
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

    double* x = new double[DIM];
    for (int i = 0; i < DIM; ++i)
        x[i] = drand48();

    CSA::Solver<double, double> solver;
    solver.m = 2;
    solver.minimize(DIM, x, f, step, progress, nullptr);

    for (int i = 0; i < DIM; ++i)
        std::cout << 500 * x[i] << " ";
    std::cout << std::endl;

    delete[] x;

    return EXIT_SUCCESS;
}
