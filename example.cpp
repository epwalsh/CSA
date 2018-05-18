#include <iostream>

#include <csa.hpp>


const int n = 10;


double f(void* instance, double* x)
{
    return 0.;
}


void step(void* instance, double* x)
{
    for (int i = 0; i < n; ++i)
        x[i] = 0.;
}


void progress(void* instance)
{
    return ;
}


int main()
{
    CSA::Solver<double, double> solver;
    double* x = new double[n];

    int res = solver.minimize(n, x, f, step, progress, nullptr);

    for (int i = 0; i < n; ++i)
        std::cout << x[i] << " ";
    std::cout << std::endl;
    std::cout << res << std::endl;

    delete[] x;

    return EXIT_SUCCESS;
}
