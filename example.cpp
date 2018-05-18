#include <iostream>

#include <csa.hpp>


int main()
{
    CSA::Solver<double> solver;
    int n = 10;
    double* x = new double[10];

    int res = solver.minimize(n, x);

    for (int i = 0; i < n; ++i)
        std::cout << x[i] << " ";
    std::cout << std::endl;
    std::cout << res << std::endl;

    delete[] x;

    return EXIT_SUCCESS;
}
