#include <cmath>
#include <iostream>

#include <csa.hpp>


#define PI 3.14159265358979323846264338327
#define DIM 10


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


void progress(
    void* instance, double cost, float tgen, float tacc, int opt_id, int iter)
{
    printf(
        "bestcost=%1.3e \t tgen=%1.3e \t tac=%1.3e "
        "thread=%d \t iter=%d\n",
        cost,
        tgen,
        tacc,
        opt_id,
        iter);
    return ;
}


int main()
{
    srand(0);

    // Create initial solution.
    double* x = new double[DIM];
    for (int i = 0; i < DIM; ++i)
        x[i] = drand48();
    double cost = f(nullptr, x);
    printf("Initial cost: %f\n", cost);

    // Initialize CSA solver.
    CSA::Solver<double, double> solver;
    solver.m = 2; // number of threads = number of solvers = 2

    // Start the annealing process.
    solver.minimize(DIM, x, f, step, progress, nullptr);

    cost = f(nullptr, x);
    printf("Best cost: %f\nx =\n", cost);
    for (int i = 0; i < DIM; ++i)
        std::cout << 500 * x[i] << " ";
    std::cout << std::endl;

    // Clean up.
    delete[] x;

    return EXIT_SUCCESS;
}
