/*
 *     Example of using CSA minimizing the Schwefel function, `f`. The domain of 
 *     `f` is [-500, 500]^d, where `d` is the dimension. The true minimum of `f` 
 *     is 0.0 at x = (420.9687, ..., 420.9687).
 *
 *
 * Copyright (c) 2009 Samuel Xavier-de-Souza, Johan A.K. Suykens,
 *                    Joos Vandewalle, De ́sire ́ Bolle ́
 * Copyright (c) 2018 Evan Pete Walsh 
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#include <cmath>
#include <iostream>

#include <csa.hpp>


#define PI 3.14159265358979323846264338327
#define DIM 10


// This defines the Schwefel function. The pointer `instance` can be used to
// pass data to `f`, such as the dimension.
// In this case we don't need to pass anything, so it will be given as NULL.
double f(void* instance, double* x)
{
    double sum = 0.;
    for (int i = 0; i < DIM; ++i)
        sum += 500 * x[i] * std::sin(std::sqrt(std::fabs(500 * x[i])));
    return 418.9829 * DIM - sum;
}

// This function will take a random step from `x` to `y`. The value `tgen`,
// the "generation temperature", determines the variance of the distribution of
// the step. `tgen` will decrease according to fixed a schedule throughout the
// annealing process, which corresponds to a decrease in the variance of steps.
void step(void* instance, double* y, const double* x, float tgen)
{
    int i;
    double tmp;
    for (i = 0; i < DIM; ++i) {
        tmp = std::fmod(x[i] + tgen * std::tan(PI * (drand48() - 0.5)), 1.);
        y[i] = tmp;
    }
}

// This will receive progress updates from the CSA process and print updates
// to the terminal.
void progress(
    void* instance, double cost, float tgen, float tacc, int opt_id, int iter)
{
    printf(
        "bestcost=%1.3e \t tgen=%1.3e \t tacc=%1.3e \t thread=%d\n",
        cost,
        tgen,
        tacc,
        opt_id);
    return ;
}


int main()
{
    srand(0);

    // Create initial solution from uniform distribution.
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

    // The array `x` now holds the best solution.
    cost = f(nullptr, x);
    printf("Best cost: %f\nx =\n", cost);
    for (int i = 0; i < DIM; ++i)
        std::cout << 500 * x[i] << " ";
    std::cout << std::endl;

    // Clean up.
    delete[] x;

    return EXIT_SUCCESS;
}
