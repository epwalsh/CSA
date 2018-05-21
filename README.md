# CSA

**CSA** is a C++ header-only library for
[Coupled Simulated Annealing](ftp://ftp.esat.kuleuven.be/sista/sdesouza/papers/CSA2009accepted.pdf).
The code is derived and modified from the original implementation by the authors of the paper.

## Features

- **Global optimization of arbitrary functions:** CSA represents a class of global optimization
algorithms that do not make use of the gradient.
- **Highly parallel:** The source code is implemented with OpenMP.
- **Callback interface:** The optimization method recieves the current cost
and next step via a callback interface. Progress updates can also be given
through a callback interface.

## Requirements

A compiler with OpenMP support.

## Quick example

In this example we will try to minimize the
[Schwefel function](https://www.sfu.ca/~ssurjano/schwef.html).
We start by defining the function itself:

```.cpp
#include <cmath>

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
```

The next thing we need is the random step function:

```.cpp
#define PI 3.14159265358979323846264338327


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
```

Optionally, we can define a progress function, which will print updates
to the terminal every time a new minimum cost is found:

```.cpp
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
```

The minimization process is then ran like this:

```.cpp
#include <iostream>

#include <csa.hpp>


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
```

The last few lines of output will look something like this:

```
...
bestcost=1.338e-04       tgen=1.161e-05          tac=7.006e-44        thread=1
bestcost=1.332e-04       tgen=9.740e-06          tac=7.006e-44        thread=1
bestcost=1.314e-04       tgen=9.240e-06          tac=7.006e-44        thread=1
bestcost=1.285e-04       tgen=5.911e-06          tac=7.006e-44        thread=1
bestcost=1.283e-04       tgen=3.647e-06          tac=7.006e-44        thread=1
bestcost=1.280e-04       tgen=3.523e-06          tac=7.006e-44        thread=1
bestcost=1.279e-04       tgen=3.337e-06          tac=7.006e-44        thread=1
bestcost=1.277e-04       tgen=2.924e-06          tac=7.006e-44        thread=1
bestcost=1.276e-04       tgen=2.187e-06          tac=7.006e-44        thread=1
bestcost=1.275e-04       tgen=2.129e-06          tac=7.006e-44        thread=1
bestcost=1.275e-04       tgen=1.880e-06          tac=7.006e-44        thread=1
bestcost=1.275e-04       tgen=1.507e-06          tac=7.006e-44        thread=1
bestcost=1.275e-04       tgen=1.470e-06          tac=7.006e-44        thread=1
Best cost: 0.000128
x =
420.969 420.969 420.968 420.969 420.968 420.969 420.968 420.969 420.969 420.969
```

This is pretty good, since the true minimum is `0 = f(420.9687, ..., 420.9687)`.
Check out [example.cpp](https://github.com/epwalsh/CSA/blob/master/example.cpp)
for the full source code to this example.

## Documentation

You can find the full API reference on the project website:
[epwalsh.github.io/CSA/doc/](https://epwalsh.github.io/CSA/api/).
