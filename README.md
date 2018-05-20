# CSA

**CSA** is a C++ header-only library for
[Coupled Simulated Annealing](ftp://ftp.esat.kuleuven.be/sista/sdesouza/papers/CSA2009accepted.pdf).
The code is derived and modified from the original implementation by the authors of the paper.

## Features

- **Global optimization of arbitrary functions:** CSA represents a class of global optimization
algorithms that do not make use of the derivative.
- **Highly parallel:** The source code is implemented with OpenMP.
- **Callback interface:** The optimization method recieves the current cost
and next step via a callback interface. Progress updates can also be given
through a callback interface.

## Requirements

A compiler with OpenMP support.

## Quick example

See `example.cpp` for an example of finding the global minimum
of the [Schwefel function](https://www.sfu.ca/~ssurjano/schwef.html).

## Documentation

You can find the full API reference on the project website: [epwalsh.github.io/CSA/](https://epwalsh.github.io/CSA/).
