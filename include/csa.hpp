/*
 *     C++ library for Coupled Simulated Annealing.
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

#ifndef CSA_H
#define CSA_H

#ifndef CSA_ITER_MONO
#define CSA_ITER_MONO 0
#endif

#include <cmath>
#include <omp.h>
#include <vector>


namespace CSA {

template<typename Scalar_x, typename Scalar_fx>
class State
{
public:
    ///
    /// The current solution.
    ///
    std::vector<Scalar_x> x;
    ///
    /// The solution with the best cost so far.
    ///
    std::vector<Scalar_x> best_x;
    ///
    /// The current cost at `x`.
    ///
    Scalar_fx cost;
    ///
    /// The cost associated with `best_x`.
    ///
    Scalar_fx best_cost;

    ///
    /// Contruct a state from an initial solution.
    /// 
    /// \param n   The dimension of `x0`.
    /// \param x0  The initial solution guess.
    /// \param fx0 The value of the cost function associated with `x0`.
    ///
    State(int n, const Scalar_x* x0, Scalar_fx fx0)
    {
        this->x = std::vector<Scalar_x>(x0, x0 + n);
        this->best_x = std::vector<Scalar_x>(x0, x0 + n);
        this->cost = fx0;
        this->best_cost = fx0;
    }

    ///
    /// Move the current solution to `y`. Internally, this just swaps the
    /// pointers of `x` and `y`.
    /// 
    /// \param y      The new solution.
    /// \param y_cost The value of the cost function associated with `y`.
    ///
    void step(std::vector<Scalar_x> &y, Scalar_fx y_cost)
    {
        this->cost = y_cost;
        this->x.swap(y);
    }
};  // class State


template<typename Scalar_x, typename Scalar_fx>
class SharedStates
{
public:
    ///
    /// The number of shared states.
    ///
    const int m;
    ///
    /// The dimension of the shared states.
    ///
    const int n;
    ///
    /// Vector of `States`.
    ///
    std::vector<State<Scalar_x, Scalar_fx>> states;

    State<Scalar_x, Scalar_fx>& operator[](int i) {
        return this->states[i];
    }

    const State<Scalar_x, Scalar_fx>& operator[](int i) const {
        return this->states[i];
    }

    ///
    /// Construct shared states from an initial solution.
    /// 
    /// \param m   The number of threads/shared states.
    /// \param n   The dimension of `x0`.
    /// \param x0  The initial solution guess. Each thread will start from the
    ///            same initial solution.
    /// \param fx0 The value of the cost function associated with `x0`.
    ///
    SharedStates(int m, int n, const Scalar_x* x0, Scalar_fx fx0)
        : m(m), n(n), states(m, State<Scalar_x, Scalar_fx>(n, x0, fx0))
    {
    }
};  // class SharedStates


template<typename Scalar_x, typename Scalar_fx>
class Solver
{
public:
    ///
    /// The number of threads and coupled annealing processes.
    ///
    int m = 4;
    ///
    /// The maximum number of iterations/steps.
    ///
    int max_iterations = 1000000;
    ///
    /// The initial value of the generation temperature.
    ///
    float tgen_initial = 0.01;
    ///
    /// Determines the factor that `tgen` is multiplied by during each update.
    ///
    float tgen_schedule = 0.99999;
    ///
    /// The initial value of the acceptance temperature.
    ///
    float tacc_initial = 0.9;
    ///
    /// Determines the factor by which `tacc` is increased or decreased during
    /// each update.
    ///
    float tacc_schedule = 0.01;
    ///
    /// The desired variance of the acceptance probabilities.
    ///
    float desired_variance = 0.99;

    ///
    /// Default constructor.
    ///
    Solver() {  };

    ///
    /// Run the CSA process to minimize the target function.
    /// 
    /// \param n        The size of the input array `x`.
    /// \param x        The input array, representing an initial guess of the
    ///                 solution. This array will be modified with the
    ///                 best solution found.
    /// \param fx       The function to minimize. The first argument to the
    ///                 function is a pointer (possibly NULL) to the callback
    ///                 object `instance`.
    /// \param step     The step function. Should populate the array `y` with a
    ///                 a random step based on the current position `x` and the
    ///                 generation schedule `tgen`. Like `fx` the first argument
    ///                 to the function should be a pointer to the callback
    ///                 object `instance`.
    /// \param progress An optional function that receives updates when a new
    ///                 best solution is found. Like `fx` and `step`, the first
    ///                 argument is a pointer to the callback object `instance`.
    ///                 The next arguments, in order, are the current cost, the
    ///                 current generation temperature, the current acceptance
    ///                 temperature, the thread ID, and the iteration number.
    /// \param instance An optional pointer to a callback object.
    ///
    /// \warning By default, the iterations of the main loop are not processed
    /// monotonically by the threads, so the `iter` parameter that the
    /// `progress` function receives will appear to be random. This behavior
    /// can be changed by setting the macro CSA_ITER_MONO to a non-zero value.
    /// However, this may result in a significant slow-down.
    ///
    inline int minimize(
        int n,
        Scalar_x* x,
        Scalar_fx (*fx)(void*, Scalar_x*),
        void (*step)(void*, Scalar_x* y, const Scalar_x*, float tgen),
        void (*progress)(void*,
                         Scalar_fx cost,
                         float tgen,
                         float tacc,
                         int opt_id,
                         int iter),
        void* instance)
    {
        Scalar_fx fx0 = fx(instance, x);

        // Initialize shared values.
        SharedStates<Scalar_x, Scalar_fx> shared_states(this->m, n, x, fx0);
        float tacc = this->tacc_initial;
        float tgen = this->tgen_initial;
        float tmp, sum_a, prob_var, gamma = m;

        omp_lock_t lock;
        omp_init_lock(&lock);

        #pragma omp parallel shared(n, shared_states, tacc, tgen, gamma) num_threads(this->m) default(none)
        {
            int k, opt_id = omp_get_thread_num();

            Scalar_fx max_cost = shared_states[0].cost;
            Scalar_fx cost;
            std::vector<Scalar_x> y(n, Scalar_x(0));
            float unif, prob;

#if CSA_ITER_MONO == 0
            #pragma omp for
#else
            #pragma omp for schedule(monotonic:static)
#endif
            for (int iter = 0; iter < this->max_iterations; ++iter) {
                // Generate a new solution.
                step(instance, y.data(), shared_states[opt_id].x.data(), tgen);
                cost = fx(instance, y.data());

                // Decide if we should take a step and update best solution so far.
                if (cost < shared_states[opt_id].cost) {
                    omp_set_lock(&lock);
                    // Update best solution for this thread.
                    if (cost < shared_states[opt_id].best_cost) {
                        shared_states[opt_id].best_cost = cost;
                        shared_states[opt_id].best_x = y;
                        if (progress != nullptr)
                            progress(instance, cost, tgen, tacc, opt_id, iter);
                    }

                    // Take the step.
                    shared_states[opt_id].step(y, cost);
                    omp_unset_lock(&lock);
                } else {
                    // We accept the "worse" solution with probability `prob`:
                    unif = drand48();
                    prob = std::exp((shared_states[opt_id].cost - max_cost) / tacc) / gamma;
                    if (prob > unif) {
                        omp_set_lock(&lock);
                        shared_states[opt_id].step(y, cost);
                        omp_unset_lock(&lock);
                    }
                }

                // Update temperatures.
                if (omp_test_lock(&lock)) {

                    // Get maximum cost.
                    max_cost = shared_states[0].cost;
                    for (k = 0; k < this->m; ++k)
                        if (shared_states[k].cost > max_cost)
                            max_cost = shared_states[k].cost;

                    // Update `gamma` and `prob_var`;
                    gamma = sum_a = 0.;
                    for (k = 0; k < this->m; ++k) {
                        tmp = (shared_states[k].cost - max_cost) / tacc;
                        gamma += std::exp(tmp);
                        sum_a += std::exp(2.0 * tmp);
                    }
                    prob_var = (this->m * (sum_a / (gamma * gamma)) - 1.) /
                               (this->m * this->m);

                    // Update the acceptance temperature.
                    if (prob_var > this->desired_variance)
                        tacc += this->tacc_schedule * tacc;
                    else
                        tacc -= this->tacc_schedule * tacc;

                    // Update the generation temperature.
                    tgen = this->tgen_schedule * tgen;

                    // Unset the lock.
                    omp_unset_lock(&lock);
                }
            }  // End of "#pragma omp for"
        }  // End of "#pragma omp parallel"

        // Get best result.
        int best_ind = 0;
        Scalar_fx best_cost = shared_states[0].best_cost;
        for (int k = 0; k < this->m; ++k) {
            if (shared_states[k].best_cost < best_cost) {
                best_cost = shared_states[k].best_cost;
                best_ind = k;
            }
        }
        State<Scalar_x,Scalar_fx> best_state = shared_states[best_ind];
        for (int i = 0; i < n; ++i)
            x[i] = best_state.best_x[i];

        // Clean up.
        omp_destroy_lock(&lock);

        return 0;
    }
};  // class Solver

}  // namespace CSA

#endif
