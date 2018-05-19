#include <cmath>
#include <omp.h>
#include <vector>


namespace CSA {

template<typename Scalar_x, typename Scalar_fx>
class State
{
public:
    /*
     * =========================================================================
     * Public member variables.
     * =========================================================================
     */

    std::vector<Scalar_x> x;
    std::vector<Scalar_x> best_x;
    Scalar_fx cost;
    Scalar_fx best_cost;

    /*
     * =========================================================================
     * Constructors.
     * =========================================================================
     */

    State(int n, const Scalar_x* x0, Scalar_fx fx0)
    {
        this->x = std::vector<Scalar_x>(x0, x0 + n);
        this->best_x = std::vector<Scalar_x>(x0, x0 + n);
        this->cost = fx0;
        this->best_cost = fx0;
    }

    /*
     * =========================================================================
     * Public member functions.
     * =========================================================================
     */

    void step(std::vector<Scalar_x> &y)
    {
        this->x.swap(y);
    }
};


template<typename Scalar_x, typename Scalar_fx>
class SharedStates
{
public:
    /*
     * =========================================================================
     * Public member variables.
     * =========================================================================
     */

    // The number of shared states.
    const int m;
    // The dimension of the shared states.
    const int n;
    // Vector of states.
    std::vector<State<Scalar_x, Scalar_fx>> states;

    /*
     * =========================================================================
     * Public member functions and operators.
     * =========================================================================
     */

    State<Scalar_x, Scalar_fx>& operator[](int i) {
        return this->states[i];
    }

    const State<Scalar_x, Scalar_fx>& operator[](int i) const {
        return this->states[i];
    }

    /*
     * =========================================================================
     * Constructors.
     * =========================================================================
     */

    SharedStates(int m, int n, const Scalar_x* x0, Scalar_fx fx0)
        : m(m), n(n), states(m, State<Scalar_x, Scalar_fx>(n, x0, fx0))
    {
    }
};


template<typename Scalar_x, typename Scalar_fx>
class Solver
{
public:
    /*
     * =========================================================================
     * Public member variables.
     * =========================================================================
     */

    // The number of coupled annealing processes.
    int m = 4;

    // The maximum number of iterations/steps.
    int max_iterations = 10000;

    // The initial value of the generation temperature.
    float tgen_initial = 0.01;

    // Determines the factor that tgen is multiplied by during each update.
    float tgen_schedule = 0.99999;

    // The initial value of the acceptance temperature.
    float tacc_initial = 0.9;

    // Determines the factor by which `tacc` is increased or decreased during each
    // update.
    float tacc_schedule = 0.01;

    // The desired variance of the acceptance probabilities.
    float desired_variance = 0.99;

    /*
     * =========================================================================
     * Constructors.
     * =========================================================================
     */

    Solver() {  };

    /*
     * =========================================================================
     * Public member functions.
     * =========================================================================
     */

    inline int minimize(int n,
                        Scalar_x* x,
                        Scalar_fx (*fx)(void*, Scalar_x*),
                        void (*step)(void*, Scalar_x* y, const Scalar_x*),
                        void (*progress)(void*),
                        void* instance)
    {
        Scalar_fx fx0 = fx(instance, x);
        printf("Initial cost: %f\n", fx0);

        // Initialize shared values.
        SharedStates<Scalar_x, Scalar_fx> shared_states(this->m, n, x, fx0);
        float tacc = this->tacc_initial;
        float tgen = this->tgen_initial;
        float gamma, tmp, sum_a, prob_var;

        omp_lock_t lock;
        omp_init_lock(&lock);

        #pragma omp parallel shared(n, shared_states, tacc, tgen, gamma) num_threads(this->m) default(none)
        {
            int k, opt_id = omp_get_thread_num();

            Scalar_fx max_cost = shared_states[0].cost;
            Scalar_fx cost;
            std::vector<Scalar_x> y(n, Scalar_x(0));

            #pragma omp for
            for (int iter = 0; iter < this->max_iterations; ++iter) {
                // Generate a new solution.
                step(instance, y.data(), shared_states[opt_id].x.data());
                cost = fx(instance, y.data());

                // Take step and update best solution so far.
                // TODO

                // Update temperatures.
                if (omp_test_lock(&lock)) {

                    // Get maximum cost.
                    max_cost = shared_states[0].cost;
                    for (k = 0; k < this->m; ++k)
                        if (shared_states[k].cost > max_cost)
                            max_cost = shared_states[k].cost;

                    // Update gamma and prob_var;
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
            }
        }  // End of #pragma omp parallel

        // Get best result.
        int best_ind = 0;
        Scalar_fx best_cost = shared_states[0].best_cost;
        for (int k = 0; k < this->m; ++k) {
            if (shared_states[k].best_cost < best_cost) {
                best_cost = shared_states[k].best_cost;
                best_ind = k;
            }
        }
        printf("Best cost: %f\n", best_cost);
        State<Scalar_x,Scalar_fx> best_state = shared_states[best_ind];
        for (int i = 0; i < n; ++i)
            x[i] = best_state.best_x[i];

        // Clean up.
        omp_destroy_lock(&lock);

        return 0;
    }
};

}  // namespace CSA
