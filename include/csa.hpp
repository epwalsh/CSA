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

    State(int n, const Scalar_x* x0, Scalar_fx fx0) {
        this->x = std::vector<Scalar_x>(x0, x0 + n);
        this->best_x = std::vector<Scalar_x>(x0, x0 + n);
        this->cost = fx0;
        this->best_cost = fx0;
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

    const State<Scalar_x, Scalar_fx>& operator[](int i) const {
        return this->states[i];
    }

    /*
     * =========================================================================
     * Constructors.
     * =========================================================================
     */

    SharedStates(int m, int n, const Scalar_x* x0, Scalar_fx fx0) : m(m), n(n)
    {
        this->states.resize(m, State<Scalar_x, Scalar_fx>(n, x0, fx0));
    }
};


template<typename Scalar_x, typename Scalar_fx>
class Solver
{
private:
public:
    /*
     * =========================================================================
     * Public member variables.
     * =========================================================================
     */

    // The number of coupled annealing processes.
    int m;

    // The maximum number of iterations/steps.
    int max_iterations = 10000;

    // The initial value of the generation temperature.
    Scalar_fx tgen_initial = 0.01;

    // Determines the factor that tgen is multiplied by during each update.
    Scalar_fx tgen_schedule = 0.99999;

    // The initial value of the acceptance temperature.
    Scalar_fx tacc_initial = 0.9;

    // Determines the factor that `tacc` is multiplied by during each update.
    Scalar_fx tagg_schedule = 0.95;

    // The desired variance of the acceptance probabilities.
    Scalar_fx desired_variance = 0.99;

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

    inline int minimize(const int n,
                        Scalar_x* x,
                        Scalar_fx (*fx)(void*, Scalar_x*),
                        void (*step)(void*, Scalar_x*),
                        void (*progress)(void*),
                        void* instance)
    {
        Scalar_fx fx0 = fx(instance, x);
        printf("Initial cost: %f\n", fx0);

        // Initialize shared values.
        SharedStates<Scalar_x, Scalar_fx> shared_states(this->m, n, x, fx0);

        omp_lock_t lock;
        omp_init_lock(&lock);

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
