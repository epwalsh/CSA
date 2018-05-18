#include <omp.h>
#include <vector>


namespace CSA {

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
    int n_process;
    // The maximum number of iterations/steps.
    int max_iterations = 10000;
    // Specifies how many steps in between updates to the generation and
    // acceptance temperatures.
    int update_interval = 5;
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
        Scalar_fx res = fx(instance, x);
        printf("fx: %f\n", res);

        for (int i = 0; i < n; ++i)
            x[i] = 0;

        return 0;
    }
};

}  // namespace CSA
