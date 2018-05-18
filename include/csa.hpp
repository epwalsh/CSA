#include <omp.h>
#include <vector>


namespace CSA {

template<typename Scalar>
class Solver
{
private:
public:
    inline int minimize(const int n, Scalar* x) {
        for (int i = 0; i < n; ++i)
            x[i] = 0;
        return 0;
    }
};

}  // namespace CSA
