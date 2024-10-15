#ifdef _OPENMP
#include <omp.h>
#endif

[[cpp11::register]]
int openmp_get_thread_id() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return -1;
#endif
}
