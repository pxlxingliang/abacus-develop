#include "module_basis/module_ao/parallel_2d.h"
#include "module_base/macros.h"

#ifdef __MPI
#include <mpi.h>
#endif

namespace hsolver
{


#ifdef __MPI

/**
 * @brief Parallel do the generalized eigenvalue problem
 * 
 * @tparam T double or complex<double> or float or complex<float>
 * @param H the hermitian matrix H in H x=lambda S x.
 * @param S the overlap matrix S in H x=lambda S x.
 * @param lda the leading dimension of H and S
 * @param nband the number of bands to be calculated
 * @param ekb to store the eigenvalues.
 * @param wfc to store the eigenvectors
 * @param comm the communicator
 * @param diago_dav_method the method to solve the generalized eigenvalue problem
 * 
 * @note 1. h and s should be full matrix in rank 0 of the communicator, and the other ranks is not concerned.
 */
template <typename T>
void Diago_HS_para(
                T* h, 
                T* s, 
                 const int lda, 
                 const int nband,
                 typename GetTypeReal<T>::type *const ekb,
                 T *const wfc,
                 const MPI_Comm& comm,
                 const int block_size,
                 const int diago_dav_method);
#endif

} // namespace hsolver                 
              