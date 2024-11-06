#ifndef HSOLVER_DIAGO_PXXXGVX_H
#define HSOLVER_DIAGO_PXXXGVX_H
#include <complex>

namespace hsolver
{

#ifdef __MPI
void pxxxgvx_diag(const int* const desc,
                  const int ncol,
                  const int nrow,
                  const int ndim_global,
                  const int nbands,
                  const double* const h_mat,
                  const double* const s_mat,
                  double* const ekb,
                  double* const wfc_2d);

void pxxxgvx_diag(const int *const desc,
                      const int ncol,
                      const int nrow,
                      const int ndim_global,
                      const int nbands,
                      const std::complex<double> *const h_mat,
                      const std::complex<double> *const s_mat,
                      double *const ekb,
                      std::complex<double> *const wfc_2d);

void pxxxgvx_diag(const int *const desc,
                      const int ncol,
                      const int nrow,
                      const int ndim_global,
                      const int nbands,
                      const float *const h_mat,
                      const float *const s_mat,
                      float *const ekb,
                      float *const wfc_2d);

void pxxxgvx_diag(const int *const desc,
                        const int ncol,
                        const int nrow,
                        const int ndim_global,
                        const int nbands,
                        const std::complex<float> *const h_mat,
                        const std::complex<float> *const s_mat,
                        float *const ekb,
                        std::complex<float> *const wfc_2d);
                                              
#endif 

} // namespace hsolver

#endif // HSOLVER_DIAGO_PXXXGVX_H