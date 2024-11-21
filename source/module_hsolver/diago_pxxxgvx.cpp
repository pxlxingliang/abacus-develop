#include <complex>
#include <vector>
#include <cstring>
#include <iostream>

#include "diago_pxxxgvx.h"
#include "module_base/scalapack_connector.h"
#include "module_base/blacs_connector.h"

namespace hsolver
{

#ifdef __MPI

void pxxxgvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
		const int* n, double* A, const int* ia, const int* ja, const int*desca, double* B, const int* ib, const int* jb, const int*descb,
		const double* vl, const double* vu, const int* il, const int* iu,
		const double* abstol, int* m, int* nz, double* w, const double*orfac, double* Z, const int* iz, const int* jz, const int*descz,
		double* work, int* lwork, double* rwork, int* lrwork, int*iwork, int*liwork, int* ifail, int*iclustr, double*gap, int* info)
{
    // for a uniform interface, rwork and lrwork are input arguments, but not used in pdsygvx_/pssygvx_
    pdsygvx_(itype, jobz, range, uplo,
		n, A, ia, ja, desca, B, ib, jb, descb,
		vl, vu, il,  iu,
		abstol, m, nz, w, orfac, Z, iz, jz, descz,
		work, lwork, iwork, liwork, ifail, iclustr, gap, info);
}

void pxxxgvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
		const int* n, std::complex<double>* A, const int* ia, const int* ja, const int*desca, std::complex<double>* B, const int* ib, const int* jb, const int*descb,
		const double* vl, const double* vu, const int* il, const int* iu,
		const double* abstol, int* m, int* nz, double* w, const double*orfac, std::complex<double>* Z, const int* iz, const int* jz, const int*descz,
		std::complex<double>* work, int* lwork, double* rwork, int* lrwork, int*iwork, int*liwork, int* ifail, int*iclustr, double*gap, int* info)
{
    pzhegvx_(itype, jobz, range, uplo,
		n, A, ia, ja, desca, B, ib, jb, descb,
		vl, vu, il,  iu,
		abstol, m, nz, w, orfac, Z, iz, jz, descz,
		work, lwork,rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info);
}

void pxxxgvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
		const int* n, float* A, const int* ia, const int* ja, const int*desca, float* B, const int* ib, const int* jb, const int*descb,
		const float* vl, const float* vu, const int* il, const int* iu,
		const float* abstol, int* m, int* nz, float* w, const float*orfac, float* Z, const int* iz, const int* jz, const int*descz,
		float* work, int* lwork, float* rwork, int* lrwork, int*iwork, int*liwork, int* ifail, int*iclustr, float*gap, int* info)
{
    pssygvx_(itype, jobz, range, uplo,
		n, A, ia, ja, desca, B, ib, jb, descb,
		vl, vu, il,  iu,
		abstol, m, nz, w, orfac, Z, iz, jz, descz,
		work, lwork, iwork, liwork, ifail, iclustr, gap, info);
}

void pxxxgvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
		const int* n, std::complex<float>* A, const int* ia, const int* ja, const int*desca, std::complex<float>* B, const int* ib, const int* jb, const int*descb,
		const float* vl, const float* vu, const int* il, const int* iu,
		const float* abstol, int* m, int* nz, float* w, const float*orfac, std::complex<float>* Z, const int* iz, const int* jz, const int*descz,
		std::complex<float>* work, int* lwork, float* rwork, int* lrwork, int*iwork, int*liwork, int* ifail, int*iclustr, float*gap, int* info)
{
    pchegvx_(itype, jobz, range, uplo,
		n, A, ia, ja, desca, B, ib, jb, descb,
		vl, vu, il,  iu,
		abstol, m, nz, w, orfac, Z, iz, jz, descz,
		work, lwork,rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info);
}

// post processing for pdsygvx/pzhegvx/pdsygvx/pzhegvx
void pxxxgvx_post_processing(const int info, 
                    const std::vector<int>& ifail,
                    const std::vector<int>& iclustr,
                    const int M,
                    const int NZ,
                    const int nbands, 
                    int& degeneracy_max)
{
    const std::string str_info = "Scalapack diagonalization: \n    info = " + std::to_string(info) + ".\n";

    if (info == 0)
    {
        return;
    }
    else if (info < 0)
    {
        const int info_negative = -info;
        const std::string str_index
            = (info_negative > 100)
                  ? std::to_string(info_negative / 100) + "-th argument "
                        + std::to_string(info_negative % 100) + "-entry is illegal.\n"
                  : std::to_string(info_negative) + "-th argument is illegal.\n";
        throw std::runtime_error(str_info + str_index);
    }
    else if (info % 2)
    {
        std::string str_ifail = "ifail = ";
        for (const int i: ifail)
        {
            str_ifail += std::to_string(i) + " ";
        }
        throw std::runtime_error(str_info + str_ifail);
    }
    else if (info / 2 % 2)
    {
        int degeneracy_need = 0;
        for (int irank = 0; irank < iclustr.size()/2; ++irank)
        {
            degeneracy_need = std::max(degeneracy_need, iclustr[2 * irank + 1] - iclustr[2 * irank]);
        }
        const std::string str_need = "degeneracy_need = " + std::to_string(degeneracy_need) + ".\n";
        const std::string str_saved
            = "degeneracy_saved = " + std::to_string(degeneracy_max) + ".\n";
        if (degeneracy_need <= degeneracy_max)
        {
            throw std::runtime_error(str_info + str_need + str_saved);
        }
        else
        {
            std::cout << str_need << str_saved;
            degeneracy_max = degeneracy_need;
            return;
        }
    }
    else if (info / 4 % 2)
    {
        const std::string str_M = "M = " + std::to_string(M) + ".\n";
        const std::string str_NZ = "NZ = " + std::to_string(NZ) + ".\n";
        const std::string str_NBANDS
            = "Number of eigenvalues solved = " + std::to_string(nbands) + ".\n";
        throw std::runtime_error(str_info + str_M + str_NZ + str_NBANDS);
    }
    else if (info / 16 % 2)
    {
        const std::string str_npos = "Not positive definite = " + std::to_string(ifail[0]) + ".\n";
        throw std::runtime_error(str_info + str_npos);
    }
    else
    {
        throw std::runtime_error(str_info);
    }
}

void get_lwork(int& lwork, std::vector<double>& work)
{
    lwork = work[0];
}

void get_lwork(int& lwork, std::vector<float>& work)
{
    lwork = work[0];
}

void get_lwork(int& lwork, std::vector<std::complex<double>>& work)
{
    lwork = work[0].real();
}

void get_lwork(int& lwork, std::vector<std::complex<float>>& work)
{
    lwork = work[0].real();
}

template <typename T>
void pxxxgvx_diag(const int *const desc,
                      const int ncol,
                      const int nrow,
                      const int nbands,
                      const T *const h_mat,
                      const T *const s_mat,
                      typename GetTypeReal<T>::type *const ekb,
                      T *const wfc_2d)
{
    int nprow, npcol, myprow, mypcol;
    Cblacs_gridinfo(desc[1], &nprow, &npcol, &myprow, &mypcol);
    int dsize = nprow * npcol;

    int degeneracy_max = 12; // only used for complex<float> and complex<double>
    while (true)
    {
        std::vector<T> h_tmp(ncol * nrow);
        std::vector<T> s_tmp(ncol * nrow);
        memcpy(h_tmp.data(), h_mat, sizeof(T) * ncol * nrow);
        memcpy(s_tmp.data(), s_mat, sizeof(T) * ncol * nrow);

        int ndim_global = desc[2];
        const char jobz = 'V', range = 'I', uplo = 'U';
        const int itype = 1, il = 1, iu = nbands, one = 1;
        int M = 0, NZ = 0, lwork = -1, lrwork = -1, liwork = -1, info = 0;
        const typename GetTypeReal<T>::type abstol = 0, orfac = -1;
        const typename GetTypeReal<T>::type vl = 0, vu = 0;
        std::vector<T> work(1, 0);
        std::vector<typename GetTypeReal<T>::type> rwork(3, 0); // only used for complex<float> and complex<double>
        std::vector<int> iwork(1, 0);
        std::vector<int> ifail(ndim_global, 0);
        std::vector<int> iclustr(2 * dsize);
        std::vector<typename GetTypeReal<T>::type> gap(dsize);

        pxxxgvx_(&itype,
                 &jobz,
                 &range,
                 &uplo,
                 &ndim_global,
                 h_tmp.data(),
                 &one,
                 &one,
                 desc,
                 s_tmp.data(),
                 &one,
                 &one,
                 desc,
                 &vl,
                 &vu,
                 &il,
                 &iu,
                 &abstol,
                 &M,
                 &NZ,
                 ekb,
                 &orfac,
                 wfc_2d,
                 &one,
                 &one,
                 desc,
                 work.data(),
                 &lwork, // is not used for real data type
                 rwork.data(), // is not used for real data type
                 &lrwork,
                 iwork.data(),
                 &liwork,
                 ifail.data(),
                 iclustr.data(),
                 gap.data(),
                 &info);

        if (info)
        {
            throw std::runtime_error("Scalapack diagonalization: \n    info = " + std::to_string(info) + ".\n");
        }

        if (std::is_same<T, std::complex<float>>::value || std::is_same<T, std::complex<double>>::value)
        {
            get_lwork(lwork, work);
            work.resize(lwork, 0);
            liwork = iwork[0];
            iwork.resize(liwork, 0);
            lrwork = rwork[0] + degeneracy_max * ndim_global;
            int maxlrwork = std::max(lrwork, 3);
            rwork.resize(maxlrwork, 0);
        }
        else
        {
            get_lwork(lwork, work);
            work.resize(std::max(lwork, 3), 0);
            liwork = iwork[0];
            iwork.resize(liwork, 0);
        }

        pxxxgvx_(&itype,
                 &jobz,
                 &range,
                 &uplo,
                 &ndim_global,
                 h_tmp.data(),
                 &one,
                 &one,
                 desc,
                 s_tmp.data(),
                 &one,
                 &one,
                 desc,
                 &vl,
                 &vu,
                 &il,
                 &iu,
                 &abstol,
                 &M,
                 &NZ,
                 ekb,
                 &orfac,
                 wfc_2d,
                 &one,
                 &one,
                 desc,
                 work.data(),
                 &lwork,
                 rwork.data(), // is not used for real data type
                 &lrwork,      // is not used for real data type
                 iwork.data(),
                 &liwork,
                 ifail.data(),
                 iclustr.data(),
                 gap.data(),
                 &info);

        if (info == 0)
        {
            return;
        }
        pxxxgvx_post_processing(info, ifail,iclustr, M,NZ,nbands,degeneracy_max);

        // break the loop for real data type
        if (std::is_same<T, float>::value || std::is_same<T,double>::value)
        {
            return;
        }
    }
}

// template instantiation
template void pxxxgvx_diag(const int *const desc,
                           const int ncol,
                           const int nrow,
                           const int nbands,
                           const double *const h_mat,
                           const double *const s_mat,
                           double *const ekb,
                           double *const wfc_2d);
template void pxxxgvx_diag(const int *const desc,
                           const int ncol,
                           const int nrow,
                           const int nbands,
                           const std::complex<double> *const h_mat,
                           const std::complex<double> *const s_mat,
                           double *const ekb,
                           std::complex<double> *const wfc_2d);
template void pxxxgvx_diag(const int *const desc,
                            const int ncol,
                            const int nrow,
                            const int nbands,
                            const float *const h_mat,
                            const float *const s_mat,
                            float *const ekb,
                            float *const wfc_2d);
template void pxxxgvx_diag(const int *const desc,
                            const int ncol,
                            const int nrow,
                            const int nbands,
                            const std::complex<float> *const h_mat,
                            const std::complex<float> *const s_mat,
                            float *const ekb,
                            std::complex<float> *const wfc_2d);                           

#endif

} // namespace