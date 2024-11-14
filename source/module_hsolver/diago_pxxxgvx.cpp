#include <complex>
#include <vector>
#include "string.h"
#include <iostream>

#include "diago_pxxxgvx.h"
#include "module_base/scalapack_connector.h"
#include "module_base/blacs_connector.h"

namespace hsolver
{

#ifdef __MPI
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

void pxxxgvx_diag(const int* const desc,
                  const int ncol,
                  const int nrow,
                  const int nbands,
                  const double* const h_mat,
                  const double* const s_mat,
                  double* const ekb,
                  double* const wfc_2d)
{
    std::vector<double> h_tmp(ncol * nrow);
    std::vector<double> s_tmp(ncol * nrow);
    memcpy(h_tmp.data(), h_mat, sizeof(double) * ncol * nrow);
    memcpy(s_tmp.data(), s_mat, sizeof(double) * ncol * nrow);

    int ndim_global = desc[2];
    int nprow, npcol, myprow, mypcol;
    Cblacs_gridinfo(desc[1], &nprow, &npcol, &myprow, &mypcol);
    int dsize = nprow * npcol;

    const char jobz = 'V', range = 'I', uplo = 'U';
    const int itype = 1, il = 1, iu = nbands, one = 1;
    int M = 0, NZ = 0, lwork = -1, liwork = -1, info = 0;
    double vl = 0, vu = 0;
    const double abstol = 0, orfac = -1;
    std::vector<double> work(3, 0);
    std::vector<int> iwork(1, 0);
    std::vector<int> ifail(ndim_global, 0);
    std::vector<int> iclustr(2 * dsize);
    std::vector<double> gap(dsize);

    pdsygvx_(&itype,
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

    lwork = work[0];
    work.resize(std::max(lwork, 3), 0);
    liwork = iwork[0];
    iwork.resize(liwork, 0);

    pdsygvx_(&itype,
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
             iwork.data(),
             &liwork,
             ifail.data(),
             iclustr.data(),
             gap.data(),
             &info);


    if (info == 0)
    {
        return ;
    }
    int degeneracy_max = 12;
    pxxxgvx_post_processing(info, ifail,iclustr, M,NZ,nbands,degeneracy_max);
}

void pxxxgvx_diag(const int* const desc,
                  const int ncol,
                  const int nrow,
                  const int nbands,
                  const std::complex<double>* const h_mat,
                  const std::complex<double>* const s_mat,
                  double* const ekb,
                  std::complex<double>* const wfc_2d)
{
    int degeneracy_max = 12;

    int nprow, npcol, myprow, mypcol;
    Cblacs_gridinfo(desc[1], &nprow, &npcol, &myprow, &mypcol);
    int dsize = nprow * npcol;

    while (true)
    {
        std::vector<std::complex<double>> h_tmp(ncol * nrow);
        std::vector<std::complex<double>> s_tmp(ncol * nrow);
        memcpy(h_tmp.data(), h_mat, sizeof(std::complex<double>) * ncol * nrow);
        memcpy(s_tmp.data(), s_mat, sizeof(std::complex<double>) * ncol * nrow);

        int ndim_global = desc[2];
        const char jobz = 'V', range = 'I', uplo = 'U';
        const int itype = 1, il = 1, iu = nbands, one = 1;
        int M = 0, NZ = 0, lwork = -1, lrwork = -1, liwork = -1, info = 0;
        const double abstol = 0, orfac = -1;
        const double vl = 0, vu = 0;
        std::vector<std::complex<double>> work(1, 0);
        std::vector<double> rwork(3, 0);
        std::vector<int> iwork(1, 0);
        std::vector<int> ifail(ndim_global, 0);
        std::vector<int> iclustr(2 * dsize);
        std::vector<double> gap(dsize);

        pzhegvx_(&itype,
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
                 rwork.data(),
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

        lwork = work[0].real();
        work.resize(lwork, 0);
        lrwork = rwork[0] + degeneracy_max * ndim_global;

        int maxlrwork = std::max(lrwork, 3);
        rwork.resize(maxlrwork, 0);
        liwork = iwork[0];
        iwork.resize(liwork, 0);

        pzhegvx_(&itype,
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
                 rwork.data(),
                 &lrwork,
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
    }
}


void pxxxgvx_diag(const int *const desc,
                      const int ncol,
                      const int nrow,
                      const int nbands,
                      const float *const h_mat,
                      const float *const s_mat,
                      float *const ekb,
                      float *const wfc_2d)
{
    int nprow, npcol, myprow, mypcol;
    Cblacs_gridinfo(desc[1], &nprow, &npcol, &myprow, &mypcol);
    int dsize = nprow * npcol;

    std::vector<float> h_tmp(ncol * nrow);
    std::vector<float> s_tmp(ncol * nrow);
    memcpy(h_tmp.data(), h_mat, sizeof(float) * ncol * nrow);
    memcpy(s_tmp.data(), s_mat, sizeof(float) * ncol * nrow);

    int ndim_global = desc[2];
    const char jobz = 'V', range = 'I', uplo = 'U';
    const int itype = 1, il = 1, iu = nbands, one = 1;
    int M = 0, NZ = 0, lwork = -1, liwork = -1, info = 0;
    float vl = 0, vu = 0;
    const float abstol = 0, orfac = -1;
    std::vector<float> work(3, 0);
    std::vector<int> iwork(1, 0);
    std::vector<int> ifail(ndim_global, 0);
    std::vector<int> iclustr(2 * dsize);
    std::vector<float> gap(dsize);

    pssygvx_(&itype,
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

    lwork = work[0];
    work.resize(std::max(lwork, 3), 0);
    liwork = iwork[0];
    iwork.resize(liwork, 0);

    pssygvx_(&itype,
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
             iwork.data(),
             &liwork,
             ifail.data(),
             iclustr.data(),
             gap.data(),
             &info);


    if (info == 0)
    {
        return ;
    }
    int degeneracy_max = 12;
    pxxxgvx_post_processing(info, ifail,iclustr, M,NZ,nbands,degeneracy_max);
}    

void pxxxgvx_diag(const int *const desc,
                      const int ncol,
                      const int nrow,
                      const int nbands,
                      const std::complex<float> *const h_mat,
                      const std::complex<float> *const s_mat,
                      float *const ekb,
                      std::complex<float> *const wfc_2d)
{
    int nprow, npcol, myprow, mypcol;
    Cblacs_gridinfo(desc[1], &nprow, &npcol, &myprow, &mypcol);
    int dsize = nprow * npcol;

    int degeneracy_max = 12;
    while (true)
    {
        std::vector<std::complex<float>> h_tmp(ncol * nrow);
        std::vector<std::complex<float>> s_tmp(ncol * nrow);
        memcpy(h_tmp.data(), h_mat, sizeof(std::complex<float>) * ncol * nrow);
        memcpy(s_tmp.data(), s_mat, sizeof(std::complex<float>) * ncol * nrow);

        int ndim_global = desc[2];
        const char jobz = 'V', range = 'I', uplo = 'U';
        const int itype = 1, il = 1, iu = nbands, one = 1;
        int M = 0, NZ = 0, lwork = -1, lrwork = -1, liwork = -1, info = 0;
        const float abstol = 0, orfac = -1;
        const float vl = 0, vu = 0;
        std::vector<std::complex<float>> work(1, 0);
        std::vector<float> rwork(3, 0);
        std::vector<int> iwork(1, 0);
        std::vector<int> ifail(ndim_global, 0);
        std::vector<int> iclustr(2 * dsize);
        std::vector<float> gap(dsize);

        pchegvx_(&itype,
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
                    rwork.data(),
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

        lwork = work[0].real();
        work.resize(lwork, 0);
        lrwork = rwork[0] + degeneracy_max * ndim_global;

        int maxlrwork = std::max(lrwork, 3);
        rwork.resize(maxlrwork, 0);
        liwork = iwork[0];
        iwork.resize(liwork, 0);

        pchegvx_(&itype,
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
                 rwork.data(),
                 &lrwork,
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
    }
}


#endif

} // namespace