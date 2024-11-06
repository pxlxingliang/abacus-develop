#ifndef TEST_DIAGO_PXXXGVX_H
#define TEST_DIAGO_PXXXGVX_H

#include "../diag_hs_para.h"
#include <gtest/gtest.h>
#include <vector>
#include <complex>
#include <random>
#include <mpi.h>

template <typename T>
typename std::enable_if<std::is_same<T, double>::value || std::is_same<T, float>::value>::type
generate_random_hs_impl(int d, std::mt19937& gen, std::uniform_real_distribution<typename GetTypeReal<T>::type>& dis, std::vector<T>& h_mat, std::vector<T>& s_mat) {
    // For S matrix, firstly we generate a random symmetric matrix s_tmp, then we set S = s_tmp * s_tmp^T + n * I
    std::vector<T> s_tmp(d*d);
    for (int i = 0; i < d; ++i) {
        for (int j = i; j < d; ++j) {
            typename GetTypeReal<T>::type value1 = static_cast<typename GetTypeReal<T>::type>(dis(gen));
            h_mat[i * d + j] = value1;
            h_mat[j * d + i] = value1;

            // construct a random overlap matrix
            typename GetTypeReal<T>::type value2 = static_cast<typename GetTypeReal<T>::type>(dis(gen));
            s_tmp[i * d + j] = value2;
            s_tmp[j * d + i] = value2;
        }
    }

    // set S = s_tmp * s_tmp^T + n * I
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            s_mat[i * d + j] = 0;
            for (int k = 0; k < d; ++k) {
                s_mat[i * d + j] += s_tmp[i * d + k] * s_tmp[j * d + k];
            }
            if (i == j) {
                s_mat[i * d + j] += 2.0;
            }
        }
    }
}

template <typename T>
typename std::enable_if<std::is_same<T, std::complex<double>>::value || std::is_same<T, std::complex<float>>::value>::type
generate_random_hs_impl(int d, std::mt19937& gen, std::uniform_real_distribution<typename GetTypeReal<T>::type>& dis, std::vector<T>& h_mat, std::vector<T>& s_mat) {
    std::vector<T> s_tmp(d*d);
    for (int i = 0; i < d; ++i) {
        for (int j = i; j < d; ++j) {
            typename GetTypeReal<T>::type value1 = static_cast<typename GetTypeReal<T>::type>(dis(gen));
            typename GetTypeReal<T>::type value2 = static_cast<typename GetTypeReal<T>::type>(dis(gen));
            h_mat[i * d + j] = T(value1, value2);
            if (i != j)
            {
                h_mat[j * d + i] = T(value1, -value2);
            }
            else{
                h_mat[j * d + i] = T(value1, 0);
            }

            // construct a random overlap matrix
            value1 = static_cast<typename GetTypeReal<T>::type>(dis(gen));
            value2 = static_cast<typename GetTypeReal<T>::type>(dis(gen));
            s_tmp[i * d + j] = T(value1, value2);
            s_tmp[j * d + i] = T(value1, -value2);
        }
    }

    // set S = s_tmp * s_tmp^T + n * I
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            s_mat[i * d + j] = T(0, 0);
            for (int k = 0; k < d; ++k) {
                s_mat[i * d + j] += s_tmp[i * d + k] * std::conj(s_tmp[j * d + k]);
            }
            if (i == j) {
                s_mat[i * d + j] += T(2.0, 0);
            }
        }
    }
}


template <typename T>
void generate_random_hs(int d, int random_seed ,std::vector<T>& h_mat, std::vector<T>& s_mat) {
    std::mt19937 gen(random_seed);
    std::uniform_real_distribution<typename GetTypeReal<T>::type> dis(-1.0,1.0);

    h_mat.resize(d * d);
    s_mat.resize(d * d);
    generate_random_hs_impl(d, gen, dis, h_mat, s_mat);
}


template <typename T>
typename std::enable_if<std::is_same<T, double>::value || std::is_same<T, float>::value>::type
verify_results(const std::vector<T>& h_psi, const std::vector<T>& s_psi, const std::vector<typename GetTypeReal<T>::type>& ekb, int lda, int nbands, double threshold) {
    for (int i = 0; i < lda; ++i) {
        for (int j = 0; j < nbands; ++j) {
            ASSERT_NEAR(h_psi[j * lda + i], ekb[j] * s_psi[j * lda + i], threshold);
        }
    }
}

template <typename T>
typename std::enable_if<std::is_same<T, std::complex<double>>::value || std::is_same<T, std::complex<float>>::value>::type
verify_results(const std::vector<T>& h_psi, const std::vector<T>& s_psi, const std::vector<typename GetTypeReal<T>::type>& ekb, int lda, int nbands, double threshold) {
    for (int i = 0; i < lda; ++i) {
        for (int j = 0; j < nbands; ++j) {
            ASSERT_NEAR(h_psi[j * lda + i].real(), ekb[j] * s_psi[j * lda + i].real(), threshold);
            ASSERT_NEAR(h_psi[j * lda + i].imag(), ekb[j] * s_psi[j * lda + i].imag(), threshold);
        }
    }
}

template <typename T>
void test_diago_hs(int lda, int nb, int random_seed, int nbands, int diag_type, MPI_Comm comm) {
    // diag_type should be 1 (for elpa) or 2 (for scalapack)
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);

    std::vector<T> h_mat, s_mat, wfc, h_psi, s_psi;
    std::vector<typename GetTypeReal<T>::type> ekb(lda);
    if (my_rank==0)
    {
        h_mat.resize(lda * lda);
        s_mat.resize(lda * lda);
        wfc.resize(lda * lda);
        generate_random_hs(lda, random_seed, h_mat, s_mat);
    }
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    hsolver::Diago_HS_para<T>(h_mat.data(), s_mat.data(), lda, nbands,ekb.data(), wfc.data(), comm, nb, diag_type);
    MPI_Barrier(comm);
    std::cout << __FILE__ << " " << __LINE__ << std::endl;

    // Verify results
    if (my_rank == 0){
        double threshold = 1e-6;
        if (std::is_same<T, std::complex<double>>::value || std::is_same<T, double>::value) {
            threshold = 1e-12;
        }

        h_psi.resize(lda * nbands, 0);
        s_psi.resize(lda * nbands, 0);

        for (int i = 0; i < lda; ++i) {
            for (int j = 0; j < nbands; ++j) {
                for (int k = 0; k < lda; ++k) {
                    h_psi[j * lda + i] += h_mat[k * lda + i] * wfc[j * lda + k];
                    s_psi[j * lda + i] += s_mat[k * lda + i] * wfc[j * lda + k];
                }
            }
        }
        verify_results<T>(h_psi, s_psi, ekb, lda, nbands, threshold);
    }
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
}

//test_diago_hs(int lda, int nb, int random_seed, int nbands, int diag_type, MPI_Comm comm)
//TEST(DiagoPxxxgvxElpaTest, Double) {
//    test_diago_hs<double>(16, 4, 0, 10, 1,MPI_COMM_WORLD);
//}
//
//TEST(DiagoPxxxgvxElpaTest, ComplexDouble) {
//    test_diago_hs<std::complex<double>>(16, 4, 0, 10, 1, MPI_COMM_WORLD);
//}

TEST(DiagoPxxxgvxScalapackTest, Double) {

    test_diago_hs<double>(16, 4, 0, 10, 2,MPI_COMM_WORLD);
}

//TEST(DiagoPxxxgvxScalapackTest, ComplexDouble) {
//    test_diago_hs<std::complex<double>>(16, 4, 0, 10, 2, MPI_COMM_WORLD);
//}

//TEST(DiagoPxxxgvxScalapackTest, Float) {
//    test_diago_hs<float>(16, 4, 0, 10,2,MPI_COMM_WORLD);
//}
//
//TEST(DiagoPxxxgvxScalapackTest, ComplexFloat) {
//    test_diago_hs<std::complex<float>>(16, 4, 0, 10,2,MPI_COMM_WORLD);
//}

int main(int argc, char** argv) {
    //::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

    if (myrank != 0) {
        //::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
        delete listeners.Release(listeners.default_result_printer());
    }

    int result = RUN_ALL_TESTS();
    if (myrank == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        MPI_Finalize();
        return result;
	}

    MPI_Finalize();

	return 0;
}

#endif // TEST_DIAGO_PXXXGVX_H