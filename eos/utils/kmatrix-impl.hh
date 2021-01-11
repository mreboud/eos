#ifndef EOS_GUARD_EOS_UTILS_KMATRIX_IMPL_HH
#define EOS_GUARD_EOS_UTILS_KMATRIX_IMPL_HH 1

#include <eos/utils/kmatrix.hh>
#include <eos/utils/complex.hh>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

//#include <iostream>

namespace eos

{
    template <unsigned nchannels_, unsigned nresonances_>
    KMatrix<nchannels_, nresonances_>::KMatrix(std::initializer_list<std::shared_ptr<KMatrix::Channel>> channels,
					       std::initializer_list<std::shared_ptr<KMatrix::Resonance>> resonances,
					       const std::string & prefix) :
	_channels(channels),
	_resonances(resonances),
	_prefix(prefix),
	_Khat(gsl_matrix_complex_calloc(nchannels_, nchannels_)),
	_That(gsl_matrix_complex_calloc(nchannels_, nchannels_)),
	_tmp_1(gsl_matrix_complex_calloc(nchannels_, nchannels_)),
	_tmp_2(gsl_matrix_complex_calloc(nchannels_, nchannels_)),
	_perm(gsl_permutation_calloc(nchannels_))
    {
	//Perform size checks
	if (channels.size() != nchannels_)
	    {
		throw InternalError("The size of the channels array does not match nchannels_");
	    }
	if (resonances.size() != nresonances_)
	    {
		throw InternalError("The size of the resonances array does not match nresonances_");
	    }


	//Perform pointer checks
	if (_Khat == nullptr)
	    {
		throw std::bad_alloc();
	    }
	if (_That == nullptr)
	    {
		throw std::bad_alloc();
	    }
	if (_tmp_1 == nullptr)
	    {
		throw std::bad_alloc();
	    }
	if (_tmp_2 == nullptr)
	    {
		throw std::bad_alloc();
	    }
	if (_perm == nullptr)
	    {
		throw std::bad_alloc();
	    }

    }

    //Return row n° rowindex of the T matrix defined as T = (1-i*rho*K)^(-1)*K
    template <unsigned nchannels_, unsigned nresonances_>
    std::array<complex<double>, nchannels_>
    KMatrix<nchannels_, nresonances_>::tmatrix_row(unsigned rowindex, double s) const
    {

	std::array<complex<double>, nchannels_> tmatrixrow;
	auto channels = this->_channels;
	auto resonances = this->_resonances;


	// std::cout << "Entering tmatrix row evaluation, nbr channels = " << channels.size() << " nbr resonances = " << resonances.size() << std::endl;
	// std::cout << " --- Working at s = " << s << std::endl;


	///////////////////
	//1. Fill tmp1 = rho matrix
	///////////////////

	// zero out tmp_1;
	gsl_matrix_complex_set_zero(_tmp_1);

	for (size_t i = 0 ; i < nchannels_ ; ++i)
	    {
		complex<double> rhoentry = channels[i]->rho(s);
		// std::cout << " --- Channel "<< i << " m1 " << channels[i]->_m1 << " m2 " << channels[i]->_m2 << std::endl;
		// std::cout << " ------" << " rho(s) = " << rhoentry.real() << " + i*" <<  rhoentry.imag() << std::endl;
		gsl_matrix_complex_set(_tmp_1,  i,  i, gsl_complex_rect(rhoentry.real(), rhoentry.imag()));
	    }

	///////////////////
	// 2. Fill Khat
	///////////////////

	// row index
	for (size_t i = 0 ; i < nchannels_ ; ++i)
	    {
		// column index
		for (size_t j = 0 ; j < nchannels_ ; ++j)
		    {

			complex<double> entry = 0.0;

			// resonance index
			for (size_t a = 0 ; a < nresonances_ ; ++a)
			    {

				double mres = resonances[a]->_m;

				Parameter g0rci = channels[i]->_g0s[a];
				Parameter g0rcj = channels[j]->_g0s[a];

				// std::cout << " --- Resonance "<< a << " m = " << mres << " g0rci = " << g0rci << " g0rcj = " << g0rcj << std::endl;

				// DO SOMETHING TO AVOID THE POLES (DIVISION BY 0)
				// * A first solution is to add an small imaginary part
				// -> It may bias the fit and forces Khat to be complex...
				complex<double> mres2 = { mres * mres, 1.0e-7 };
				entry += g0rci * g0rcj / (mres2 - s);
				// * One could also avoid this by averaging s-eps and s+eps at the pole
				// TODO...

			    }

			// std::cout << " --- Khat "<< i << " " << j << " = " << entry << std::endl;
			gsl_matrix_complex_set(_Khat, i, j, gsl_complex_rect(entry.real(), entry.imag()));
		    }
	    }
	// std::cout << " --- KMatrix computed " << std::endl;



	///////////////////
	// 3. Compute That
	///////////////////

	// std::cout << " --- Entering KMatrix inversion: " << std::endl;

	// For the moment I am computing the full T matrix, while we only need one row
	static const gsl_complex one  = gsl_complex_rect(1.0, 0.0);
	static const gsl_complex minusi  = gsl_complex_rect(0.0, -1.0);
	static const gsl_complex zero = gsl_complex_rect(0.0, 0.0);

	// 3a. multiply -i*rho with Khat from the right
	// -> tmp_2 = -i * tmp_1 * Khat + zero * tmp_2 ; no transpose anywhere
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, minusi, _tmp_1, _Khat, zero, _tmp_2);

	// 3b. add 1.0 to the diagonal elements of tmp_2
	for (unsigned i = 0 ; i < nchannels_ ; ++i)
	    {
		gsl_complex value = gsl_matrix_complex_get(_tmp_2, i, i);
		// std::cout << " --- (1-i*rho*K)_{" << i << "," << i << "} = " << GSL_REAL(gsl_complex_add_real(value, 1.0)) << " + i*" << GSL_IMAG(gsl_complex_add_real(value, 1.0)) << std::endl;
		gsl_matrix_complex_set(_tmp_2, i, i, gsl_complex_add_real(value, 1.0));
	    }

	// 3c. invert tmp_2 and assign the result to tmp_1
	int signum = 0;
	gsl_permutation_init(_perm);
	gsl_linalg_complex_LU_decomp(_tmp_2, _perm, &signum);
	gsl_linalg_complex_LU_invert(_tmp_2, _perm, _tmp_1);

	// 3d. calculate That = Khat / (1 - i rho * Khat) = Khat * tmp_1
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, _Khat, _tmp_1, zero, _That);


	///////////////////
	// 4. Extract T matrix row
	///////////////////

	for (size_t i = 0 ; i < nchannels_ ; ++i)
	    {
		gsl_complex value = gsl_matrix_complex_get(_That, rowindex, i);
		complex<double> cvalue(GSL_REAL(value), GSL_IMAG(value));
		tmatrixrow[i] = cvalue;
		// std::cout << " --- That_{"<< rowindex << "," << i << "} = " << GSL_REAL(value) << " + i*" << GSL_IMAG(value) << std::endl << std::endl;
	    }
	// std::cout << std::endl;
	return tmatrixrow;
    }


    template <unsigned nchannels_, unsigned nresonances_>
    KMatrix<nchannels_, nresonances_>::~KMatrix()
    {
	// unalloc gsl...

	if (_perm)
	    gsl_permutation_free(_perm);
	_perm = nullptr;

	if (_tmp_1)
	    gsl_matrix_complex_free(_tmp_1);
	_tmp_1 = nullptr;

	if (_tmp_2)
	    gsl_matrix_complex_free(_tmp_2);
	_tmp_2 = nullptr;

	if (_That)
	    gsl_matrix_complex_free(_That);
	_That = nullptr;

	if (_Khat)
	    gsl_matrix_complex_free(_Khat);
	_Khat = nullptr;
    }

}

#endif