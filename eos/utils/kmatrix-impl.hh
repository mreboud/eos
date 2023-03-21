/*
 * Copyright (c) 2019 Stephan Kürten
 * Copyright (c) 2019 Danny van Dyk
 * Copyright (c) 2021 Méril Reboud
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef EOS_GUARD_EOS_UTILS_KMATRIX_IMPL_HH
#define EOS_GUARD_EOS_UTILS_KMATRIX_IMPL_HH 1

#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/kmatrix.hh>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include <limits>

namespace eos
{
    using std::abs;
    using std::sqrt;

    template <unsigned nchannels_, unsigned nresonances_, unsigned order_>
    KMatrix<nchannels_, nresonances_, order_>::KMatrix(std::array<std::shared_ptr<KMatrix::Channel>, nchannels_> channels,
                                                       std::array<std::shared_ptr<KMatrix::Resonance>, nresonances_> resonances,
                                                       std::array<std::array<std::array<Parameter, nchannels_>, nchannels_>, order_ + 1> bkgpol,
                                                       const std::string & prefix) :
        _channels(channels),
        _resonances(resonances),
        _bkgpol(bkgpol),
        _prefix(prefix),
        _Khat(gsl_matrix_complex_calloc(nchannels_, nchannels_)),
        _That(gsl_matrix_complex_calloc(nchannels_, nchannels_)),
        _tmp_1(gsl_matrix_complex_calloc(nchannels_, nchannels_)),
        _tmp_2(gsl_matrix_complex_calloc(nchannels_, nchannels_)),
        _tmp_3(gsl_matrix_complex_calloc(nchannels_, nchannels_)),
        _perm(gsl_permutation_calloc(nchannels_))
    {
        // Perform size checks
        if (channels.size() != nchannels_)
            throw InternalError("The size of the channels array does not match nchannels_.");
        if (bkgpol.size() != order_ + 1)
            throw InternalError("The size of the array of background constants does not match order_ + 1.");
        if (bkgpol.size() != 0 and bkgpol[0].size() != nchannels_)
            throw InternalError("The array of background constants is not valid");
        if (bkgpol.size() != 0 and bkgpol[0].size() != 0 and bkgpol[0][0].size() != nchannels_)
            throw InternalError("The array of background constants is not valid");
        if (resonances.size() != nresonances_)
            throw InternalError("The size of the resonances array does not match nresonances_.");

        // Perform pointer checks
        if (_Khat == nullptr)
            throw std::bad_alloc();
        if (_That == nullptr)
            throw std::bad_alloc();
        if (_tmp_1 == nullptr)
            throw std::bad_alloc();
        if (_tmp_2 == nullptr)
            throw std::bad_alloc();
        if (_tmp_3 == nullptr)
            throw std::bad_alloc();
        if (_perm == nullptr)
            throw std::bad_alloc();
    }

    template <unsigned nchannels_, unsigned nresonances_, unsigned order_>
    KMatrix<nchannels_, nresonances_, order_>::~KMatrix()
    {
        if (_perm)
            gsl_permutation_free(_perm);
        _perm = nullptr;

        if (_tmp_1)
            gsl_matrix_complex_free(_tmp_1);
        _tmp_1 = nullptr;

        if (_tmp_2)
            gsl_matrix_complex_free(_tmp_2);
        _tmp_2 = nullptr;

        if (_tmp_3)
            gsl_matrix_complex_free(_tmp_3);
        _tmp_3 = nullptr;

        if (_That)
            gsl_matrix_complex_free(_That);
        _That = nullptr;

        if (_Khat)
            gsl_matrix_complex_free(_Khat);
        _Khat = nullptr;
    }


    // Adapt s to avoid resonances masses
    template <unsigned nchannels_, unsigned nresonances_, unsigned order_>
    complex<double>
    KMatrix<nchannels_, nresonances_, order_>::adapt_s(const complex<double> s) const
    {
        // Disallowed range around resonance masses
        const double minimal_distance = 1.0e-7;
        const auto & resonances = this->_resonances;

        complex<double> adapted_s = s;

        for (size_t a = 0 ; a < nresonances_ ; ++a)
        {
            const double mres_a = resonances[a]->_m;

            for (size_t b = 0 ; b < a ; ++b)
            {
                const double mres_b = resonances[b]->_m;

                if (abs(mres_a * mres_a - mres_b * mres_b) < minimal_distance)
                    throw InternalError("The resonances masses are degenerate.");
            }

            if (abs(mres_a * mres_a - s) < minimal_distance) [[unlikely]]
            {
                if (real(s) > mres_a * mres_a)
                    adapted_s = mres_a * mres_a + minimal_distance;
                else
                    adapted_s = mres_a * mres_a - minimal_distance;
            }
        }

        return adapted_s;
    }


    // Return the row corresponding to the index rowindex of the T matrix defined as T = (1-i*rho*K)^(-1)*K
    template <unsigned nchannels_, unsigned nresonances_, unsigned order_>
    std::array<complex<double>, nchannels_>
    KMatrix<nchannels_, nresonances_, order_>::tmatrix_row(unsigned rowindex, complex<double> _s) const
    {
        std::array<complex<double>, nchannels_> tmatrixrow;
        const auto & channels = this->_channels;
        const auto & resonances = this->_resonances;
        const auto & bkgpol = this->_bkgpol;
        // Adapt s to avoid pole in the K matrix
        const complex<double> s = adapt_s(_s);

        ///////////////////
        // 1. Fill tmp1 = rho matrix
        ///////////////////
        gsl_matrix_complex_set_zero(_tmp_1);
        gsl_matrix_complex_set_zero(_tmp_2);
        gsl_matrix_complex_set_zero(_tmp_3);

        for (size_t i = 0 ; i < nchannels_ ; ++i)
        {
            complex<double> rhoentry = channels[i]->rho(s);
            gsl_matrix_complex_set(_tmp_1,  i,  i, gsl_complex_rect(rhoentry.real(), rhoentry.imag()));
        }

        ///////////////////
        // 2. Fill tmp2 =  n matrix
        ///////////////////
        for (size_t i = 0 ; i < nchannels_ ; ++i)
        {
            const unsigned li = channels[i]->_l_orbital;

            // TODO: the BW-Factors should be independent of the resonance, shouldn't they?
            // If I want to make them dependent on the resonance like before I need to access the information about the
            // resonances here. That is not intended in the PDG, though.

            //const double qres = resonances[a]->_q;
            const double qres = 2; //I think the breakup momentum for the blattweisskopf needs to be universal here?

            //const complex<double> mres_2 = power_of<2>((double)resonances[a]->_m);
            const complex<double> mres_2 = 4;

            const complex<double> mi1_2 = power_of<2>((double)channels[i]->_m1);
            const complex<double> mi2_2 = power_of<2>((double)channels[i]->_m2);

            const double BFi = kmatrix_utils::blatt_weisskopf_factor(li, 0.5 * sqrt(abs(eos::lambda(s, mi1_2, mi2_2) / s)), qres) /
                    kmatrix_utils::blatt_weisskopf_factor(li, 0.5 * sqrt(abs(eos::lambda(mres_2, mi1_2, mi2_2) / mres_2)), qres);

            //TODO: I can not put li in the eos power_of Function. What to do here?
            //const double QFi = power_of<li>(0.5 * sqrt(abs(eos::lambda(mres_2, mi1_2, mi2_2) / mres_2)) / qres);
            const double QFi = power_of<2>(0.5 * sqrt(abs(eos::lambda(mres_2, mi1_2, mi2_2) / mres_2)) / qres);

            complex<double> nentry = QFi * BFi;
            gsl_matrix_complex_set(_tmp_2,  i,  i, gsl_complex_rect(nentry.real(), nentry.imag()));
        }

        ///////////////////
        // 3. Fill Khat
        ///////////////////
        for (size_t i = 0 ; i < nchannels_ ; ++i)
        {
            for (size_t j = 0 ; j < nchannels_ ; ++j)
            {
                complex<double> entry = 0.;

                // Fast evaluation of the polynomial describing the non-resonant contribution
                for (size_t k = 0 ; k <= order_ ; ++k)
                {
                    // The square root is evaluated with a branch cut along the negative real axis
                    entry = 1.0 / sqrt(s) * entry + complex<double>(bkgpol[order_ - k][i][j].evaluate(), 0.0);
                }

                for (size_t a = 0 ; a < nresonances_ ; ++a)
                {
                    const complex<double> mres_2 = power_of<2>((double)resonances[a]->_m);

                    Parameter g0rci = channels[i]->_g0s[a];
                    Parameter g0rcj = channels[j]->_g0s[a];

                    entry += g0rci * g0rcj / (mres_2 - s);
                }

                gsl_matrix_complex_set(_Khat, i, j, gsl_complex_rect(entry.real(), entry.imag()));
            }
        }

        ///////////////////
        // 4. Compute That
        ///////////////////
        static const gsl_complex one  = gsl_complex_rect(1.0, 0.0);
        static const gsl_complex minusi  = gsl_complex_rect(0.0, -1.0);
        static const gsl_complex zero = gsl_complex_rect(0.0, 0.0);

        // 4a. multiply Khat with -i*rho and assign to tmp_3
        // -> tmp_3 = -i * Khat * tmp_1; tmp_1=rho
        // This is the only time I need tmp_1 for the rho matrix
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, minusi, _Khat, _tmp_1, zero, _tmp_3);

        // 4b. Set 1.0 to the diagonal elements of tmp_1
        // tmp_1 = 1
        for (unsigned i = 0 ; i < nchannels_ ; ++i)
        {
            gsl_matrix_complex_set(_tmp_1, i, i, one);
        }

        // 4c. multiply tmp_3 with n
        // -> tmp_2 = -i * Khat * rho * n 
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, _tmp_3, _tmp_2, zero, _tmp_3);

        // 4d. multiply tmp_3 with n, add unit matrix and assign to tmp_1
        // -> tmp_1 = -i * Khat * rho * n^2 + 1 
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, _tmp_3, _tmp_2, one, _tmp_1);

        // 4e. invert tmp_1 and assign the result to tmp_3
        // -> tmp_3 = (tmp_1)^-1
        int signum = 0;
        gsl_permutation_init(_perm);
        gsl_linalg_complex_LU_decomp(_tmp_1, _perm, &signum);
        gsl_linalg_complex_LU_invert(_tmp_1, _perm, _tmp_3);

        // 4f. multiply tmp_2 (n-matrix) with tmp_3 and assign to tmp_3
        // -> tmp_3 = n * tmp_3
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, _tmp_2, _tmp_3, zero, _tmp_3);

        // 4g. multiply tmp_3 with Khat and assign to tmp_3
        // -> tmp_3 = tmp_3 * Khat
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, _tmp_3, _Khat, zero, _tmp_3);

        // 4h. multiply tmp_3 with tmp_2 (n-matrix) and assign to That
        // -> tmp_3 = tmp_3 * n
        // now the T matrix is That = n * Inverse[ 1- Khat * i * rho *  n * n] * Khat * n
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, _tmp_3, _tmp_2, zero, _That);

        ///////////////////
        // 4. Extract T matrix row
        ///////////////////
        for (size_t i = 0 ; i < nchannels_ ; ++i)
        {
            gsl_complex value = gsl_matrix_complex_get(_That, rowindex, i);
            complex<double> cvalue(GSL_REAL(value), GSL_IMAG(value));
            tmatrixrow[i] = cvalue;
        }

        return tmatrixrow;
    }

    template <unsigned nchannels_, unsigned nresonances_, unsigned order_>
    double
    KMatrix<nchannels_, nresonances_, order_>::partial_width(unsigned resonance, unsigned channel) const
    {
        double mres = this->_resonances[resonance]->_m;
        double rho  = real(this->_channels[channel]->rho(mres*mres));

        return power_of<2>((double)this->_channels[channel]->_g0s[resonance]) / mres * rho;
    }

    template <unsigned nchannels_, unsigned nresonances_, unsigned order_>
    double
    KMatrix<nchannels_, nresonances_, order_>::width(unsigned resonance) const
    {
        double result = 0.0;

        for (size_t i = 0; i < nchannels_; i++)
        {
            result += KMatrix<nchannels_, nresonances_, order_>::partial_width(resonance, i);
        }

        return result;
    }
}

#endif
