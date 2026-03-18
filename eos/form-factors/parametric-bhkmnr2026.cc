/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Fatemeh Nouri
 * Copyright (c) 2026 Méril Reboud
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

#include <eos/form-factors/parametric-bhkmnr2026.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/stringify.hh>
#include <eos/maths/integrate.hh>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>


#include <functional>
#include <numeric>

namespace eos
{
    using namespace std::literals::string_literals;

    BHKMNR2026FormFactors<VacuumToPiPi>::BHKMNR2026FormFactors(const Parameters & p, const Options & o) :
        _a_fp_I1{{
            UsedParameter(p[_par_name("+", "1", "4")], *this),
            UsedParameter(p[_par_name("+", "1", "5")], *this),
            UsedParameter(p[_par_name("+", "1", "6")], *this),
            UsedParameter(p[_par_name("+", "1", "7")], *this),
            UsedParameter(p[_par_name("+", "1", "8")], *this),
            UsedParameter(p[_par_name("+", "1", "9")], *this),
            UsedParameter(p[_par_name("+", "1", "10")], *this),
            UsedParameter(p[_par_name("+", "1", "11")], *this),
            UsedParameter(p[_par_name("+", "1", "12")], *this)
        }},
        _M_fp_I1{{
            UsedParameter(p["0->pipi::M_(+,1,0)@BHKMNR2026"], *this),
            UsedParameter(p["0->pipi::M_(+,1,1)@BHKMNR2026"], *this),
            UsedParameter(p["0->pipi::M_(+,1,2)@BHKMNR2026"], *this)
        }},
        _G_fp_I1{{
            UsedParameter(p["0->pipi::Gamma_(+,1,0)@BHKMNR2026"], *this),
            UsedParameter(p["0->pipi::Gamma_(+,1,1)@BHKMNR2026"], *this),
            UsedParameter(p["0->pipi::Gamma_(+,1,2)@BHKMNR2026"], *this)
        }},
        _n_resonances(o, option_specifications, "n-resonances"_ok),
        _m_pi(p["mass::pi^+"], *this),
        _s_0(p["0->pipi::s_0@BHKMNR2026"], *this),
        _s_in(p["0->pipi::s_in@BHKMNR2026"], *this),
        _hbar(p["QM::hbar"], *this),
        _M(gsl_matrix_alloc(4, 4)),
        _inv_M(gsl_matrix_alloc(4, 4)),
        _L(gsl_vector_alloc(4)),
        _perm(gsl_permutation_calloc(4)),
        _constrained_coefficents(gsl_vector_alloc(4))

    {
        // Perform pointer checks
        if (_M == nullptr)
            throw std::bad_alloc();
        if (_inv_M == nullptr)
            throw std::bad_alloc();
        if (_L == nullptr)
            throw std::bad_alloc();
        if (_perm == nullptr)
            throw std::bad_alloc();
        if (_constrained_coefficents == nullptr)
            throw std::bad_alloc();
    }

    BHKMNR2026FormFactors<VacuumToPiPi>::~BHKMNR2026FormFactors()
    {
        if (_perm)
        {
            gsl_permutation_free(_perm);
        }
        _perm = nullptr;
        if (_constrained_coefficents)
        {
            gsl_vector_free(_constrained_coefficents);
        }
        _constrained_coefficents = nullptr;
        if (_L)
        {
            gsl_vector_free(_L);
        }
        _L = nullptr;
        if (_inv_M)
        {
            gsl_matrix_free(_inv_M);
        }
        _inv_M = nullptr;
        if (_M)
        {
            gsl_matrix_free(_M);
        }
        _M = nullptr;
    }

    FormFactors<VacuumToPP> *
    BHKMNR2026FormFactors<VacuumToPiPi>::make(const Parameters & p, const Options & o)
    {
        return new BHKMNR2026FormFactors<VacuumToPiPi>(p, o);
    }


    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::psi(const complex<double> & s) const
    {
        return this->_s_to_psi_11(s);
    }


    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::P(const complex<double> & psi) const
    {
        return this->_P(psi);
    }


    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::dPdpsi(const complex<double> & psi) const
    {
        return this->_dPdpsi(psi);
    }


    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::dfdpsi_terms(const unsigned k, const complex<double> & psi) const
    {
        return this->_dfdpsi_terms(k, psi);
    }


    complex<double> BHKMNR2026FormFactors<VacuumToPiPi>::series(const complex<double> & psi, const std::array<double, 13> & a) const
    {
        complex<double> series    = 0.0;
        complex<double> psi_power = 1.0;

        for (std::size_t k = 0; k < a.size(); ++k)
        {
            series    += a[k] * psi_power;
            psi_power *= psi;
        }

        return series;
    }


    std::array<double, 4u>
    BHKMNR2026FormFactors<VacuumToPiPi>::constrained_a_fp_I1() const
    {

        const complex<double> psi_p  = _s_to_psi_11(_s_p());
        const complex<double> psi_in = _s_to_psi_11(_s_in());
        const complex<double> psi_0  = _s_to_psi_11(_s_0());
        const complex<double> P0     = _P(psi_0);

        //Fill M
        complex<double> psi0_pow = 1.0;

        for (unsigned k = 0; k < 4; ++k)
        {
            const complex<double> val_p  = _dfdpsi_terms(k, psi_p);
            const complex<double> val_in = _dfdpsi_terms(k, psi_in);

            gsl_matrix_set(_M, 0, k, val_p.real());
            gsl_matrix_set(_M, 1, k, val_in.real());
            gsl_matrix_set(_M, 2, k, val_in.imag());
            gsl_matrix_set(_M, 3, k, (P0 * psi0_pow).real());

            psi0_pow *= psi_0;
        }

        //Fill L
        const unsigned n = 3 + _a_fp_I1.size();

        complex<double> sum_p    = 0.0;
        complex<double> sum_in   = 0.0;
        complex<double> sum_0    = 0.0;

        for (unsigned k = 4; k <= n; ++k)
        {
            const double ak = _a_fp_I1[k - 4]();

            sum_p  += ak * _dfdpsi_terms(k, psi_p);
            sum_in += ak * _dfdpsi_terms(k, psi_in);
            sum_0  += ak * psi0_pow;

            psi0_pow *= psi_0;
        }

        gsl_vector_set(_L, 0, -sum_p.real());
        gsl_vector_set(_L, 1, -sum_in.real());
        gsl_vector_set(_L, 2, -sum_in.imag());

        const complex<double> entry = 1.0 - P0  * sum_0;
        gsl_vector_set(_L, 3, entry.real());

        // Invert M and solve for the constrained coefficients
        int signum = 0;
        gsl_permutation_init(_perm);
        gsl_linalg_LU_decomp(_M, _perm, &signum);
        gsl_linalg_LU_invert(_M, _perm, _inv_M);

        gsl_blas_dgemv(CblasNoTrans, 1.0, _inv_M, _L, 0.0, _constrained_coefficents);

        std::array<double, 4u> result;
        for (unsigned i = 0; i < 4; ++i)
        {
            result[i] = gsl_vector_get(_constrained_coefficents, i);
        }
        return result;
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_p_of_psi(const complex<double> & psi) const
    {
        // prepare expansion coefficients
        std::array<double, 13> a;
        const auto constrained_a = this->constrained_a_fp_I1();
        std::copy(constrained_a.cbegin(), constrained_a.cend(), a.begin()); // copy constrained coefficients
        std::copy(_a_fp_I1.cbegin(), _a_fp_I1.cend(), a.begin() + 4);       // copy unconstrained coefficients

        const auto series = this->series(psi, a);
        return this->_P(psi) * series;
    }


    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_p(const complex<double> & s) const
    {
        const complex<double> psi  = this->_s_to_psi_11(s);
        const complex<double> P    = this->_P(psi);

        // prepare expansion coefficients
        std::array<double, 13> a;
        const auto constrained_a = this->constrained_a_fp_I1();
        std::copy(constrained_a.cbegin(), constrained_a.cend(), a.begin()); // copy constrained coefficients
        std::copy(_a_fp_I1.cbegin(), _a_fp_I1.cend(), a.begin() + 4);       // copy unconstrained coefficients

        const auto series = this->series(psi, a);
        return P * series;
    }


    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_p(const double & s) const
    {
        static const double eps = 1.0e-14;
        return f_p(complex<double>(s, eps));
    }


    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_p_21(const complex<double> & s) const
    {
        const complex<double> psi  = this->_s_to_psi_21(s);
        const complex<double> P    = this->_P(psi);

        // prepare expansion coefficients
        std::array<double, 13> a;
        const auto constrained_a = this->constrained_a_fp_I1();
        std::copy(constrained_a.cbegin(), constrained_a.cend(), a.begin()); // copy constrained coefficients
        std::copy(_a_fp_I1.cbegin(), _a_fp_I1.cend(), a.begin() + 4);       // copy unconstrained coefficients

        return P * this->series(psi, a);
    }



    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_p_21(const double & s) const
    {
        static const double eps = 1.0e-14;
        return f_p_21(complex<double>(s, eps));
    }



    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::partial_wave(const complex<double> & s) const
    {
        const complex<double> s_p = this->_s_p();
        const double  eps = 1e-14;
        if (std::abs(s - s_p ) < eps)
            return 0.0;

        const complex<double> f_p_11 = this->f_p(s);
        const complex<double> f_p_21 = this->f_p_21(s);

        return std::sqrt(s)/ (2.0 * std::sqrt(s_p - s)) * (1.0 - f_p_11 / f_p_21);
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::partial_wave(const double & s) const
    {
        static const double eps = 1.0e-14;
        return partial_wave(complex<double>(s, eps));
    }


    std::array<complex<double>, 2>
    BHKMNR2026FormFactors<VacuumToPiPi>::scattering_lenght_parameters() const
    {
        const complex<double> s_p   = this->_s_p();
        const complex<double> f_p   = this->f_p(s_p);
        const complex<double> psi_p = this->_s_to_psi_11(s_p);
        const complex<double> P     = this->_P(psi_p);

        const auto derivs = this->_P_derivatives(psi_p);
        const std::complex<double> P1 = derivs.P1;
        const std::complex<double> P2 = derivs.P2;
        const std::complex<double> P3 = derivs.P3;

        // prepare expansion coefficients
        std::array<double, 13> a;
        const auto constrained_a = this->constrained_a_fp_I1();
        std::copy(constrained_a.cbegin(), constrained_a.cend(), a.begin()); // copy constrained coefficients
        std::copy(_a_fp_I1.cbegin(), _a_fp_I1.cend(), a.begin() + 4);       // copy unconstrained coefficients

        std::complex<double> psi_k   = 1.0;
        std::complex<double> psi_km1 = 0.0;
        std::complex<double> psi_km2 = 0.0;
        std::complex<double> psi_km3 = 0.0;

        std::complex<double> S  = 0.0; // S = \sum a_k psi^k
        std::complex<double> S1 = 0.0; // S1 = \sum a_k k psi^(k-1)
        std::complex<double> S2 = 0.0; // S2 = \sum a_k k (k-1) psi^(k-2)
        std::complex<double> S3 = 0.0; // S3 = \sum a_k k (k-1) (k-2) psi^(k-3)

        for (std::size_t k = 0; k < a.size(); ++k)
        {
            S  += a[k] * psi_k;
            S1 += a[k] * k * psi_km1;
            S2 += a[k] * k * (k - 1) * psi_km2;
            S3 += a[k] * k * (k - 1) * (k - 2) * psi_km3;

            psi_km3 = psi_km2;
            psi_km2 = psi_km1;
            psi_km1 = psi_k;
            psi_k  *= psi_p;
        }

        const complex<double> f2 = P * S2 + 2.0 * P1 * S1 + P2 * S;                 //d2f/dpsi2 at psi_p
        const complex<double> f3 = P * S3 + 3.0 * P1 * S2 + 3.0 * P2 * S1 + P3 * S; //d3f/dpsi3 at psi_p

        return {f2 / f_p, f3 / f_p};
    }


    //This function will be used to compute charge radius of pion
    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::dfdpsi_11(const complex<double> & s) const
    {
        const complex<double> psi  = this->_s_to_psi_11(s);
        const complex<double> P    = this->_P(psi);


        const std::complex<double> P1 = this->_dPdpsi(psi); //first derivative of P with respect to psi

        // prepare expansion coefficients
        std::array<double, 13> a;
        const auto constrained_a = this->constrained_a_fp_I1();
        std::copy(constrained_a.cbegin(), constrained_a.cend(), a.begin()); // copy constrained coefficients
        std::copy(_a_fp_I1.cbegin(), _a_fp_I1.cend(), a.begin() + 4);       // copy unconstrained coefficients

        std::complex<double> psi_k    = 1.0;
        std::complex<double> psi_km1  = 0.0;

        std::complex<double> S  = 0.0; // S = \sum a_k psi^k
        std::complex<double> S1 = 0.0; // S1 = \sum a_k k psi^(k-1)


        for (std::size_t k = 0; k < a.size(); ++k)
        {
            S  += a[k] * psi_k;
            S1 += a[k] * k * psi_km1;

            psi_km1 = psi_k;
            psi_k  *= psi;
        }

        return P1 * S + P * S1;
    }

    double BHKMNR2026FormFactors<VacuumToPiPi>::dispersive_integrand(const double & x) const
    {
        // change of variable s = s_p + x / (1 - x), then integral over s from s_p to infinity will become integral over x from 0 to 1
        const double Q2         = 1.0;
        const double chi        = 0.00683918; // GeV^-2, at Q^2 = 1 GeV^2 using [BL:1998A] Sec VI.A
        const double s_p        = real(this->_s_p());
        const double denom      = 48.0 * power_of<2>(M_PI) * chi;

        complex<double> f_p = this->f_p(s_p + x / (1.0 - x));

        return std::pow(x, 1.5) * std::norm(f_p) / std::sqrt(s_p * (1 - x) + x) / power_of<3>((s_p + Q2) * (1 - x) + x) / denom;
    }

    double BHKMNR2026FormFactors<VacuumToPiPi>::saturation() const
    {
        std::function<double (const double &)> f = [this](const double & x) -> double { return this->dispersive_integrand(x); };
        return integrate<1, 1>(f, 0, 1, cubature::Config().epsrel(1.0e-5));
    }

     complex<double> BHKMNR2026FormFactors<VacuumToPiPi>::residue(const unsigned & k) const
    {
        // prepare expansion coefficients
        std::array<double, 13> a;
        const auto constrained_a = this->constrained_a_fp_I1();
        std::copy(constrained_a.cbegin(), constrained_a.cend(), a.begin()); // copy constrained coefficients
        std::copy(_a_fp_I1.cbegin(), _a_fp_I1.cend(), a.begin() + 4);       // copy unconstrained coefficients

        const complex<double> psi_r     = this->_psi_r(this->_M_fp_I1[k](), this->_G_fp_I1[k]());
        const complex<double> series_r  = series(psi_r, a);
        const complex<double> P_r       = this->_P_residue(k);

        return P_r * series_r;
    }

    double BHKMNR2026FormFactors<VacuumToPiPi>::re_residue_rho() const
    {
        return std::real(this->residue(0u));
    }

    double BHKMNR2026FormFactors<VacuumToPiPi>::im_residue_rho() const
    {
        return std::imag(this->residue(0u));
    }

    //complex<double> BHKMNR2026FormFactors<VacuumToPiPi>::residue_rho_s() const
    //{
    //    const complex<double> s_rho = power_of<2>(complex<double>(this->_M_fp_I1[0u](), -this->_G_fp_I1[0u]()/2));
    //   return this->residue(0u) ;
    //}

    //double BHKMNR2026FormFactors<VacuumToPiPi>::re_residue_rho_s() const
    //{
    //    return std::real(this->residue_rho_s());
    //}

    //double BHKMNR2026FormFactors<VacuumToPiPi>::im_residue_rho_s() const
    //{
    //    return std::imag(this->residue_rho_s());
    //}



    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_0(const double & /*s*/) const
    {
        return 0.0; // vanishes in our approximation
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_0(const complex<double> & /*s*/) const
    {
        return 0.0; // vanishes in our approximation
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_t(const double & /*s*/) const
    {
        throw InternalError("Not implemented!");
        return 0.0;
    }

    complex<double>
    BHKMNR2026FormFactors<VacuumToPiPi>::f_t(const complex<double> & /*s*/) const
    {
        throw InternalError("Not implemented!");
        return 0.0;
    }

    const std::vector<OptionSpecification>
    BHKMNR2026FormFactors<VacuumToPiPi>::option_specifications
    {
        { "n-resonances"_ok, { "1"s, "2"s, "3"s }, "1"s }
    };

    std::vector<OptionSpecification>::const_iterator
    BHKMNR2026FormFactors<VacuumToPiPi>::begin_options()
    {
        return option_specifications.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BHKMNR2026FormFactors<VacuumToPiPi>::end_options()
    {
        return option_specifications.cend();
    }

    const std::set<ReferenceName>
    BHKMNR2026FormFactors<VacuumToPiPi>::references
    {
        "BHKMR:2025A"_rn
    };

}