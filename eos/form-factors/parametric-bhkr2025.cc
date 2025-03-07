/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Méril Reboud
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

#include <eos/form-factors/parametric-bhkr2025.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/stringify.hh>
#include <eos/maths/integrate.hh>

#include <functional>
#include <numeric>

namespace eos
{
    /* Vacuum -> pi pi */
    BHKR2025FormFactors<VacuumToPiPi>::BHKR2025FormFactors(const Parameters & p, const Options & o) :
        n_resonances_1m(o, option_specifications, "n-resonances-1m"_ok),
        _b_fp_I1{{
            UsedParameter(p[_par_name("+", "1", "5")], *this),
            UsedParameter(p[_par_name("+", "1", "6")], *this),
            UsedParameter(p[_par_name("+", "1", "7")], *this),
            UsedParameter(p[_par_name("+", "1", "8")], *this),
            UsedParameter(p[_par_name("+", "1", "9")], *this)
        }},
        _M_fp_I1{{
            UsedParameter(p["0->pipi::M_(+,0)@BHKR2025"], *this),
            UsedParameter(p["0->pipi::M_(+,1)@BHKR2025"], *this),
            UsedParameter(p["0->pipi::M_(+,2)@BHKR2025"], *this)
        }},
        _G_fp_I1{{
            UsedParameter(p["0->pipi::Gamma_(+,0)@BHKR2025"], *this),
            UsedParameter(p["0->pipi::Gamma_(+,1)@BHKR2025"], *this),
            UsedParameter(p["0->pipi::Gamma_(+,2)@BHKR2025"], *this),
        }},
        _m_pi(p["mass::pi^+"], *this),
        _t_in(p["0->pipi::s_th@BHKR2025"], *this),
        _hbar(p["QM::hbar"], *this),
        _tmp(gsl_matrix_alloc(5, 5)),
        _tmp_inverse(gsl_matrix_alloc(5, 5)),
        _perm(gsl_permutation_alloc(5)),
        _tmp_vec(gsl_vector_alloc(5)),
        _b_vector(gsl_vector_alloc(5))
    {
    }

    BHKR2025FormFactors<VacuumToPiPi>::~BHKR2025FormFactors() = default;

    FormFactors<VacuumToPP> *
    BHKR2025FormFactors<VacuumToPiPi>::make(const Parameters & p, const Options & o)
    {
        return new BHKR2025FormFactors<VacuumToPiPi>(p, o);
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::zeta(const double & q2) const
    {
        return _zeta(complex<double>(q2, -std::numeric_limits<double>::min()));
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::zeta(const complex<double> & q2) const
    {
        return _zeta(q2);
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::q2(const complex<double> & zeta) const
    {
        return _q2(zeta);
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::dzetadq2(const complex<double> & q2) const
    {
        const double t_p  = this->_t_p();
        const double t_in = this->_t_in();

        return zeta(q2) * (t_in - t_p) / (2.0 * (q2 - t_p) * sqrt((q2 - t_in) * (t_p - t_in)));
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::dzetadq2_21(const complex<double> & q2) const
    {
        // dzeta/dq2 on the 21 Riemann sheet
        return -dzetadq2(q2);
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::w(const complex<double> & zeta) const
    {
        return power_of<2>(zeta * zeta - 1.0);
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::xi(const complex<double> & zeta) const
    {
        // Position of the branch point on the second Riemann sheet
        const complex<double> zeta_L = -this->zeta(0);
        if (real(zeta) < 1.0)
        {
            return (sqrt(zeta * power_of<2>(zeta_L - 1.0) - zeta_L * power_of<2>(zeta - 1.0)) - sqrt( - zeta_L * power_of<2>(zeta - 1.0)))
                / (sqrt(zeta * power_of<2>(zeta_L - 1.0) - zeta_L * power_of<2>(zeta - 1.0)) + sqrt( - zeta_L * power_of<2>(zeta - 1.0)));
        }
        else
        {
            return (sqrt(zeta * power_of<2>(zeta_L - 1.0) - zeta_L * power_of<2>(zeta - 1.0)) + sqrt( - zeta_L * power_of<2>(zeta - 1.0)))
                / (sqrt(zeta * power_of<2>(zeta_L - 1.0) - zeta_L * power_of<2>(zeta - 1.0)) - sqrt( - zeta_L * power_of<2>(zeta - 1.0)));

        }
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::resonance_product_p(const complex<double> & zeta) const
    {
        const std::size_t num_resonances = stoi(n_resonances_1m.value());
        complex<double> result = 1.0;

        for (auto i = 0u; i < num_resonances; i++)
        {
            result *= this->_resonance_poles(_M_fp_I1[i], _G_fp_I1[i], zeta);
        }

        return result;
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::w_prime(const complex<double> & zeta) const
    {
        return 4.0 * zeta * (zeta * zeta - 1.0);
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::resonance_product_p_prime(const complex<double> & zeta) const
    {
        const std::size_t num_resonances = stoi(n_resonances_1m.value());
        complex<double> tmp = 0.0;

        for (auto i = 0u; i < num_resonances; i++)
        {
            tmp += this->_resonance_poles_prime(_M_fp_I1[i], _G_fp_I1[i], zeta) / this->_resonance_poles(_M_fp_I1[i], _G_fp_I1[i], zeta);
        }

        return resonance_product_p(zeta) * tmp;
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::series_m(const complex<double> & zeta, const std::array<double, 10u> & b) const
    {
        std::array<complex<double>, 10> xi_values;
        complex<double> current_xi_values = 1.0;
        complex<double> xi = this->xi(zeta);

        for (complex<double> & xiv: xi_values)
        {
            xiv = current_xi_values;
            current_xi_values *= xi;
        }

        return std::inner_product(b.cbegin(), b.cend(), xi_values.cbegin(), complex<double>(0.0, 0.0));
    }

    std::array<double, 5u>
    BHKR2025FormFactors<VacuumToPiPi>::_fixed_b_fp_I1() const
    {
        const complex<double> zeta_0 = zeta(0.0);
        const double xi_0 = real(xi(zeta(0.0))); // Real by construction
        const double w_p_zeta_0 = real(w(zeta_0) * resonance_product_p(zeta_0)); // Real by construction
        gsl_matrix_set(_tmp, 0, 0, w_p_zeta_0);
        gsl_matrix_set(_tmp, 0, 1, w_p_zeta_0 * xi_0);
        gsl_matrix_set(_tmp, 0, 2, w_p_zeta_0 * xi_0 * xi_0);
        gsl_matrix_set(_tmp, 0, 3, w_p_zeta_0 * power_of<3>(xi_0));
        gsl_matrix_set(_tmp, 0, 4, w_p_zeta_0 * power_of<4>(xi_0));

        const complex<double> w_0 = w(complex<double>(0.0, 0.0)),
                              p_0 = resonance_product_p(complex<double>(0.0, 0.0)),
                              w_prime_0 = w_prime(complex<double>(0.0, 0.0)),
                              p_prime_0 = resonance_product_p_prime(complex<double>(0.0, 0.0)),
                              xi_prime_0 = power_of<2>(zeta_0 + 1.0) / (4.0 * zeta_0);
        gsl_matrix_set(_tmp, 1, 0, real(p_0 * w_prime_0 + w_0 * p_prime_0));
        gsl_matrix_set(_tmp, 1, 1, real(p_0 * w_0 * xi_prime_0));
        gsl_matrix_set(_tmp, 1, 2, 0.0);
        gsl_matrix_set(_tmp, 1, 3, 0.0);
        gsl_matrix_set(_tmp, 1, 4, 0.0);

        const complex<double> w_I = w(complex<double>(0.0, 1.0)),
                              p_I = resonance_product_p(complex<double>(0.0, 1.0)),
                              xi_I = xi(complex<double>(0.0, 1.0)),
                              w_prime_I = w_prime(complex<double>(0.0, 1.0)),
                              p_prime_I = resonance_product_p_prime(complex<double>(0.0, 1.0)),
                              xi_prime_I = complex<double>(0.0, -1.0) * xi_I * std::sqrt(2.0 * zeta_0 / (zeta_0 * zeta_0 + 1.0));
        gsl_matrix_set(_tmp, 2, 0, 2.0 * real(p_I * w_prime_I + w_I * p_prime_I));
        gsl_matrix_set(_tmp, 2, 1, 2.0 * real(p_I * w_I * xi_prime_I + xi_I * (p_I * w_prime_I + w_I * p_prime_I)));
        gsl_matrix_set(_tmp, 2, 2, 2.0 * real(2.0 * p_I * w_I * xi_I * xi_prime_I + xi_I * xi_I * (p_I * w_prime_I + w_I * p_prime_I)));
        gsl_matrix_set(_tmp, 2, 3, 2.0 * real(3.0 * p_I * w_I * xi_I * xi_I * xi_prime_I + power_of<3>(xi_I) * (p_I * w_prime_I + w_I * p_prime_I)));
        gsl_matrix_set(_tmp, 2, 4, 2.0 * real(4.0 * p_I * w_I * power_of<3>(xi_I) * xi_prime_I + power_of<4>(xi_I) * (p_I * w_prime_I + w_I * p_prime_I)));

        gsl_matrix_set(_tmp, 3, 0, 2.0 * imag(p_I * w_prime_I + w_I * p_prime_I));
        gsl_matrix_set(_tmp, 3, 1, 2.0 * imag(p_I * w_I * xi_prime_I + xi_I * (p_I * w_prime_I + w_I * p_prime_I)));
        gsl_matrix_set(_tmp, 3, 2, 2.0 * imag(2.0 * p_I * w_I * xi_I * xi_prime_I + xi_I * xi_I * (p_I * w_prime_I + w_I * p_prime_I)));
        gsl_matrix_set(_tmp, 3, 3, 2.0 * imag(3.0 * p_I * w_I * xi_I * xi_I * xi_prime_I + power_of<3>(xi_I) * (p_I * w_prime_I + w_I * p_prime_I)));
        gsl_matrix_set(_tmp, 3, 4, 2.0 * imag(4.0 * p_I * w_I * power_of<3>(xi_I) * xi_prime_I + power_of<4>(xi_I) * (p_I * w_prime_I + w_I * p_prime_I)));

        gsl_matrix_set(_tmp, 4, 0,  0.0);
        gsl_matrix_set(_tmp, 4, 1,  1.0);
        gsl_matrix_set(_tmp, 4, 2, -2.0);
        gsl_matrix_set(_tmp, 4, 3,  3.0);
        gsl_matrix_set(_tmp, 4, 4, -4.0);

        int signum = 0;
        gsl_permutation_init(_perm);
        gsl_linalg_LU_decomp(_tmp, _perm, &signum);
        gsl_linalg_LU_invert(_tmp, _perm, _tmp_inverse);

        std::array<double, 10> b = { 0.0 };
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin() + 5);

        // compute the sums {n * b[n] * xi_I^(n-1)} and {n * b[n] * (-1)^(n-1)}
        complex<double> xi_I_sum = 0.0;
        complex<double> current_xi_value = 1.0;
        double xi_m1_sum = 0.0;
        int current_m1_value = -1.0;
        for (unsigned i = 1u; i < b.size(); i++)
        {
            xi_I_sum += i * b[i] * current_xi_value;
            current_xi_value *= xi_I;
            xi_m1_sum += i * b[i] * current_m1_value;
            current_m1_value *= -1.0;
        }

        gsl_vector_set(_tmp_vec, 0, 1.0 - w_p_zeta_0 * real(series_m(zeta_0, b)));
        gsl_vector_set(_tmp_vec, 1, 0.0);
        gsl_vector_set(_tmp_vec, 2, -2.0 * real((p_I * w_prime_I + w_I * p_prime_I) * series_m(complex<double>(0.0, 1.0), b) + xi_prime_I * w_I * p_I * xi_I_sum));
        gsl_vector_set(_tmp_vec, 3, -2.0 * imag((p_I * w_prime_I + w_I * p_prime_I) * series_m(complex<double>(0.0, 1.0), b) + xi_prime_I * w_I * p_I * xi_I_sum));
        gsl_vector_set(_tmp_vec, 4, xi_m1_sum);

        gsl_blas_dgemv(CblasNoTrans, 1.0, _tmp_inverse, _tmp_vec, 0.0, _b_vector);

        return
        {
            gsl_vector_get(_b_vector, 0),
            gsl_vector_get(_b_vector, 1),
            gsl_vector_get(_b_vector, 2),
            gsl_vector_get(_b_vector, 3),
            gsl_vector_get(_b_vector, 4)
        };
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::f_p_zeta(const complex<double> & zeta) const
    {
        // prepare expansion coefficients
        std::array<double, 10> b;
        // Fix b[0] to b[5] by enforcing F(q2=0) = 1 and F'(q2=t+) = F'(q2=t_th) = F'(q2=t_yh) = F'_21(q2=0) = 0
        const std::array<double, 5u> b_vector = _fixed_b_fp_I1();
        std::copy(b_vector.cbegin(), b_vector.cend(), b.begin());
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin() + 5);

        return w(zeta) * resonance_product_p(zeta) * this->series_m(zeta, b);
    }

    double
    BHKR2025FormFactors<VacuumToPiPi>::re_f_p_zeta(const double & re_zeta, const double & im_zeta) const
    {
        return real(this->f_p_zeta(complex<double>(re_zeta, im_zeta)));
    }

    double
    BHKR2025FormFactors<VacuumToPiPi>::im_f_p_zeta(const double & re_zeta, const double & im_zeta) const
    {
        return imag(this->f_p_zeta(complex<double>(re_zeta, im_zeta)));
    }

    double
    BHKR2025FormFactors<VacuumToPiPi>::abs2_f_p_zeta(const double & re_zeta, const double & im_zeta) const
    {
        return norm(this->f_p_zeta(complex<double>(re_zeta, im_zeta)));
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::f_p(const double & q2) const
    {
        return f_p(complex<double>(q2, std::numeric_limits<double>::min()));
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::f_p(const complex<double> & q2) const
    {
        const auto zeta = this->zeta(q2);

        return f_p_zeta(zeta);
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::f_t(const double & /*q2*/) const
    {
        throw InternalError("Not implemented!");
        return 0.0;
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::f_t(const complex<double> & /*q2*/) const
    {
        throw InternalError("Not implemented!");
        return 0.0;
    }

    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::f_0(const double & /*q2*/) const
    {
        return 0.0; // vanishes in our approximation
    }


    complex<double>
    BHKR2025FormFactors<VacuumToPiPi>::f_0(const complex<double> & /*q2*/) const
    {
        return 0.0; // vanishes in our approximation
    }

    double
    BHKR2025FormFactors<VacuumToPiPi>::b_0() const
    {
        return _fixed_b_fp_I1()[0];
    }

    double
    BHKR2025FormFactors<VacuumToPiPi>::b_1() const
    {
        return _fixed_b_fp_I1()[1];
    }

    double
    BHKR2025FormFactors<VacuumToPiPi>::b_2() const
    {
        return _fixed_b_fp_I1()[2];
    }

    double
    BHKR2025FormFactors<VacuumToPiPi>::b_3() const
    {
        return _fixed_b_fp_I1()[3];
    }

    double
    BHKR2025FormFactors<VacuumToPiPi>::b_4() const
    {
        return _fixed_b_fp_I1()[4];
    }

    // double
    // BHKR2025FormFactors<VacuumToPiPi>::dFdq2_q2eq0() const
    // {
    //     const double z0 = std::real(this->z(0.0));
    //     const double chi = 0.00683918; // GeV^-2, at Q^2 = 1 GeV^2 using [BL:1998A] Sec VI.A

    //     const double phitilde_z0      = std::real(this->phitilde_p(z0, chi));
    //     const double phitildeprime_z0 = std::real(this->phitildeprime_p(z0, chi));

    //     // Super-threshold pole location
    //     const auto zr =  this->_zr(this->_M_fp_I1(), this->_G_fp_I1());

    //     // prepare expansion coefficients
    //     std::array<double, 10> b;
    //     std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);
    //     // Fix b[0] and b[1] to enforce F(q2=0) = 1 and F'(q2=t+) = 0
    //     b[0] = _b0_fp_I1(chi, zr);
    //     b[1] = _b1_fp_I1(chi, zr);

    //     const double sum1 = std::real(series_m(z0, b));
    //     double sum2 = 0.0;
    //     for (auto i = 0u; i < b.size(); i++)
    //     {
    //         sum2 += b[i] * i * std::pow(z0, i-1);
    //     };

    //     const double xprime_z0 = ( 2.0 * (std::real(zr) - z0) * phitilde_z0 - std::norm(z0 - zr) * phitildeprime_z0 ) / power_of<2>(std::norm(z0 - zr) * phitilde_z0);

    //     const double dFdz_z0 = sum1 * xprime_z0 + sum2 / (phitilde_z0 * std::norm(z0 - zr));
    //     return dFdz_z0 * std::real(this->dzdq2(0.0));
    // }

    // double
    // BHKR2025FormFactors<VacuumToPiPi>::r_pi_squared() const
    // {
    //     return 6.0 * this->dFdq2_q2eq0() * power_of<2>(this->_hbarc());
    // }

    // complex<double> BHKR2025FormFactors<VacuumToPiPi>::residue_rho() const
    // {
    //     // Super-threshold pole location
    //     const auto zr       = this->_zr(this->_M_fp_I1(), this->_G_fp_I1());
    //     const auto chi      = 0.00683918; // GeV^-2, at Q^2 = 1 GeV^2 using [BL:1998A] Sec VI.A
    //     const auto phitilde = this->phitilde_p(zr, chi);

    //     // prepare expansion coefficients
    //     std::array<double, 10> b;
    //     std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);
    //     // Fix b[0] and b[1] to enforce F(q2=0) = 1 and F'(q2=t+) = 0
    //     b[0] = _b0_fp_I1(chi, zr);
    //     b[1] = _b1_fp_I1(chi, zr);

    //     const auto series = this->series_m(zr, b);

    //     return series / (zr - std::conj(zr)) /  phitilde;
    // }

    // double BHKR2025FormFactors<VacuumToPiPi>::re_residue_rho() const
    // {
    //     return std::real(this->residue_rho());
    // }

    // double BHKR2025FormFactors<VacuumToPiPi>::im_residue_rho() const
    // {
    //     return std::imag(this->residue_rho());
    // }

    // complex<double> BHKR2025FormFactors<VacuumToPiPi>::residue_rho_q2() const
    // {
    //     const auto s_rho = power_of<2>(complex<double>(this->_M_fp_I1(), -this->_G_fp_I1()/2));
    //     return this->residue_rho() / this->dzdq2_II(s_rho);
    // }

    // double BHKR2025FormFactors<VacuumToPiPi>::re_residue_rho_q2() const
    // {
    //     return std::real(this->residue_rho_q2());
    // }

    // double BHKR2025FormFactors<VacuumToPiPi>::im_residue_rho_q2() const
    // {
    //     return std::imag(this->residue_rho_q2());
    // }

    double BHKR2025FormFactors<VacuumToPiPi>::dispersive_integrand(const complex<double> & zeta) const
    {
        const double t_p     = this->_t_p();
        const double Q2      = 1.0;
        const double chi     = 0.00683918; // GeV^-2, at Q^2 = 1 GeV^2 using [BL:1998A] Sec VI.A
        const double t       = real(q2(zeta)); // Real by construction

        const double prefactor = 1.0 / (48.0 * M_PI * M_PI * chi) / power_of<3>(t + Q2) * pow(1.0 - t_p / t, 1.5) / abs(dzetadq2(t));

        return prefactor * abs2_f_p_zeta(real(zeta), imag(zeta));
    }

    double BHKR2025FormFactors<VacuumToPiPi>::saturation() const
    {
        // The integral is first performed over the positive imaginary zeta axis
        std::function<double (const double &)> f_first = [this](const double & x) -> double
        {
            return this->dispersive_integrand(complex<double>(0.0, x));
        };
        // The integral is then performed over the [0, pi/2] zeta quadrant
        std::function<double (const double &)> f_second = [this](const double & alpha) -> double
        {
            return this->dispersive_integrand(complex<double>(cos(alpha), sin(alpha)));
        };

        return integrate<GSL::QAGS>(f_first, 0.0, 1.0) + integrate<GSL::QAGS>(f_second, 0.0, M_PI_2);
    }

    const std::set<ReferenceName>
    BHKR2025FormFactors<VacuumToPiPi>::references
    {
    };

    const std::vector<OptionSpecification>
    BHKR2025FormFactors<VacuumToPiPi>::option_specifications
    {
        { "n-resonances-1m"_ok, { "1", "2", "3" }, "1" },
    };

    std::vector<OptionSpecification>::const_iterator
    BHKR2025FormFactors<VacuumToPiPi>::begin_options()
    {
        return option_specifications.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BHKR2025FormFactors<VacuumToPiPi>::end_options()
    {
        return option_specifications.cend();
    }

    template class BHKR2025FormFactors<VacuumToPiPi>;
}
