/* vim: set sw=4 sts=4 et tw=120 foldmethod=syntax : */

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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BHKR2025_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BHKR2025_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/reference-name.hh>

#include <array>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

namespace eos
{
    template <typename Process_> class BHKR2025FormFactors;

    template <> class BHKR2025FormFactors<VacuumToPiPi> :
        public FormFactors<VacuumToPP>
    {
        private:
            // number of resonances in the 1- channel
            SpecifiedOption n_resonances_1m;

            // parameters for form factor f_+ (I=1 projection)
            std::array<UsedParameter, 5u> _b_fp_I1;
            std::array<UsedParameter, 3u> _M_fp_I1;
            std::array<UsedParameter, 3u> _G_fp_I1;

            // hadron masses
            UsedParameter _m_pi;

            // squared energy of the inelastic threshold
            UsedParameter _t_in;

            UsedParameter _hbar;

            // A temporary matrix and a permutation are needed to evaluate the first four coefficients
            // of the series expansion of the form factor f_+
            gsl_matrix * _tmp;
            gsl_matrix * _tmp_inverse;
            gsl_permutation * _perm;
            gsl_vector * _tmp_vec;
            gsl_vector * _b_vector;

            inline std::string _par_name(const std::string & ff, const std::string & isospin, const std::string & index) const
            {
                return "0->pipi::b_(" + ff + "," + isospin + ")^" + index + "@BHKR2025";
            }

            inline double _hbarc() const
            {
                return _hbar() * 299792458 * 1e15; // GeV fm
            }

            inline double _t_p() const
            {
                return 4.0 * _m_pi() * _m_pi();
            }

            inline complex<double> _zeta(const complex<double> & q2) const
            {
                const auto t_p  = _t_p();
                const auto t_in = _t_in();

                return sqrt((q2 + t_p + 2.0 * sqrt((q2 - t_in) * (t_p - t_in)) - 2.0 * t_in) / (q2 - t_p));
            }

            inline complex<double> _q2(const complex<double> & zeta) const
            {
                const auto t_p  = _t_p();
                const auto t_in = _t_in();

                return (-4.0 * t_in * zeta * zeta + t_p * power_of<2>(1.0 + zeta * zeta)) / power_of<2>(zeta * zeta - 1.0);
            }

            inline complex<double> _resonance_poles(const double & M, const double & Gamma, const complex<double> & zeta) const
            {
                if (M * M < _t_p())
                {
                    // The pole is on the 11 Riemann sheet
                    complex<double> zeta_r = _zeta(power_of<2>(complex<double>(M, -0.5 * Gamma)));
                    return 1.0 / (zeta - zeta_r) / (zeta - std::conj(zeta_r));
                }
                else if (M * M < _t_in())
                {
                    // The pole is on the 21 Riemann sheet
                    complex<double> zeta_r = -_zeta(power_of<2>(complex<double>(M, -0.5 * Gamma)));
                    return 1.0 / (zeta - zeta_r) / (zeta - std::conj(zeta_r));
                }
                else
                {
                    // The poles are on the 22 and 21 Riemann sheets, we neglect the one on the 12 Riemann sheet
                    complex<double> zeta_r = _zeta(power_of<2>(complex<double>(M, -0.5 * Gamma)));
                    return 1.0 / (zeta - 1.0 / zeta_r) / (zeta - 1.0 / std::conj(zeta_r)) / (zeta + zeta_r) / (zeta + std::conj(zeta_r));
                }
            }

            inline complex<double> _resonance_poles_prime(const double & M, const double & Gamma, const complex<double> & zeta) const
            {
                if (M * M < _t_p())
                {
                    // The pole is on the 11 Riemann sheet
                    complex<double> zeta_r = _zeta(power_of<2>(complex<double>(M, -0.5 * Gamma)));
                    return -2.0 * (zeta - std::real(zeta_r)) / power_of<2>((zeta - zeta_r) * (zeta - std::conj(zeta_r)));
                }
                else if (M * M < _t_in())
                {
                    // The pole is on the 21 Riemann sheet
                    complex<double> zeta_r = -_zeta(power_of<2>(complex<double>(M, -0.5 * Gamma)));
                    return -2.0 * (zeta - std::real(zeta_r)) / power_of<2>((zeta - zeta_r) * (zeta - std::conj(zeta_r)));
                }
                else
                {
                    // The poles are on the 22 and 21 Riemann sheets, we neglect the one on the 12 Riemann sheet
                    complex<double> zeta_r = _zeta(power_of<2>(complex<double>(M, -0.5 * Gamma)));
                    return (-2.0 * (zeta - std::real(1.0 / zeta_r)) / (zeta - 1.0 / zeta_r) / (zeta - 1.0 / std::conj(zeta_r)) \
                        -2.0 * (zeta + std::real(zeta_r)) / (zeta + zeta_r) / (zeta + std::conj(zeta_r))) \
                        / (zeta - 1.0 / zeta_r) / (zeta - 1.0 / std::conj(zeta_r)) / (zeta + zeta_r) / (zeta + std::conj(zeta_r));
                }
            }

            std::array<double, 5u> _fixed_b_fp_I1() const;

        public:
            BHKR2025FormFactors(const Parameters & p, const Options & o);
            ~BHKR2025FormFactors();

            static FormFactors<VacuumToPP> * make(const Parameters & p, const Options & o);

            /* auxiliary functions */
            complex<double> zeta(const double & q2) const;
            complex<double> zeta(const complex<double> & q2) const;
            complex<double> q2(const complex<double> & zeta) const;
            complex<double> dzetadq2(const complex<double> & q2) const;
            complex<double> dzetadq2_21(const complex<double> & q2) const;
            complex<double> w(const complex<double> & zeta) const;
            complex<double> xi(const complex<double> & zeta) const;
            complex<double> resonance_product_p(const complex<double> & zeta) const;
            complex<double> series_m(const complex<double> & zeta, const std::array<double, 10u> & c) const;
            // Derivatives with respect to zeta
            complex<double> w_prime(const complex<double> & zeta) const;
            complex<double> resonance_product_p_prime(const complex<double> & zeta) const;

            /* form factors as a function of zeta */
            complex<double> f_p_zeta(const complex<double> & q2) const;
            double re_f_p_zeta(const double & q2_real, const double & q2_imag) const;
            double im_f_p_zeta(const double & q2_real, const double & q2_imag) const;
            double abs2_f_p_zeta(const double & q2_real, const double & q2_imag) const;

            /* form factors on the real axis */
            virtual complex<double> f_p(const double & q2) const override;
            virtual complex<double> f_t(const double & q2) const override;
            virtual complex<double> f_0(const double & q2) const override;

            /* form factor in the complex q2 plane */
            virtual complex<double> f_p(const complex<double> & q2) const override;
            virtual complex<double> f_t(const complex<double> & q2) const override;
            virtual complex<double> f_0(const complex<double> & q2) const override;

            /* auxiliary observables */
            double r_pi_squared() const; // squared charge radius of the pion

            /* auxiliary pseudo observables */
            double b_0() const; // value of the series coefficient b_0
            double b_1() const; // value of the series coefficient b_1
            double b_2() const; // value of the series coefficient b_2
            double b_3() const; // value of the series coefficient b_3
            double b_4() const; // value of the series coefficient b_4
            // complex<double> residue_rho() const; // complex-valued residue of the form factor on the rho pole
            // double re_residue_rho() const;
            // double im_residue_rho() const;
            // complex<double> residue_rho_q2() const; // complex-valued residue in q2 of the form factor on the rho pole
            // double re_residue_rho_q2() const;
            // double im_residue_rho_q2() const;

            /* saturation of the dispersive bound */
            double dispersive_integrand(const complex<double> & zeta) const;
            double saturation() const;

            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
            static const std::vector<OptionSpecification> option_specifications;

            static const std::set<ReferenceName> references;
    };

    extern template class BHKR2025FormFactors<VacuumToPiPi>;
}

#endif
