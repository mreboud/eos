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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BHKMNR2026_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BHKMNR2026_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/reference-name.hh>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include <array>

namespace eos
{
    template <typename Process_> class BHKMNR2026FormFactors;

    template <> class BHKMNR2026FormFactors<VacuumToPiPi> :
        public FormFactors<VacuumToPP>
    {
        private:
            // parameters for form factor f_+ (I=1 projection)
            std::array<UsedParameter, 9u> _a_fp_I1; // unconstrained expansion coefficients
            std::array<UsedParameter, 3u> _M_fp_I1; // masses of the rho, rho', etc.
            std::array<UsedParameter, 3u> _G_fp_I1; // widths of the rho, rho', etc.

            RestrictedOption _n_resonances; // number of resonances

            UsedParameter _m_pi; // pion mass

            UsedParameter _s_0;  // free parameter
            UsedParameter _s_in; // inelastic threshold

            UsedParameter _hbar;

            // matrices and vectors needed for the linear system of equations to determine the constrained coefficients
            gsl_matrix * _M;
            gsl_matrix * _inv_M;
            gsl_vector * _L;
            gsl_permutation * _perm;
            gsl_vector * _constrained_coefficents;

            inline std::string _par_name(const std::string & ff, const std::string & isospin, const std::string & index) const
            {
                return "0->pipi::a_(" + ff + "," + isospin + ")^" + index + "@BHKMNR2026";
            }

            inline complex<double> _s_p() const
            {
                return complex<double>(4.0 * power_of<2>(_m_pi()), 0.0); // pair-production threshold s_plus
            }


            inline complex<double> _s_m() const
            {
                return complex<double>(0.0, 0.0); // start of the left-hand cut s_minus
            }


            inline complex<double> _phi_to_s(const complex<double> & phi, const complex<double> & s_in) const
            {
                const complex<double> s_p = _s_p();

                return (s_p * power_of<2>(1.0 + power_of<2>(phi)) - 4.0 * s_in * power_of<2>(phi)) / power_of<2>(1.0 - power_of<2>(phi));
            }


            inline complex<double> _s_to_phi_11(const complex<double> & s, const complex<double> & s_in) const
            {
                const complex<double> s_p = _s_p();
                const double eps          = 1e-14;

                if (std::abs(s - s_p) < eps)
                {
                    return complex<double>(0.0, 0.0);
                }

                return (std::sqrt(s_in - s) - std::sqrt(s_in - s_p)) / std::sqrt(s_p - s);
            }


            inline complex<double> _s_to_phi_21(const complex<double> & s, const complex<double> & s_in) const
            {
                const complex<double> s_p = _s_p();
                const double eps          = 1e-14;

                if (std::abs(s - s_p) < eps)
                {
                    return complex<double>(0.0, 0.0);
                }

                return (- std::sqrt(s_in - s) + std::sqrt(s_in - s_p)) / std::sqrt(s_p - s);
            }



            inline complex<double> _s_to_phi_22(const complex<double> & s, const complex<double> & s_in) const
            {
                const complex<double> s_p = _s_p();
                const double eps          = 1e-14;

                if (std::abs(s - s_p) < eps)
                {
                    return complex<double>(0.0, 0.0);
                }

                return (std::sqrt(s_in - s) + std::sqrt(s_in - s_p)) / std::sqrt(s_p - s);
            }



            inline complex<double> _s_to_phi_12(const complex<double> & s, const complex<double> & s_in) const
            {
                const complex<double> s_p = _s_p();
                const double eps          = 1e-14;

                if (std::abs(s - s_p) < eps)
                {
                    return complex<double>(0.0, 0.0);
                }

                return (- std::sqrt(s_in - s) - std::sqrt(s_in - s_p)) / std::sqrt(s_p - s);
            }


            // The name chi is chosen to be consistent with eq. (4.3) of arxiv:2510.25584
            inline complex<double> _chi(const complex<double> & x, const complex<double> & x_L, const complex<double> & x_0) const
            {
                const complex<double> A = (x * power_of<2>(x_L - 1.0) - x_L * power_of<2>(x - 1.0)) * power_of<2>(x_0 - 1.0);
                const complex<double> B = (x_0 * power_of<2>(x_L - 1.0) - x_L * power_of<2>(x_0 - 1.0)) * power_of<2>(x - 1.0);

                return (std::sqrt(A) - std::sqrt(B)) / (std::sqrt(A) + std::sqrt(B));
            }

            inline complex<double> _s_to_psi_11(const complex<double> & s) const
            {
                const complex<double> s_0   = _s_0();
                const complex<double> s_in  = _s_in();
                const complex<double> phi_L = _s_to_phi_21(_s_m(), s_in);

                return _chi(_s_to_phi_11(s, s_in), phi_L, _s_to_phi_11(s_0, s_in));
            }

            inline complex<double> _s_to_psi_21(const complex<double> & s) const
            {
                const complex<double> s_0   = _s_0();
                const complex<double> s_in  = _s_in();
                const complex<double> phi_L = _s_to_phi_21(_s_m(), s_in);

                return _chi(_s_to_phi_21(s, s_in), phi_L, _s_to_phi_11(s_0, s_in));
            }

            inline complex<double> _s_to_psi_22(const complex<double> & s) const
            {
                const complex<double> s_0   = _s_0();
                const complex<double> s_in  = _s_in();
                const complex<double> phi_L = _s_to_phi_21(_s_m(), s_in);

                return _chi(_s_to_phi_22(s, s_in), phi_L, _s_to_phi_11(s_0, s_in));
            }

            inline complex<double> _s_to_psi_12(const complex<double> & s) const
            {
                const complex<double> s_0   = _s_0();
                const complex<double> s_in  = _s_in();
                const complex<double> phi_L = _s_to_phi_21(_s_m(), s_in);

                return _chi(_s_to_phi_12(s, s_in), phi_L, _s_to_phi_11(s_0, s_in));
            }


            inline complex<double> _psi_r(const double & M, const double & Gamma) const
            {
                if (M * M < _s_in()) // the resonance is below the inelastic threshold, so we are on the 21 Riemann sheet
                {
                    return _s_to_psi_21(power_of<2>(complex<double>(M, -Gamma / 2.0)));
                }
                else // the resonance is above the inelastic threshold, so we are on the 22 Riemann sheet
                {
                    return _s_to_psi_22(power_of<2>(complex<double>(M, -Gamma / 2.0)));
                }
            }



            inline complex<double> _P(const complex<double> & psi) const
            {
                complex<double> psi_r;
                const std::size_t num_resonances = stoi(_n_resonances.value());
                complex<double> result           = power_of<2>(psi - 1.0);

                for (auto i = 0u; i < num_resonances; i++)
                {
                    psi_r = _psi_r(_M_fp_I1[i](), _G_fp_I1[i]());
                    result /= (psi - psi_r) * (psi - std::conj(psi_r));
                }

                return result;
            }


            inline complex<double> _P_residue(const unsigned & k) const
            {
                const std::size_t num_resonances = stoi(_n_resonances.value());

                if (k > num_resonances)
                    throw InternalError("The residue index must be smaller than the number of used resonances.");

                complex<double> psi_residue = _psi_r(_M_fp_I1[k](), _G_fp_I1[k]());
                complex<double> result      = power_of<2>(psi_residue - 1.0) / (psi_residue - std::conj(psi_residue));

                for (auto i = 0u; i < num_resonances; i++)
                {
                    if (i != k)
                    {
                        complex<double> psi_r = _psi_r(_M_fp_I1[i](), _G_fp_I1[i]());
                        result /= (psi_residue - psi_r) * (psi_residue - std::conj(psi_r));
                    }
                }

                return result;
            }



            inline complex<double> _dPdpsi(const complex<double> & psi) const
            {
                const std::size_t num_resonances = stoi(_n_resonances.value());
                const complex<double> P_val      = _P(psi);
                complex<double> sum              = complex<double>(0.0, 0.0);

                for (auto i = 0u; i < num_resonances; ++i)
                {
                    const complex<double> psi_r = _psi_r(_M_fp_I1[i](), _G_fp_I1[i]());
                    const complex<double> denom = (psi - psi_r) * (psi - std::conj(psi_r));

                    sum += (2.0 * psi - psi_r - std::conj(psi_r)) / denom;
                }

                return P_val * (2.0 / (psi - 1.0) - sum);
            }



            inline complex<double> _dfdpsi_terms(const unsigned k, const complex<double> & psi) const
            {
                const complex<double> dP_val = _dPdpsi(psi);

                switch (k)
                {
                    case 0:
                        return dP_val;
                    case 1:
                        return dP_val * psi + _P(psi);
                    default:
                        complex<double> psi_km1 = std::pow(psi, k - 1);
                        complex<double> psi_k   = psi_km1 * psi;
                        return dP_val * psi_k + static_cast<double>(k) * _P(psi) * psi_km1;
                }
            }


            //This function will be used to find scattering lenght
            struct PDerivatives
            {
                complex<double> P1;  // first derivative respect to psi
                complex<double> P2;  // second derivative respect to psi
                complex<double> P3;  // third derivative respect to psi
            };

            inline PDerivatives _P_derivatives(const complex<double>& psi) const
            {
                const std::size_t num_resonances = stoi(_n_resonances.value());

                const complex<double> P_val = _P(psi);

                complex<double> L  =  2.0 / (psi - 1.0);
                complex<double> L1 = -2.0 / power_of<2>(psi - 1.0);
                complex<double> L2 =  4.0 / power_of<3>(psi - 1.0);

                for (auto i = 0u; i < num_resonances; ++i)
                {
                    const complex<double> psi_r = _psi_r(_M_fp_I1[i](), _G_fp_I1[i]());

                    const complex<double> a = psi - psi_r;
                    const complex<double> b = psi - std::conj(psi_r);

                    const complex<double> D  = a * b;
                    const complex<double> D2 = D * D;
                    const complex<double> D3 = D2 * D;

                    const complex<double> N  = 2.0 * psi - psi_r - std::conj(psi_r);
                    const complex<double> N2 = N * N;
                    const complex<double> N3 = N2 * N;

                    L -= N / D;

                    L1 -= 2.0 / D;
                    L1 += N2 / D2;

                    L2 += 6.0 * N / D2;
                    L2 -= 2.0 * N3 / D3;
                }

                PDerivatives out;

                const complex<double> term2 = L*L + L1;
                const complex<double> term3 = L*L*L + 3.0*L*L1 + L2;

                out.P1 = P_val * L;
                out.P2 = P_val * term2;
                out.P3 = P_val * term3;

                return out;
            }


        public:
            BHKMNR2026FormFactors(const Parameters & p, const Options & o);
            ~BHKMNR2026FormFactors();

            static FormFactors<VacuumToPP> * make(const Parameters & p, const Options & o);

            /* auxiliary functions */
            std::array<double, 4u> constrained_a_fp_I1() const;
            complex<double> psi(const complex<double> & s) const;
            complex<double> P(const complex<double> & psi) const;
            complex<double> dPdpsi(const complex<double> & psi) const;
            complex<double> dfdpsi_terms(const unsigned k, const complex<double> & psi) const;
            complex<double> series(const complex<double> & psi, const std::array<double, 13> & a) const;
            complex<double> f_p_of_psi(const complex<double> & psi) const;
            double abs2_f_p_of_psi(const double & re_psi, const double & im_psi) const
            {
                return std::norm(f_p_of_psi(complex<double>(re_psi, im_psi)));
            }
            double arg_f_p_of_psi(const double & re_psi, const double & im_psi) const
            {
                return std::arg(f_p_of_psi(complex<double>(re_psi, im_psi)));
            }

            /* form factors on the real axis */
            virtual complex<double> f_p(const double & s) const override;
            virtual complex<double> f_0(const double & s) const override;
            virtual complex<double> f_t(const double & s) const override;

            /* form factor in the complex s plane */
            virtual complex<double> f_p(const complex<double> & s) const override;
            virtual complex<double> f_0(const complex<double> & s) const override;
            virtual complex<double> f_t(const complex<double> & s) const override;

            /* form factors on the 21 Rieman sheet */
            complex<double> f_p_21(const double & s) const;
            complex<double> f_p_21(const complex<double> & s) const;

            /* Isospin 1, P-wave partial wave */
            complex<double> partial_wave(const double & s) const;
            complex<double> partial_wave(const complex<double> & s) const;
            double re_partial_wave(const double & s) const
            {
                return std::real(partial_wave(s));
            }
            double im_partial_wave(const double & s) const
            {
                return std::imag(partial_wave(s));
            }
            std::array<complex<double>, 2> scattering_lenght_parameters() const;
            double d2fdpsi2_over_f() const
            {
                return std::real(scattering_lenght_parameters()[0]);
            }
            double d3fdpsi3_over_f() const
            {
                return std::real(scattering_lenght_parameters()[1]);
            }

            complex<double> dfdpsi_11(const complex<double> & s) const;
            double dfdpsi_11_at_0() const
            {
                return std::real(dfdpsi_11(0.0));
            }
            double dispersive_integrand(const double & s) const;
            double saturation() const;

            //residue functions
            complex<double> residue(const unsigned & k) const;
            double re_residue_rho() const;
            double im_residue_rho() const;
            //complex<double> residue_rho_s() const
            //double re_residue_rho_s() const
            //double im_residue_rho_s() const


            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
            static const std::vector<OptionSpecification> option_specifications;

            static const std::set<ReferenceName> references;
    };

    extern template class BHKMNR2026FormFactors<VacuumToPiPi>;
}

#endif
