/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2020 Danny van Dyk
 * Copyright (c) 2020 Nico Gubernari
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

#include <eos/rare-b-decays/nonlocal-formfactors.hh>
#include <eos/form-factors/b-lcdas.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/integrate-impl.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>

#include <cmath>
#include <map>

#include <gsl/gsl_sf.h>

namespace eos
{
    namespace nc
    {
        struct BToK
        {
            constexpr static const char * label = "B->K";
        };
        constexpr const char * BToK::label;
    }

    namespace nc_p_to_p
    {
        class Naive :
            public NonlocalFormFactor<nc::PToP>
        {
            public:
                Naive(const Parameters &, const Options &)
                {
                }

                ~Naive() = default;

                virtual complex<double> H_plus(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> Hhat_plus(const double & q2) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_A(const double & q2) const
                {
                    return 0.0;
                }

                static NonlocalFormFactorPtr<nc::PToP> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nc::PToP>(new Naive(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    return {};
                }
        };

        class LCSR :
            public NonlocalFormFactor<nc::PToP>
        {

                std::shared_ptr<Model> model;

                // B-meson parameters
                UsedParameter m_B;
                UsedParameter f_B;

                // final state meson parameters
                UsedParameter m_P;
                UsedParameter f_P;

                // sum rule parameters
                UsedParameter M2;
                UsedParameter s0_0_A;
                UsedParameter s0_1_A;

                // properties of the virtual quark in the sum rule's formfactor
                std::function<double ()> m_v;

                // renormalization scale for the sum rule
                UsedParameter mu_sr;

                // renormalization scale for the Wilson coefficients
                UsedParameter mu_ren;

                // renormalization scale for the virtual charm quark's mass
                UsedParameter mu_c;

                BMesonLCDAs b_lcdas;

            public:
                LCSR(const Parameters & p, const Options & o) :
                    model(Model::make(o.get("model", "SM"), p, o)),
                    m_B(p["mass::B_d"], *this),
                    f_B(p["decay-constant::B_d"], *this),
                    m_P(p["mass::K_d"], *this),
                    f_P(p["decay-constant::K_d"], *this),
                    M2(p["B->K::M^2@B-LCSR"], *this),
                    s0_0_A(p["B->K::s_0^+,0@B-LCSR"], *this),
                    s0_1_A(p["B->K::s_0^+,1@B-LCSR"], *this),
                    mu_sr(p["B->K::mu@B-LCSR"], *this),
                    mu_ren(p["b->sccbar::mu"], *this),
                    mu_c(p["b->sccbar::mu_c"], *this),
                    b_lcdas(p, o + Options{ { "q", "d" } })
                {
                    this->uses(b_lcdas);

                    m_v = std::bind(&LCSR::m_virtual_s, this);
                }

                ~LCSR() = default;

                // mass of the virtual quark in the LCSR setup; v = s
                double m_virtual_s() const
                {
                    return model->m_s_msbar(mu_sr());
                }

                // mass of the virtual charm quark
                double m_c() const
                {
                    return model->m_c_msbar(mu_c());
                }

                /* forwarding the LCDAs */
                inline
                double phi_3(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.phi_3(omega_1, omega_2);
                }

                inline
                double phi_bar_3(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.phi_bar_3(omega_1, omega_2);
                }

                inline
                double phi_double_bar_3(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.phi_double_bar_3(omega_1, omega_2);
                }

                inline
                double phi_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.phi_4(omega_1, omega_2);
                }

                inline
                double phi_bar_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.phi_bar_4(omega_1, omega_2);
                }

                inline
                double phi_double_bar_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.phi_double_bar_4(omega_1, omega_2);
                }

                inline
                double psi_bar_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.psi_bar_4(omega_1, omega_2);
                }

                inline
                double psi_double_bar_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.psi_double_bar_4(omega_1, omega_2);
                }

                inline
                double chi_bar_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.chi_bar_4(omega_1, omega_2);
                }

                inline
                double chi_double_bar_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.chi_double_bar_4(omega_1, omega_2);
                }

                /* auxilliary functions */

                double s(const double & sigma, const double & q2) const
                {
                    const double sigmabar = 1.0 - sigma;

                    return sigma * pow(m_B(), 2) + (pow(m_v(), 2) - sigma * q2) / sigmabar;
                }

                double sigma(const double & s, const double & q2) const
                {
                    const double m_B2 = pow(m_B(), 2);
                    const double m_v2 = pow(m_v(), 2);

                    return (m_B2 - q2 + s - std::sqrt(4.0 * (m_v2 - s) * m_B2 + pow(m_B2 - q2 + s, 2))) / (2.0 * m_B2);
                }

                double sigma_0(const double & q2, const double & s0_0, const double & s0_1) const
                {
                    const double s0 = s0_0 + s0_1 * q2;

                    return sigma(s0, q2);
                }

                /* Integrands */

                inline
                double I1_A_phi_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_3 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_3 = this->phi_3(omega_1, omega_2);

                    const double C_1 = (pow(m_B,-2.0)*pow(sigmabar,-2.0)*(m_B*(4.0*m_v3*t*tbar*(-2.0 + u) -
                                       4.0*m_B2*m_v*sigmabar2*t*tbar*(-2.0 + u) - 4.0*m_B3*sigmabar3*t*tbar*(-3.0 + 2.0*u) +
                                       4.0*m_v*q2*(-sigmabar + 2.0*u + sigma*(-2.0 - 6.0*t2 + 6.0*t)*u + t*(2.0 - 7.0*u) +
                                       t2*(-2.0 + 7.0*u)) - m_B*sigmabar*
                                       (-(4.0*m_v2*t*tbar*(-3.0 + 2.0*u)) + q2*(-1 + 2.0*sigma)*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u)))
                                       - 4.0*m_B*u2*sigmabar*pow(omega,2.0)*
                                       (-(4.0*m_v*(1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u)) +
                                       2.0*m_B*sigma*(-3.0 + 6.0*u + t*(6.0 - 20.0*u) + t2*(-6.0 + 20.0*u)) -
                                       3.0*m_B*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)) +
                                       2.0*omega*u*(4.0*m_v*(m_v2 - q2)*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_B3*sigmabar2*(3.0 - 6.0*u - 8.0*t2*(-3.0 + 4.0*u) + 8.0*t*(-3.0 + 4.0*u) +
                                       2.0*sigma*(-3.0 + 6.0*u + t*(6.0 - 20.0*u) + t2*(-6.0 + 20.0*u))) +
                                       4.0*m_B2*m_v*sigmabar*(1 - 2.0*u + t2*(5.0 - 8.0*u) + t*(-5.0 + 8.0*u) +
                                       sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       m_B*(q2*(-1 + 2.0*sigma)*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) +
                                       m_v2*(-(2.0*sigma*(-3.0 + 6.0*u + t*(6.0 - 20.0*u) + t2*(-6.0 + 20.0*u))) +
                                       3.0*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))))))/4.0;

                    return C_1 * phi_3;
                }

                inline
                double I1_A_phi_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_3 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*(-(m_B*sigmabar*(q2 + q2*(-2.0 + 8.0*sigma*t*tbar)*u +
                                       4.0*t*tbar*(5.0*m_B2*sigmabar2*(-1 + u) - m_v2*(1 + u) + 2.0*m_B*m_v*sigmabar*(-2.0 + u)))) +
                                       4.0*omega*u*(2.0*m_v2*t*(-tbar + 2.0*u - 2.0*t*u) +
                                       m_B2*sigmabar2*(1 - 2.0*u + t2*(2.0 - 8.0*u) + t*(-2.0 + 8.0*u) +
                                       sigma*(-5.0 + t*(10.0 - 32*u) + 10.0*u + 2.0*t2*(-5.0 + 16.0*u))) +
                                       2.0*m_B*m_v*sigmabar*(-((1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u)) +
                                       2.0*sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       sigma*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u +
                                       m_v2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)))) -
                                       4.0*u2*sigmabar*(4.0*m_v*t*tbar + m_B*
                                       (1 - 2.0*u - 4.0*t2*(4.0 - u + 2.0*sigma*(-2.0 + u)) + 4.0*t*(4.0 - u + 2.0*sigma*(-2.0 + u))))*
                                       pow(omega,2.0))*pow(sigmabar,-3.0))/4.0;

                    return C_1 * phi_bar_3;
                }

                inline
                double I2_A_phi_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_3 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

                    const double C_2 = (pow(m_B,-2.0)*pow(sigmabar,-4.0)*(2.0*omega*u*
                                       (m_B4*sigmabar4*(1 + t2*(40 - 36*u) - 2.0*u + 4.0*t*(-10.0 + 9.0*u) +
                                       2.0*sigma*(-5.0 + t*(10.0 - 32*u) + 10.0*u + 2.0*t2*(-5.0 + 16.0*u))) +
                                       4.0*m_B3*m_v*sigmabar3*(1 - 2.0*u + t2*(9.0 - 10.0*u) + t*(-9.0 + 10.0*u) +
                                       2.0*sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) -
                                       4.0*m_B*m_v*sigmabar*(q2*(-1 + t2 + 2.0*sigma - t + 2.0*u + 6.0*t2*u - 6.0*t*u -
                                       4.0*sigma*(1 + 3.0*t2 - 3.0*t)*u) +
                                       m_v2*(-((1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u)) +
                                       2.0*sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)))) -
                                       4.0*m_B2*sigmabar2*(t*tbar*(m_v2*(-3.0 + 2.0*u) + 2.0*q2*(-5.0 + 4.0*u)) +
                                       sigma*(2.0*m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       q2*(2.0 - 4.0*u - t2*(5.0 + 6.0*u) + t*(5.0 + 6.0*u)))) -
                                       (m_v2 - q2)*(-(q2*(-1 + 2.0*sigma)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       m_v2*(1 - 2.0*u + t2*(4.0 - 12.0*u) + 4.0*t*(-1 + 3.0*u) +
                                       2.0*sigma*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))))) -
                                       4.0*u2*sigmabar*(-(4.0*m_v*(m_v2 - q2)*t*tbar) +
                                       4.0*m_B2*m_v*sigmabar*(1 - 2.0*u + t2*(5.0 - 8.0*u) + t*(-5.0 + 8.0*u) +
                                       sigma*(-2.0 + 4.0*u + t*(5.0 - 12.0*u) + t2*(-5.0 + 12.0*u))) +
                                       m_B3*sigmabar2*(1 - 20.0*t2*(-1 + u) + 20.0*t*(-1 + u) - 2.0*u +
                                       2.0*sigma*(-5.0 + t*(10.0 - 32*u) + 10.0*u + 2.0*t2*(-5.0 + 16.0*u))) +
                                       m_B*(q2*(-1 + 20.0*t2*(-1 + u) - 20.0*t*(-1 + u) + 2.0*u +
                                       sigma*(2.0 - 4.0*u - 4.0*t2*(-5.0 + 6.0*u) + 4.0*t*(-5.0 + 6.0*u))) +
                                       m_v2*(1 - 2.0*u + t2*(4.0 - 12.0*u) + 4.0*t*(-1 + 3.0*u) +
                                       2.0*sigma*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)))))*pow(omega,2.0) -
                                       m_B*sigmabar*(4.0*m_v4*t*tbar + 8.0*m_B3*m_v*sigmabar3*t*tbar*(-2.0 + u) +
                                       4.0*m_B4*sigmabar4*t*tbar*(-5.0 + 4.0*u) +
                                       4.0*m_B*m_v*sigmabar*(-(2.0*m_v2*t*tbar*(-2.0 + u)) +
                                       q2*(-1 + 2.0*sigma)*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) +
                                       m_B2*sigmabar2*(-(16.0*m_v2*t*tbar*(-1 + u)) +
                                       q2*(1 - 2.0*u + 4.0*t2*(-5.0 + 3.0*u) - 4.0*t*(-5.0 + 3.0*u) +
                                       2.0*sigma*(-5.0 + 2.0*(5.0 + 12.0*t2 - 12.0*t)*u))) -
                                       (-1 + 2.0*sigma)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(q2,2.0) +
                                       m_v2*q2*(1 - 2.0*u + t2*(4.0 - 12.0*u) + 4.0*t*(-1 + 3.0*u) +
                                       2.0*sigma*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/4.0;

                    return C_2 * phi_bar_3;
                }

                inline
                double I1_A_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to phi_double_bar_3 with a single pole in k2
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_B2      = pow(m_B, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_1 = -(2.0*pow(m_B,-2.0)*(m_B2*sigmabar*t*tbar +
                                       m_B*omega*u*(2.0*t*tbar*u + sigma*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))) +
                                       4.0*u2*t*(-tbar + u - t*u)*pow(omega,2.0))*pow(sigmabar,-2.0));

                    return C_1 * phi_double_bar_3;
                }

                inline
                double I2_A_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_3 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_2 = -(pow(sigmabar,-3.0)*(2.0*omega*u*(m_B2*sigma3*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       2.0*(-(m_B*m_v*t*tbar) - (m_B2 - q2)*t*tbar*(-2.0 + u) + m_v2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) -
                                       2.0*m_B*sigma2*(m_B*(-1 + (2.0 + 3.0*t2 - 3.0*t)*u) +
                                       m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       sigma*(2.0*m_B*m_v*(1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u) + m_B2*(-1 + 6.0*t2 - 6.0*t + 2.0*u) +
                                       q2*(3.0 - 6.0*u + t2*(2.0 - 16.0*u) + 2.0*t*(-1 + 8.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))))*pow(m_B,-1) -
                                       4.0*u2*(-(2.0*(3.0*m_B2 + m_B*m_v + m_v2 - 3.0*q2)*t*tbar*(-1 + u)) +
                                       m_B2*sigma2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)) +
                                       m_B*sigma*(-(2.0*m_v*t*tbar) + m_B*(1 - 2.0*u - 4.0*t2*(-3.0 + 4.0*u) + 4.0*t*(-3.0 + 4.0*u))))
                                       *pow(m_B,-2.0)*pow(omega,2.0) - sigmabar*
                                       (2.0*m_B2*sigma2*t*(-tbar + 2.0*u - 2.0*t*u) -
                                       2.0*t*tbar*(m_v2*(1 + u) - m_B*m_v*(-2.0 + u) + m_B2*(1 - 2.0*u) + q2*(-1 + 2.0*u)) +
                                       sigma*(-(2.0*m_B*t*tbar*(m_v*(-2.0 + u) + m_B*(-2.0 + 4.0*u))) + q2*(1 - 2.0*u*pow(1 - 2.0*t,2.0))))));

                    return C_2 * phi_double_bar_3;
                }

                inline
                double I3_A_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_3 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_3 = 2.0*pow(sigmabar,-4.0)*(sigmabar*(-(2.0*m_B4*sigma4*t*tbar*(-1 + u)) - m_v4*t*tbar*u +
                                       m_B*m_v3*t*tbar*(-2.0 + u) - m_B*m_v*(m_B2 - q2)*t*tbar*(-2.0 + u) +
                                       m_B2*sigma3*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u + m_B*t*tbar*(8.0*m_B*(-1 + u) + m_v*(-2.0 + u))) +
                                       m_v2*(m_B2*t*tbar*(-2.0 + 3.0*u) + q2*(1 - 2.0*u + t2*(2.0 - 5.0*u) + t*(-2.0 + 5.0*u))) +
                                       m_B*sigma2*(-(12.0*m_B3*t*tbar*(-1 + u)) - 3.0*m_B2*m_v*t*tbar*(-2.0 + u) +
                                       m_v*q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       m_B*(m_v2*t*tbar*(-2.0 + 3.0*u) + q2*(-2.0 + 4.0*u + 4.0*t2*(1 + u) - 4.0*t*(1 + u)))) +
                                       sigma*(8.0*m_B4*t*tbar*(-1 + u) + 3.0*m_B3*m_v*t*tbar*(-2.0 + u) +
                                       q2*(2.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       m_B2*(-(2.0*m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                       q2*(1 - 2.0*u + 4.0*t2*(-2.0 + u) - 4.0*t*(-2.0 + u))) +
                                       m_B*(-(m_v3*t*tbar*(-2.0 + u)) + m_v*q2*(1 - 2.0*u - t2*(2.0 + 5.0*u) + t*(2.0 + 5.0*u)))) -
                                       2.0*t*tbar*(-1 + u)*pow(m_B2 - q2,2.0)) +
                                       2.0*omega*u*pow(m_B,-1)*(m_v4 - m_v2*q2 - 4.0*m_v2*t2*q2 + 4.0*m_v2*q2*t +
                                       4.0*m_B4*t*tbar*(-1 + u) - 2.0*m_v4*u - 5.0*m_v4*t2*u + 2.0*m_v2*q2*u + 7.0*m_v2*t2*q2*u +
                                       5.0*m_v4*t*u - 7.0*m_v2*q2*t*u + 2.0*m_B3*m_v*t*tbar*(-2.0 + u) -
                                       2.0*m_B*m_v*(m_v2 + q2)*t*tbar*(-2.0 + u) +
                                       m_B4*sigma5*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_B2*sigma3*(m_B*m_v*(3.0 + 10.0*t2 - 10.0*t)*(-1 + 2.0*u) + 4.0*q2*t*(-tbar + u - t*u) +
                                       2.0*m_B2*(-3.0 + t*(14.0 - 26*u) + 6.0*u + 2.0*t2*(-7.0 + 13.0*u)) +
                                       m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       m_B3*sigma4*(m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) -
                                       4.0*m_B*(-1 + 2.0*u + t*(3.0 - 7.0*u) + t2*(-3.0 + 7.0*u))) +
                                       m_B2*(-(8.0*q2*t*tbar*(-1 + u)) + m_v2*(-1 + 2.0*u + t*(4.0 - 9.0*u) + t2*(-4.0 + 9.0*u))) +
                                       m_B*sigma2*(-(4.0*m_B3*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u))) +
                                       m_v3*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_v*q2*(1 - 2.0*u - 2.0*t2*(1 + 3.0*u) + t*(2.0 + 6.0*u)) -
                                       3.0*m_B2*m_v*(-1 + 2.0*u + t*(6.0 - 8.0*u) + t2*(-6.0 + 8.0*u)) +
                                       m_B*(-(16.0*q2*t*tbar*(-1 + u)) + m_v2*(-3.0 + t*(8.0 - 21*u) + 6.0*u + t2*(-8.0 + 21*u)))) +
                                       sigma*(m_B4*(-1 + t*(18.0 - 22*u) + 2.0*u + 2.0*t2*(-9.0 + 11.0*u)) +
                                       m_B3*m_v*(-1 + 2.0*u + 2.0*t2*(-7.0 + 6.0*u) - 2.0*t*(-7.0 + 6.0*u)) +
                                       m_B2*(20.0*q2*t*tbar*(-1 + u) + m_v2*(3.0 + t2*(10.0 - 24*u) - 6.0*u + 2.0*t*(-5.0 + 12.0*u))) +
                                       q2*(q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) +
                                       m_B*(m_v*q2*(-1 + 2.0*u - 2.0*t*(3.0 + 2.0*u) + t2*(6.0 + 4.0*u)) +
                                       m_v3*(1 - 2.0*u + t2*(6.0 - 8.0*u) + t*(-6.0 + 8.0*u)))) + 4.0*t2*pow(q2,2.0) -
                                       4.0*t*pow(q2,2.0) - 4.0*t2*u*pow(q2,2.0) + 4.0*t*u*pow(q2,2.0)) +
                                       4.0*u2*pow(m_B,-2.0)*pow(omega,2.0)*
                                       (-(2.0*m_B4*t*tbar*(-1 + u)) + 2.0*(m_v2 - q2)*q2*t*tbar*(-1 + u) - m_B3*m_v*t*tbar*(-2.0 + u) +
                                       m_B*m_v*(m_v2 + q2)*t*tbar*(-2.0 + u) +
                                       m_B4*sigma4*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_B3*sigma3*(m_B*(3.0 - 6.0*u + t2*(8.0 - 20.0*u) + 4.0*t*(-2.0 + 5.0*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       m_B2*(4.0*q2*t*(-tbar + u - t*u) + m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       m_B2*sigma2*(2.0*m_v2*t*(-tbar + u - t*u) +
                                       m_B*m_v*(-2.0 + 4.0*u + t*(6.0 - 13.0*u) + t2*(-6.0 + 13.0*u)) +
                                       3.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0) - q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)) +
                                       m_B*sigma*(-(2.0*m_v*t*tbar*(-q2 + m_v2*(-1 + u))) +
                                       m_B3*(1 - 2.0*u + t2*(8.0 - 12.0*u) + 4.0*t*(-2.0 + 3.0*u)) +
                                       m_B2*m_v*(1 - 2.0*u + t2*(6.0 - 8.0*u) + t*(-6.0 + 8.0*u)) +
                                       m_B*(q2*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_v2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))));

                    return C_3 * phi_double_bar_3;
                }

                inline
                double I3d1_A_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_3 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);
                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_3 = 2.0*pow(sigmabar,-4.0)*(m_B*sigmabar*(-(2.0*m_B4*sigma4*t*tbar*(-1 + u)) - m_v4*t*tbar*u +
                                       m_B*m_v3*t*tbar*(-2.0 + u) - m_B*m_v*(m_B2 - q2)*t*tbar*(-2.0 + u) +
                                       m_B2*sigma3*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u + m_B*t*tbar*(8.0*m_B*(-1 + u) + m_v*(-2.0 + u))) +
                                       m_v2*(m_B2*t*tbar*(-2.0 + 3.0*u) + q2*(1 - 2.0*u + t2*(2.0 - 5.0*u) + t*(-2.0 + 5.0*u))) +
                                       m_B*sigma2*(-(12.0*m_B3*t*tbar*(-1 + u)) - 3.0*m_B2*m_v*t*tbar*(-2.0 + u) +
                                       m_v*q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       m_B*(m_v2*t*tbar*(-2.0 + 3.0*u) + q2*(-2.0 + 4.0*u + 4.0*t2*(1 + u) - 4.0*t*(1 + u)))) +
                                       sigma*(8.0*m_B4*t*tbar*(-1 + u) + 3.0*m_B3*m_v*t*tbar*(-2.0 + u) +
                                       q2*(2.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       m_B2*(-(2.0*m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                       q2*(1 - 2.0*u + 4.0*t2*(-2.0 + u) - 4.0*t*(-2.0 + u))) +
                                       m_B*(-(m_v3*t*tbar*(-2.0 + u)) + m_v*q2*(1 - 2.0*u - t2*(2.0 + 5.0*u) + t*(2.0 + 5.0*u)))) -
                                       2.0*t*tbar*(-1 + u)*pow(m_B2 - q2,2.0)) +
                                       2.0*omega*u*(m_v4 - m_v2*q2 - 4.0*m_v2*t2*q2 + 4.0*m_v2*q2*t + 4.0*m_B4*t*tbar*(-1 + u) -
                                       2.0*m_v4*u - 5.0*m_v4*t2*u + 2.0*m_v2*q2*u + 7.0*m_v2*t2*q2*u + 5.0*m_v4*t*u -
                                       7.0*m_v2*q2*t*u + 2.0*m_B3*m_v*t*tbar*(-2.0 + u) - 2.0*m_B*m_v*(m_v2 + q2)*t*tbar*(-2.0 + u) +
                                       m_B4*sigma5*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_B2*sigma3*(m_B*m_v*(3.0 + 10.0*t2 - 10.0*t)*(-1 + 2.0*u) + 4.0*q2*t*(-tbar + u - t*u) +
                                       2.0*m_B2*(-3.0 + t*(14.0 - 26*u) + 6.0*u + 2.0*t2*(-7.0 + 13.0*u)) +
                                       m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       m_B3*sigma4*(m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) -
                                       4.0*m_B*(-1 + 2.0*u + t*(3.0 - 7.0*u) + t2*(-3.0 + 7.0*u))) +
                                       m_B2*(-(8.0*q2*t*tbar*(-1 + u)) + m_v2*(-1 + 2.0*u + t*(4.0 - 9.0*u) + t2*(-4.0 + 9.0*u))) +
                                       m_B*sigma2*(-(4.0*m_B3*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u))) +
                                       m_v3*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_v*q2*(1 - 2.0*u - 2.0*t2*(1 + 3.0*u) + t*(2.0 + 6.0*u)) -
                                       3.0*m_B2*m_v*(-1 + 2.0*u + t*(6.0 - 8.0*u) + t2*(-6.0 + 8.0*u)) +
                                       m_B*(-(16.0*q2*t*tbar*(-1 + u)) + m_v2*(-3.0 + t*(8.0 - 21*u) + 6.0*u + t2*(-8.0 + 21*u)))) +
                                       sigma*(m_B4*(-1 + t*(18.0 - 22*u) + 2.0*u + 2.0*t2*(-9.0 + 11.0*u)) +
                                       m_B3*m_v*(-1 + 2.0*u + 2.0*t2*(-7.0 + 6.0*u) - 2.0*t*(-7.0 + 6.0*u)) +
                                       m_B2*(20.0*q2*t*tbar*(-1 + u) + m_v2*(3.0 + t2*(10.0 - 24*u) - 6.0*u + 2.0*t*(-5.0 + 12.0*u))) +
                                       q2*(q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) +
                                       m_B*(m_v*q2*(-1 + 2.0*u - 2.0*t*(3.0 + 2.0*u) + t2*(6.0 + 4.0*u)) +
                                       m_v3*(1 - 2.0*u + t2*(6.0 - 8.0*u) + t*(-6.0 + 8.0*u)))) + 4.0*t2*pow(q2,2.0) -
                                       4.0*t*pow(q2,2.0) - 4.0*t2*u*pow(q2,2.0) + 4.0*t*u*pow(q2,2.0)) +
                                       4.0*u2*pow(m_B,-1)*pow(omega,2.0)*(-(2.0*m_B4*t*tbar*(-1 + u)) + 2.0*(m_v2 - q2)*q2*t*tbar*(-1 + u) -
                                       m_B3*m_v*t*tbar*(-2.0 + u) + m_B*m_v*(m_v2 + q2)*t*tbar*(-2.0 + u) +
                                       m_B4*sigma4*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_B3*sigma3*(m_B*(3.0 - 6.0*u + t2*(8.0 - 20.0*u) + 4.0*t*(-2.0 + 5.0*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       m_B2*(4.0*q2*t*(-tbar + u - t*u) + m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       m_B2*sigma2*(2.0*m_v2*t*(-tbar + u - t*u) +
                                       m_B*m_v*(-2.0 + 4.0*u + t*(6.0 - 13.0*u) + t2*(-6.0 + 13.0*u)) +
                                       3.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0) - q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)) +
                                       m_B*sigma*(-(2.0*m_v*t*tbar*(-q2 + m_v2*(-1 + u))) +
                                       m_B3*(1 - 2.0*u + t2*(8.0 - 12.0*u) + 4.0*t*(-2.0 + 3.0*u)) +
                                       m_B2*m_v*(1 - 2.0*u + t2*(6.0 - 8.0*u) + t*(-6.0 + 8.0*u)) +
                                       m_B*(q2*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_v2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))));

                    const double C_3d1 = 4.0*m_B4*t*tbar*(-1 + u) + 4.0*m_B3*omega*u*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                         64*u2*(m_v2 - q2)*q2*t*tbar*(-1 + u)*pow(m_B,-2.0)*pow(omega,2.0)*pow(sigmabar,-5.0) +
                                         4.0*omega*u*pow(m_B,-1)*(-(4.0*m_v3*omega*t*tbar*(3.0 + 3.0*sigma*(-1 + u) - u)*u) +
                                         4.0*m_v*omega*q2*t*tbar*u*(-3.0 + 3.0*sigma + 2.0*u) - 4.0*m_v4*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                         m_v2*q2*(-5.0 + t*(14.0 - 32*u) + 10.0*u + 2.0*t2*(-7.0 + 16.0*u) +
                                         3.0*sigma*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) +
                                         (1 - 2.0*u - 2.0*t2*(-7.0 + 9.0*u) + 2.0*t*(-7.0 + 9.0*u) +
                                         sigma*(3.0 - 6.0*u - 6.0*t2*(1 + u) + 6.0*t*(1 + u)))*pow(q2,2.0))*pow(sigmabar,-5.0) -
                                         2.0*m_B2*pow(sigmabar,-2.0)*(-(m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                         2.0*m_v*omega*u*(1 - 2.0*u - 2.0*t2*(1 + 2.0*u) + t*(2.0 + 4.0*u)) +
                                         4.0*u2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(omega,2.0) + q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)) -
                                         2.0*m_B*pow(sigmabar,-3.0)*(-(8.0*omega*q2*(-3.0 + sigma)*t*tbar*(-1 + u)*u) - 2.0*m_v3*t*tbar*(-2.0 + u) +
                                         2.0*m_v2*omega*u*(-((1 + 6.0*t2 - 6.0*t)*(-1 + 2.0*u)) +
                                         sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                         m_v*(4.0*u2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u) +
                                         sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)))*pow(omega,2.0) +
                                         q2*(sigma*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) + (-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))) -
                                         2.0*pow(sigmabar,-4.0)*(3.0*m_v4*t*tbar*u +
                                         2.0*m_v3*omega*u*(-1 + 10.0*t2 - 10.0*t + 2.0*u +
                                         2.0*sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                         2.0*m_v*omega*q2*u*(1 - 2.0*u - 2.0*t2*(-5.0 + 6.0*u) + 2.0*t*(-5.0 + 6.0*u) -
                                         2.0*sigma*(-1 + 2.0*u - 2.0*t*(1 + 3.0*u) + t2*(2.0 + 6.0*u))) +
                                         m_v2*(q2*(-3.0 + 6.0*u + t*(6.0 - 17.0*u - 4.0*sigma*u) + t2*(-6.0 + 17.0*u + 4.0*sigma*u)) -
                                         4.0*u2*(3.0 - 6.0*u + 4.0*t2*(1 + sigma*(-1 + u) - 4.0*u) - 4.0*t*(1 + sigma*(-1 + u) - 4.0*u))*
                                         pow(omega,2.0)) + q2*(q2*(1 - 2.0*u + t2*(6.0 - 10.0*u) + 2.0*t*(-3.0 + 5.0*u) +
                                         sigma*(2.0 + (-4.0 - 8.0*t2 + 8.0*t)*u)) -
                                         4.0*u2*pow(omega,2.0)*(-1 + 2.0*u - 4.0*t2*(-2.0 + u) + 4.0*t*(-2.0 + u) +
                                         2.0*sigma*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))));

                    return C_3 * phi_bar_3 + C_3d1 * phi_double_bar_3;
                }

                inline
                double I1_A_phi_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_4 = this->phi_4(omega_1, omega_2);

                    const double C_1 = (pow(sigmabar,-2.0)*(-(sigmabar*(-(4.0*(-m_v2 + m_B2*sigmabar2)*t*tbar*(-1 + u)) +
                                       q2*(-3.0 + 6.0*u + t*(8.0 - 20.0*u) + 4.0*t2*(-2.0 + 5.0*u) +
                                       sigma*(2.0 + (-4.0 - 8.0*t2 + 8.0*t)*u)))) -
                                       4.0*u2*sigmabar*pow(m_B,-1)*pow(omega,2.0)*
                                       (-(4.0*m_v*t*tbar*(-3.0 + 2.0*u)) -
                                       m_B*(-1 + 2.0*u)*(2.0*sigma*(1 + 6.0*t2 - 6.0*t) - 3.0*pow(1 - 2.0*t,2.0))) +
                                       2.0*omega*u*pow(m_B,-1)*(-(4.0*m_B*m_v*sigmabar*t*tbar*(-3.0 + 2.0*u)) +
                                       m_v2*(3.0 - 6.0*u + t2*(4.0 - 20.0*u) + 4.0*t*(-1 + 5.0*u) +
                                       2.0*sigma*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))) -
                                       m_B2*sigmabar2*(3.0 - 6.0*u - 4.0*t2*(-4.0 + 7.0*u) + 4.0*t*(-4.0 + 7.0*u) +
                                       2.0*sigma*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))) +
                                       q2*(-3.0 + 6.0*u + t*(8.0 - 20.0*u) + 4.0*t2*(-2.0 + 5.0*u) -
                                       2.0*sigma*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))))/4.0;

                    return C_1 * phi_4;
                }

                inline
                double I1_A_phi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*(-(m_B*sigmabar*(-(4.0*t*tbar*
                                       (5.0*m_B2*sigmabar2*(-1 + u) + m_B*m_v*sigmabar*(-2.0 + u) + m_v2*(1 - 2.0*u))) +
                                       q2*(-3.0 + 4.0*sigma + 6.0*u - 8.0*sigma*(1 + 2.0*t2 - 2.0*t)*u + t*(8.0 - 20.0*u) +
                                       4.0*t2*(-2.0 + 5.0*u)))) - 4.0*omega*u*
                                       (2.0*m_v2*t*tbar + 2.0*q2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_B2*sigmabar2*(1 - 2.0*u - 2.0*t2*(-5.0 + 8.0*u) + 2.0*t*(-5.0 + 8.0*u) +
                                       sigma*(-5.0 + t*(10.0 - 32*u) + 10.0*u + 2.0*t2*(-5.0 + 16.0*u))) +
                                       sigma*(-(m_v2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)) +
                                       q2*(-3.0 + 6.0*u + t*(4.0 - 16.0*u) + 4.0*t2*(-1 + 4.0*u))) +
                                       2.0*m_B*m_v*sigmabar*(-(t*tbar*(-1 + u)) + sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)))
                                       ) - 4.0*u2*sigmabar*(4.0*m_v*t*tbar*(-3.0 + 2.0*u) +
                                       m_B*(-3.0 + 4.0*sigma - 8.0*sigma*(1 + 2.0*t2 - 2.0*t)*u + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u))*
                                       pow(omega,2.0))*pow(sigmabar,-3.0))/4.0;

                    return C_1 * phi_bar_4;
                }

                inline
                double I2_A_phi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

                    const double C_2 = (pow(m_B,-2.0)*pow(sigmabar,-4.0)*(2.0*omega*u*
                                       (-(m_B4*sigmabar4*(-1 + 2.0*u - 4.0*t2*(-6.0 + 5.0*u) + 4.0*t*(-6.0 + 5.0*u) +
                                       2.0*sigma*(-3.0 + 6.0*u + t*(6.0 - 20.0*u) + t2*(-6.0 + 20.0*u)))) -
                                       4.0*m_B3*m_v*sigmabar3*(t*tbar*(-5.0 + 3.0*u) +
                                       sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       4.0*m_B*m_v*sigmabar*(-(t*tbar*(m_v2*(-1 + u) + q2*(5.0 - 3.0*u))) +
                                       sigma*(q2*(1 - 2.0*u - 2.0*t2*(1 + 2.0*u) + t*(2.0 + 4.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)))) +
                                       4.0*m_B2*sigmabar2*(q2*(6.0*t*(-tbar + u - t*u) +
                                       sigma*(2.0 - 4.0*u - t2*(3.0 + 4.0*u) + t*(3.0 + 4.0*u))) +
                                       m_v2*(-tbar + 2.0*u - 4.0*t*u + t2*(-1 + 4.0*u) +
                                       sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)))) +
                                       (m_v2 - q2)*(q2*(-1 + 2.0*sigma)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_v2*(3.0 - 6.0*u + t2*(4.0 - 20.0*u) + 4.0*t*(-1 + 5.0*u) +
                                       2.0*sigma*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))))) -
                                       4.0*u2*sigmabar*pow(omega,2.0)*(-(4.0*m_v*(m_v2 - q2)*t*tbar*(-3.0 + 2.0*u)) -
                                       4.0*m_B2*m_v*sigmabar*(sigma*(1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u) + t*tbar*(-3.0 + 2.0*u)) -
                                       m_B3*sigmabar2*(-1 + 2.0*u + 4.0*t*(-3.0 + 2.0*u) + t2*(12.0 - 8.0*u) +
                                       2.0*sigma*(-3.0 + 6.0*u + t*(6.0 - 20.0*u) + t2*(-6.0 + 20.0*u))) -
                                       m_B*(q2*(1 - 2.0*u + 4.0*t2*(-3.0 + 2.0*u) + t*(12.0 - 8.0*u) -
                                       2.0*sigma*(1 - 2.0*u + t*(6.0 - 4.0*u) + t2*(-6.0 + 4.0*u))) +
                                       m_v2*(-1 + 2.0*u)*(2.0*sigma*(1 + 6.0*t2 - 6.0*t) - 3.0*pow(1 - 2.0*t,2.0)))) +
                                       m_B*sigmabar*(4.0*m_v4*t*tbar + 12.0*m_B4*sigmabar4*t*tbar*(-1 + u) +
                                       4.0*m_B3*m_v*sigmabar3*t*tbar*(-2.0 + u) +
                                       4.0*m_B*m_v*sigmabar*(-((m_v2 + q2)*t*tbar*(-2.0 + u)) + q2*sigma*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) +
                                       m_B2*sigmabar2*(-(4.0*m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                       q2*(-1 + 2.0*u + 4.0*t2*(-3.0 + 4.0*u) - 4.0*t*(-3.0 + 4.0*u) +
                                       6.0*sigma*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))) +
                                       (-1 + 2.0*sigma)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(q2,2.0) +
                                       m_v2*q2*(3.0 - 6.0*u + t2*(4.0 - 20.0*u) + 4.0*t*(-1 + 5.0*u) +
                                       2.0*sigma*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/4.0;

                    return C_2 * phi_bar_4;
                }

                inline
                double I1_A_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to phi_double_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double sigma2    = pow(sigma, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*(m_B*sigmabar*t*tbar*(-(2.0*m_B*sigmabar) + m_v*(-2.0 + u)) +
                                       2.0*omega*u*(-(2.0*t*tbar*(m_v*(-1 + u) + m_B*u)) +
                                       m_B*sigma2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)) +
                                       sigma*(m_B*(1 - 2.0*u + t2*(2.0 - 10.0*u) + 2.0*t*(-1 + 5.0*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)))) -
                                       8.0*u2*sigmabar*t*tbar*(-1 + u)*pow(omega,2.0))*pow(sigmabar,-3.0));

                    return C_1 * phi_double_bar_4;
                }

                inline
                double I2_A_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_2 = pow(m_B,-2.0)*pow(sigmabar,-4.0)*(-(2.0*omega*u*
                                       (2.0*m_B3*t*tbar*(-2.0 + u) + m_B3*sigma4*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_B2*m_v*(-1 + 2.0*(1 + t2 - t)*u) + m_v3*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v*q2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       2.0*m_B*(-(q2*t*tbar*(-2.0 + u)) + m_v2*(1 + (-2.0 - 4.0*t2 + 4.0*t)*u)) +
                                       m_B2*sigma3*(m_B*(3.0 - 6.0*u + t2*(2.0 - 10.0*u) + 2.0*t*(-1 + 5.0*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       m_B*sigma2*(m_B2*(-3.0 + 6.0*u + 6.0*t2*(1 + u) - 6.0*t*(1 + u)) +
                                       m_B*m_v*(-3.0 + 6.0*u + t*(4.0 - 14.0*u) + 2.0*t2*(-2.0 + 7.0*u)) +
                                       q2*(3.0 - 6.0*u + t2*(2.0 - 16.0*u) + 2.0*t*(-1 + 8.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))) +
                                       sigma*(m_B3*(1 - 2.0*u + 2.0*t2*(-5.0 + u) - 2.0*t*(-5.0 + u)) +
                                       m_B2*m_v*(3.0 - 6.0*u + t2*(2.0 - 10.0*u) + 2.0*t*(-1 + 5.0*u)) +
                                       m_v3*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_v*q2*(-3.0 + 6.0*u + t*(2.0 - 14.0*u) + 2.0*t2*(-1 + 7.0*u)) +
                                       m_B*(m_v2*(-1 + 2.0*t2 - 2.0*t + 2.0*u) +
                                       q2*(-3.0 + 6.0*u + 2.0*t2*(1 + 7.0*u) - 2.0*t*(1 + 7.0*u)))))) -
                                       4.0*u2*sigmabar*(-(2.0*(3.0*m_B2 + m_B*m_v + m_v2 - 3.0*q2)*t*tbar*(-1 + u)) +
                                       m_B2*sigma2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)) +
                                       m_B*sigma*(-(m_v*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)) +
                                       m_B*(1 - 2.0*u - 4.0*t2*(-3.0 + 4.0*u) + 4.0*t*(-3.0 + 4.0*u))))*pow(omega,2.0) -
                                       m_B*sigmabar*(-(2.0*m_B3*sigma3*t*tbar*(-1 + 2.0*u)) +
                                       t*tbar*(m_B2*m_v*(-2.0 + u) + m_v3*(-2.0 + u) + m_v*q2*(2.0 - 3.0*u) + m_B3*(-2.0 + 4.0*u) -
                                       2.0*m_B*(m_v2*(1 + u) + q2*(-1 + 2.0*u))) +
                                       sigma*(m_v*q2 + m_v*q2*(-2.0 - 6.0*t2 + 6.0*t)*u + 2.0*m_B*m_v2*t*tbar*(1 + u) -
                                       2.0*m_B2*m_v*t*tbar*(-2.0 + u) - 6.0*m_B3*t*tbar*(-1 + 2.0*u) +
                                       m_B*q2*(1 - 2.0*u + t2*(2.0 - 12.0*u) + 2.0*t*(-1 + 6.0*u))) -
                                       m_B*sigma2*(-(m_B*t*tbar*(m_v*(-2.0 + u) + 6.0*m_B*(-1 + 2.0*u))) +
                                       q2*(1 - 2.0*u*pow(1 - 2.0*t,2.0)))));

                    return C_2 * phi_double_bar_4;
                }

                inline
                double I3_A_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4), m_B5 = pow(m_B, 5);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5), sigma6 = pow(sigma, 6);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(pow(m_B,-2.0)*pow(sigmabar,-5.0)*(-(m_B*sigmabar*
                                       (4.0*m_B5*t*tbar*(-1 + u) - 4.0*m_B5*sigma5*t*tbar*(-1 + u) + 2.0*m_B4*m_v*t*tbar*(-2.0 + u) +
                                       m_v*q2*(m_v2 + q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       2.0*m_B3*t*tbar*(4.0*q2*(-1 + u) + m_v2*(-2.0 + 3.0*u)) +
                                       2.0*m_B3*sigma4*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u +
                                       m_B*t*tbar*(10.0*m_B*(-1 + u) + m_v*(-2.0 + u))) +
                                       m_B2*m_v*(-(2.0*m_v2*t*tbar*(-2.0 + u)) +
                                       q2*(1 - 2.0*u - 2.0*t2*(2.0 + u) + 2.0*t*(2.0 + u))) +
                                       2.0*m_B2*sigma3*(-(20.0*m_B3*t*tbar*(-1 + u)) - 4.0*m_B2*m_v*t*tbar*(-2.0 + u) +
                                       m_v*q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       m_B*(m_v2*t*tbar*(-2.0 + 3.0*u) + q2*(-3.0 + 6.0*u - 4.0*t*(1 + 2.0*u) + t2*(4.0 + 8.0*u))))
                                       - 2.0*m_B*(-(m_v4*t*tbar*u) + m_v2*q2*(1 - 2.0*u + t2*(2.0 - 5.0*u) + t*(-2.0 + 5.0*u)) -
                                       2.0*t*tbar*(-1 + u)*pow(q2,2.0)) +
                                       m_B*sigma2*(40*m_B4*t*tbar*(-1 + u) + 4.0*m_v2*q2*t*tbar*u + 12.0*m_B3*m_v*t*tbar*(-2.0 + u) -
                                       2.0*m_B*m_v3*t*tbar*(-2.0 + u) - 6.0*m_B2*q2*(-1 + 4.0*t2 - 4.0*t + 2.0*u) -
                                       6.0*m_B2*m_v2*t*tbar*(-2.0 + 3.0*u) +
                                       m_B*m_v*q2*(5.0 - 10.0*u + t*(4.0 + 26*u) - 2.0*t2*(2.0 + 13.0*u)) +
                                       2.0*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(q2,2.0)) +
                                       2.0*sigma*(-(10.0*m_B5*t*tbar*(-1 + u)) - 4.0*m_B4*m_v*t*tbar*(-2.0 + u) -
                                       2.0*m_B2*(-(m_v3*t*tbar*(-2.0 + u)) +
                                       m_v*q2*(1 - 2.0*u - 2.0*t2*(1 + 2.0*u) + t*(2.0 + 4.0*u))) -
                                       m_B3*(-(3.0*m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                       q2*(1 - 2.0*u + 4.0*t2*(-3.0 + 2.0*u) + t*(12.0 - 8.0*u))) +
                                       m_v*(1 + (-2.0 - 4.0*t2 + 4.0*t)*u)*pow(q2,2.0) +
                                       m_B*(-(m_v4*t*tbar*u) + m_v2*q2*(1 - 2.0*u + t2*(2.0 - 3.0*u) + t*(-2.0 + 3.0*u)) +
                                       (1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u))*pow(q2,2.0))))) -
                                       2.0*omega*u*(-(8.0*m_B5*t*tbar*(-1 + u)) + m_B4*m_v*(-1 + 8.0*t2*(-1 + u) - 8.0*t*(-1 + u) + 2.0*u) +
                                       2.0*m_B5*sigma6*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       2.0*m_B4*sigma5*(m_B*(-5.0 + t*(14.0 - 34*u) + 10.0*u + 2.0*t2*(-7.0 + 17.0*u)) +
                                       m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) -
                                       2.0*m_B2*m_v*(-(2.0*q2*t*tbar*(-2.0 + u)) +
                                       m_v2*(-1 + 2.0*u + t*(4.0 - 6.0*u) + t2*(-4.0 + 6.0*u))) -
                                       2.0*m_B3*(-(8.0*q2*t*tbar*(-1 + u)) + m_v2*(-1 + 2.0*u + t*(4.0 - 9.0*u) + t2*(-4.0 + 9.0*u))) -
                                       2.0*m_B2*sigma3*(-(2.0*m_v*q2*t*tbar*(1 + u)) +
                                       4.0*m_B2*m_v*(-2.0 + 4.0*u + t*(7.0 - 13.0*u) + t2*(-7.0 + 13.0*u)) +
                                       10.0*m_B3*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)) +
                                       m_v3*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_B*(20.0*q2*t*tbar*(-1 + u) + m_v2*(4.0 + t2*(10.0 - 27*u) - 8.0*u + t*(-10.0 + 27*u)))) +
                                       2.0*m_B*sigma2*(m_B3*m_v*(-7.0 + t*(32 - 48*u) + 14.0*u + 16.0*t2*(-2.0 + 3.0*u)) +
                                       5.0*m_B4*(-1 + 2.0*u + 2.0*t2*(-5.0 + 7.0*u) - 2.0*t*(-5.0 + 7.0*u)) +
                                       q2*(q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) -
                                       3.0*m_B2*(-(12.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-2.0 + 4.0*u + t*(6.0 - 15.0*u) + 3.0*t2*(-2.0 + 5.0*u))) +
                                       m_B*(-(2.0*m_v*q2*t*tbar*(4.0 + u)) +
                                       m_v3*(3.0 - 6.0*u + t2*(8.0 - 18.0*u) + 2.0*t*(-4.0 + 9.0*u)))) +
                                       m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*(m_v4 - pow(q2,2.0)) +
                                       2.0*m_B*(m_v4*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                       m_v2*q2*(1 - 2.0*u + t2*(4.0 - 7.0*u) + t*(-4.0 + 7.0*u)) - 4.0*t*tbar*(-1 + u)*pow(q2,2.0)) -
                                       2.0*sigma*(m_v*(m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_B5*(-1 + t*(22 - 26*u) + 2.0*u + t2*(-22 + 26*u)) +
                                       m_B4*m_v*(-3.0 + t*(18.0 - 22*u) + 6.0*u + 2.0*t2*(-9.0 + 11.0*u)) +
                                       m_B3*(28*q2*t*tbar*(-1 + u) + m_v2*(4.0 + t2*(14.0 - 33*u) - 8.0*u + t*(-14.0 + 33*u))) +
                                       m_B2*(2.0*m_v*q2*t*tbar*(-5.0 + u) +
                                       m_v3*(3.0 - 6.0*u - 2.0*t2*(-5.0 + 9.0*u) + 2.0*t*(-5.0 + 9.0*u))) +
                                       m_B*(3.0*m_v2*q2*t*tbar*(-2.0 + u) + m_v4*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                       (1 - 2.0*u + 2.0*t2*(-3.0 + u) - 2.0*t*(-3.0 + u))*pow(q2,2.0))) +
                                       m_B3*sigma4*(8.0*q2*t*(-tbar + u - t*u) +
                                       2.0*m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_B*m_v*(-9.0 + 18.0*u + 8.0*t2*(-3.0 + 7.0*u) - 8.0*t*(-3.0 + 7.0*u)) +
                                       20.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))) +
                                       4.0*u2*sigmabar*pow(omega,2.0)*(-(4.0*m_B4*t*tbar*(-1 + u)) + 4.0*(m_v2 - q2)*q2*t*tbar*(-1 + u) +
                                       2.0*m_B4*sigma4*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_B3*m_v*(-1 + 2.0*u + t*(4.0 - 6.0*u) + t2*(-4.0 + 6.0*u)) -
                                       m_B*m_v*(m_v2 + q2)*(-1 + 2.0*u + t*(4.0 - 6.0*u) + t2*(-4.0 + 6.0*u)) -
                                       2.0*m_B3*sigma3*(m_B*(-3.0 + 6.0*u + t*(8.0 - 20.0*u) + 4.0*t2*(-2.0 + 5.0*u)) +
                                       m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) -
                                       2.0*m_B2*(-(4.0*q2*t*tbar*(-1 + u)) + m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       m_B2*sigma2*(4.0*m_v2*t*(-tbar + u - t*u) +
                                       m_B*m_v*(-5.0 + 10.0*u + 6.0*t2*(-2.0 + 5.0*u) - 6.0*t*(-2.0 + 5.0*u)) +
                                       6.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0) - 2.0*q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)) -
                                       2.0*m_B*sigma*(2.0*m_v3*t*tbar*(-1 + u) - m_v*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       2.0*m_B2*m_v*(1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u) +
                                       m_B3*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_B*(q2*(1 - 2.0*u + t2*(8.0 - 12.0*u) + 4.0*t*(-2.0 + 3.0*u)) -
                                       m_v2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))))));

                    return C_3 * phi_double_bar_4;
                }

                inline
                double I3d1_A_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4), m_B5 = pow(m_B, 5);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5), sigma6 = pow(sigma, 6);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);
                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(pow(sigmabar,-5.0)*(-(sigmabar*(4.0*m_B5*t*tbar*(-1 + u) - 4.0*m_B5*sigma5*t*tbar*(-1 + u) +
                                       2.0*m_B4*m_v*t*tbar*(-2.0 + u) + m_v*q2*(m_v2 + q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       2.0*m_B3*t*tbar*(4.0*q2*(-1 + u) + m_v2*(-2.0 + 3.0*u)) +
                                       2.0*m_B3*sigma4*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u +
                                       m_B*t*tbar*(10.0*m_B*(-1 + u) + m_v*(-2.0 + u))) +
                                       m_B2*m_v*(-(2.0*m_v2*t*tbar*(-2.0 + u)) +
                                       q2*(1 - 2.0*u - 2.0*t2*(2.0 + u) + 2.0*t*(2.0 + u))) +
                                       2.0*m_B2*sigma3*(-(20.0*m_B3*t*tbar*(-1 + u)) - 4.0*m_B2*m_v*t*tbar*(-2.0 + u) +
                                       m_v*q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       m_B*(m_v2*t*tbar*(-2.0 + 3.0*u) + q2*(-3.0 + 6.0*u - 4.0*t*(1 + 2.0*u) + t2*(4.0 + 8.0*u))))
                                       - 2.0*m_B*(-(m_v4*t*tbar*u) + m_v2*q2*(1 - 2.0*u + t2*(2.0 - 5.0*u) + t*(-2.0 + 5.0*u)) -
                                       2.0*t*tbar*(-1 + u)*pow(q2,2.0)) +
                                       m_B*sigma2*(40*m_B4*t*tbar*(-1 + u) + 4.0*m_v2*q2*t*tbar*u + 12.0*m_B3*m_v*t*tbar*(-2.0 + u) -
                                       2.0*m_B*m_v3*t*tbar*(-2.0 + u) - 6.0*m_B2*q2*(-1 + 4.0*t2 - 4.0*t + 2.0*u) -
                                       6.0*m_B2*m_v2*t*tbar*(-2.0 + 3.0*u) +
                                       m_B*m_v*q2*(5.0 - 10.0*u + t*(4.0 + 26*u) - 2.0*t2*(2.0 + 13.0*u)) +
                                       2.0*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(q2,2.0)) +
                                       2.0*sigma*(-(10.0*m_B5*t*tbar*(-1 + u)) - 4.0*m_B4*m_v*t*tbar*(-2.0 + u) -
                                       2.0*m_B2*(-(m_v3*t*tbar*(-2.0 + u)) +
                                       m_v*q2*(1 - 2.0*u - 2.0*t2*(1 + 2.0*u) + t*(2.0 + 4.0*u))) -
                                       m_B3*(-(3.0*m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                       q2*(1 - 2.0*u + 4.0*t2*(-3.0 + 2.0*u) + t*(12.0 - 8.0*u))) +
                                       m_v*(1 + (-2.0 - 4.0*t2 + 4.0*t)*u)*pow(q2,2.0) +
                                       m_B*(-(m_v4*t*tbar*u) + m_v2*q2*(1 - 2.0*u + t2*(2.0 - 3.0*u) + t*(-2.0 + 3.0*u)) +
                                       (1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u))*pow(q2,2.0))))) -
                                       2.0*omega*u*pow(m_B,-1)*(-(8.0*m_B5*t*tbar*(-1 + u)) +
                                       m_B4*m_v*(-1 + 8.0*t2*(-1 + u) - 8.0*t*(-1 + u) + 2.0*u) +
                                       2.0*m_B5*sigma6*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       2.0*m_B4*sigma5*(m_B*(-5.0 + t*(14.0 - 34*u) + 10.0*u + 2.0*t2*(-7.0 + 17.0*u)) +
                                       m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) -
                                       2.0*m_B2*m_v*(-(2.0*q2*t*tbar*(-2.0 + u)) +
                                       m_v2*(-1 + 2.0*u + t*(4.0 - 6.0*u) + t2*(-4.0 + 6.0*u))) -
                                       2.0*m_B3*(-(8.0*q2*t*tbar*(-1 + u)) + m_v2*(-1 + 2.0*u + t*(4.0 - 9.0*u) + t2*(-4.0 + 9.0*u))) -
                                       2.0*m_B2*sigma3*(-(2.0*m_v*q2*t*tbar*(1 + u)) +
                                       4.0*m_B2*m_v*(-2.0 + 4.0*u + t*(7.0 - 13.0*u) + t2*(-7.0 + 13.0*u)) +
                                       10.0*m_B3*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)) +
                                       m_v3*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_B*(20.0*q2*t*tbar*(-1 + u) + m_v2*(4.0 + t2*(10.0 - 27*u) - 8.0*u + t*(-10.0 + 27*u)))) +
                                       2.0*m_B*sigma2*(m_B3*m_v*(-7.0 + t*(32 - 48*u) + 14.0*u + 16.0*t2*(-2.0 + 3.0*u)) +
                                       5.0*m_B4*(-1 + 2.0*u + 2.0*t2*(-5.0 + 7.0*u) - 2.0*t*(-5.0 + 7.0*u)) +
                                       q2*(q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) -
                                       3.0*m_B2*(-(12.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-2.0 + 4.0*u + t*(6.0 - 15.0*u) + 3.0*t2*(-2.0 + 5.0*u))) +
                                       m_B*(-(2.0*m_v*q2*t*tbar*(4.0 + u)) +
                                       m_v3*(3.0 - 6.0*u + t2*(8.0 - 18.0*u) + 2.0*t*(-4.0 + 9.0*u)))) +
                                       m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*(m_v4 - pow(q2,2.0)) +
                                       2.0*m_B*(m_v4*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                       m_v2*q2*(1 - 2.0*u + t2*(4.0 - 7.0*u) + t*(-4.0 + 7.0*u)) - 4.0*t*tbar*(-1 + u)*pow(q2,2.0)) -
                                       2.0*sigma*(m_v*(m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_B5*(-1 + t*(22 - 26*u) + 2.0*u + t2*(-22 + 26*u)) +
                                       m_B4*m_v*(-3.0 + t*(18.0 - 22*u) + 6.0*u + 2.0*t2*(-9.0 + 11.0*u)) +
                                       m_B3*(28*q2*t*tbar*(-1 + u) + m_v2*(4.0 + t2*(14.0 - 33*u) - 8.0*u + t*(-14.0 + 33*u))) +
                                       m_B2*(2.0*m_v*q2*t*tbar*(-5.0 + u) +
                                       m_v3*(3.0 - 6.0*u - 2.0*t2*(-5.0 + 9.0*u) + 2.0*t*(-5.0 + 9.0*u))) +
                                       m_B*(3.0*m_v2*q2*t*tbar*(-2.0 + u) + m_v4*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                       (1 - 2.0*u + 2.0*t2*(-3.0 + u) - 2.0*t*(-3.0 + u))*pow(q2,2.0))) +
                                       m_B3*sigma4*(8.0*q2*t*(-tbar + u - t*u) +
                                       2.0*m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_B*m_v*(-9.0 + 18.0*u + 8.0*t2*(-3.0 + 7.0*u) - 8.0*t*(-3.0 + 7.0*u)) +
                                       20.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))) +
                                       4.0*u2*sigmabar*pow(m_B,-1)*pow(omega,2.0)*
                                       (-(4.0*m_B4*t*tbar*(-1 + u)) + 4.0*(m_v2 - q2)*q2*t*tbar*(-1 + u) +
                                       2.0*m_B4*sigma4*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_B3*m_v*(-1 + 2.0*u + t*(4.0 - 6.0*u) + t2*(-4.0 + 6.0*u)) -
                                       m_B*m_v*(m_v2 + q2)*(-1 + 2.0*u + t*(4.0 - 6.0*u) + t2*(-4.0 + 6.0*u)) -
                                       2.0*m_B3*sigma3*(m_B*(-3.0 + 6.0*u + t*(8.0 - 20.0*u) + 4.0*t2*(-2.0 + 5.0*u)) +
                                       m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) -
                                       2.0*m_B2*(-(4.0*q2*t*tbar*(-1 + u)) + m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       m_B2*sigma2*(4.0*m_v2*t*(-tbar + u - t*u) +
                                       m_B*m_v*(-5.0 + 10.0*u + 6.0*t2*(-2.0 + 5.0*u) - 6.0*t*(-2.0 + 5.0*u)) +
                                       6.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0) - 2.0*q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)) -
                                       2.0*m_B*sigma*(2.0*m_v3*t*tbar*(-1 + u) - m_v*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       2.0*m_B2*m_v*(1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u) +
                                       m_B3*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_B*(q2*(1 - 2.0*u + t2*(8.0 - 12.0*u) + 4.0*t*(-2.0 + 3.0*u)) -
                                       m_v2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))))));

                    const double C_3d1 = 2.0*(-(2.0*m_B4*t*tbar*(-1 + u)) - 2.0*m_B3*omega*u*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                         m_B*(24*omega*q2*t*tbar*(-1 + u)*u - 2.0*m_v3*t*tbar*(-2.0 + u) -
                                         2.0*m_v2*omega*(1 + 6.0*t2 - 6.0*t)*u*(-1 + 2.0*u) -
                                         4.0*m_v*t*tbar*(q2*(-1 + u) + 2.0*u2*pow(omega,2.0)) +
                                         sigma*(-(8.0*omega*q2*t*tbar*(-1 + u)*u) +
                                         2.0*m_v2*omega*u*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                         m_v*(q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                         4.0*u2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))*pow(omega,2.0))))*
                                         pow(sigmabar,-3.0) - pow(m_B,-1)*(-(8.0*m_v4*omega*u*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) +
                                         2.0*m_v2*omega*q2*u*(-5.0 + t*(14.0 - 32*u) + 10.0*u + 2.0*t2*(-7.0 + 16.0*u) +
                                         3.0*sigma*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) +
                                         2.0*m_v3*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u +
                                         4.0*u2*(1 - 2.0*u + t2*(3.0 + 3.0*sigma*(-1 + u) - 5.0*u) +
                                         t*(-3.0 - 3.0*sigma*(-1 + u) + 5.0*u))*pow(omega,2.0)) +
                                         m_v*q2*(q2*(-1 + 3.0*sigma)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                         4.0*u2*(1 - 2.0*u + 3.0*sigma*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) + t2*(6.0 - 8.0*u) +
                                         t*(-6.0 + 8.0*u))*pow(omega,2.0)) -
                                         2.0*omega*u*(-1 + 2.0*u + 2.0*t2*(-7.0 + 9.0*u) - 2.0*t*(-7.0 + 9.0*u) +
                                         sigma*(-3.0 + 6.0*u + 6.0*t2*(1 + u) - 6.0*t*(1 + u)))*pow(q2,2.0))*pow(sigmabar,-5.0) +
                                         omega*(m_v2 - q2)*u*(-32*omega*q2*sigmabar*t*tbar*(-1 + u)*u + 5.0*m_v3*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                         m_v*q2*(-3.0 + 8.0*sigma)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))*pow(m_B,-2.0)*pow(sigmabar,-6.0) +
                                         m_B2*pow(sigmabar,-2.0)*(-(m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                         m_v*omega*u*(1 - 2.0*u - 4.0*t2*(1 + u) + 4.0*t*(1 + u)) +
                                         4.0*u2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(omega,2.0) + q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)) +
                                         pow(sigmabar,-4.0)*(3.0*m_v4*t*tbar*u + 4.0*m_v*omega*q2*t*tbar*u*(-5.0 + 4.0*u + 2.0*sigma*(1 + u)) +
                                         4.0*m_v3*omega*u*(1 - 2.0*u + t2*(5.0 - 6.0*u) + t*(-5.0 + 6.0*u) +
                                         sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                         m_v2*(q2*(-3.0 + 6.0*u + t*(6.0 - 17.0*u - 4.0*sigma*u) + t2*(-6.0 + 17.0*u + 4.0*sigma*u)) -
                                         4.0*u2*(3.0 - 6.0*u + 4.0*t2*(1 + sigma*(-1 + u) - 4.0*u) - 4.0*t*(1 + sigma*(-1 + u) - 4.0*u))*
                                         pow(omega,2.0)) + q2*(q2*(1 - 2.0*u + t2*(6.0 - 10.0*u) + 2.0*t*(-3.0 + 5.0*u) +
                                         sigma*(2.0 + (-4.0 - 8.0*t2 + 8.0*t)*u)) -
                                         4.0*u2*pow(omega,2.0)*(-1 + 2.0*u - 4.0*t2*(-2.0 + u) + 4.0*t*(-2.0 + u) +
                                         2.0*sigma*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))));

                    return C_3 * phi_bar_4 + C_3d1 * phi_double_bar_4;
                }

                inline
                double I1_A_psi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*(-(m_B*sigmabar*(2.0*(-m_v2 + m_B2*sigmabar2)*t*tbar*(-2.0 + u) +
                                       q2*(1 - 2.0*u + t2*(4.0 - 10.0*u) + 2.0*t*(-2.0 + 5.0*u) +
                                       2.0*sigma*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)))) +
                                       4.0*omega*u*(2.0*(m_B2 - m_v2)*t*tbar*(-1 + u) +
                                       m_B2*sigma3*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       q2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) -
                                       2.0*m_B2*sigma2*(-1 + 2.0*u + t*(3.0 - 7.0*u) + t2*(-3.0 + 7.0*u)) +
                                       sigma*(2.0*q2*(-tbar + 2.0*u - 5.0*t*u + t2*(-1 + 5.0*u)) +
                                       m_B2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)) +
                                       m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)))) -
                                       4.0*m_B*u2*(-1 + 2.0*sigma)*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(omega,2.0))*
                                       pow(sigmabar,-3.0))/2.0;

                    return C_1 * psi_bar_4;
                }

                inline
                double I2_A_psi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

                    const double C_2 = ((-m_v2 + m_B2*sigmabar2 + q2*(-1 + 2.0*sigma))*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_B,-2.0)*
                                       (2.0*m_B2*omega*sigmabar2*u + 2.0*omega*(-m_v2 + q2)*u - m_B*sigmabar*(q2 + 4.0*u2*pow(omega,2.0)))*
                                       pow(sigmabar,-4.0))/2.0;

                    return C_2 * psi_bar_4;
                }

                inline
                double I1_A_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to psi_double_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double t2        = pow(t, 2);

                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_1 = m_v*(m_B*sigmabar*t*tbar*(-2.0 + u) + 2.0*omega*u*
                                       (-(2.0*t*tbar*(-1 + u)) + sigma*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))))*pow(m_B,-2.0)*
                                       pow(sigmabar,-3.0);

                    return C_1 * psi_double_bar_4;
                }

                inline
                double I2_A_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_double_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigma2    = pow(sigma, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_2 = m_v*pow(m_B,-2.0)*(-(m_B*sigmabar*(m_B2*sigma2*t*tbar*(-2.0 + u) +
                                       t*tbar*(m_B2*(-2.0 + u) - m_v2*(-2.0 + u) + q2*(-2.0 + 3.0*u)) +
                                       sigma*(-(2.0*m_B2*t*tbar*(-2.0 + u)) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)))) +
                                       2.0*omega*u*(m_B2*sigmabar2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u) +
                                       sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       m_v2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u) +
                                       sigma*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       q2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u) +
                                       sigma*(-3.0 + 6.0*u + t*(2.0 - 14.0*u) + 2.0*t2*(-1 + 7.0*u)))) -
                                       4.0*m_B*u2*sigma*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(omega,2.0))*pow(sigmabar,-4.0);

                    return C_2 * psi_double_bar_4;
                }

                inline
                double I3_A_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_B,-2.0)*
                                       (-(m_B*q2*(-m_v2 + m_B2*sigmabar2 + q2*(-1 + 2.0*sigma))*sigmabar) +
                                       2.0*omega*(m_B4*sigmabar4 + 2.0*m_B2*sigmabar2*(-m_v2 + q2*sigma) +
                                       (m_v2 - q2)*(m_v2 + q2 - 2.0*q2*sigma))*u -
                                       4.0*m_B*u2*(-m_v2 + m_B2*sigmabar2 + q2*(-1 + 2.0*sigma))*sigmabar*pow(omega,2.0))*pow(sigmabar,-5.0));

                    return C_3 * psi_double_bar_4;
                }

                inline
                double I3d1_A_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);
                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*(-(q2*(-m_v2 + m_B2*sigmabar2 + q2*(-1 + 2.0*sigma))*sigmabar) +
                                       2.0*omega*(m_B4*sigmabar4 + 2.0*m_B2*sigmabar2*(-m_v2 + q2*sigma) +
                                       (m_v2 - q2)*(m_v2 + q2 - 2.0*q2*sigma))*u*pow(m_B,-1) -
                                       4.0*u2*(-m_v2 + m_B2*sigmabar2 + q2*(-1 + 2.0*sigma))*sigmabar*pow(omega,2.0))*pow(sigmabar,-5.0));

                    const double C_3d1 = -(2.0*m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_B,-2.0)*
                                         (m_B4*omega*sigmabar4*u + 2.0*m_B2*omega*sigmabar2*(-(3.0*m_v2) + q2 + 2.0*q2*sigma)*u +
                                         omega*(m_v2 - q2)*(5.0*m_v2 + q2*(3.0 - 8.0*sigma))*u -
                                         m_B3*sigmabar3*(q2 + 4.0*u2*pow(omega,2.0)) +
                                         m_B*(2.0*m_v2 + q2 - 3.0*q2*sigma)*sigmabar*(q2 + 4.0*u2*pow(omega,2.0)))*pow(sigmabar,-6.0));

                    return C_3 * psi_bar_4 + C_3d1 * psi_double_bar_4;
                }

                inline
                double I1_A_chi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

                    const double C_1 = pow(m_B,-2.0)*(m_B*(t*tbar*(-(3.0*m_v2*u) + 2.0*m_B*m_v*sigmabar*(-2.0 + u) +
                                       m_B2*sigmabar2*(-6.0 + 7.0*u)) +
                                       q2*(1 - 2.0*u + t2*(2.0 - 5.0*u) + t*(-2.0 + 5.0*u) + sigma*(-1 + 2.0*(1 + t2 - t)*u))) +
                                       omega*u*(-(4.0*m_v2*t*tbar*u) + 2.0*q2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       2.0*m_B2*sigmabar*(1 - 2.0*u + t2*(4.0 - 10.0*u) + 2.0*t*(-2.0 + 5.0*u) +
                                       3.0*sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) -
                                       m_B*m_v*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u +
                                       4.0*sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)))) +
                                       4.0*u2*(2.0*m_v*t*tbar - m_B*sigmabar*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u)))*
                                       pow(omega,2.0))*pow(sigmabar,-2.0);

                    return C_1 * chi_bar_4;
                }

                inline
                double I2_A_chi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

                    const double C_2 = -(pow(m_B,-2.0)*pow(sigmabar,-3.0)*(2.0*omega*u*
                                       (-(16.0*m_B4*t*tbar*(-1 + u)) + m_B3*m_v*
                                       (-1 + 2.0*u + 4.0*t2*(-4.0 + 3.0*u) - 4.0*t*(-4.0 + 3.0*u)) +
                                       4.0*m_B4*sigma4*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       2.0*m_v2*(m_v2 - q2)*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)) +
                                       2.0*m_B2*(8.0*q2*t*(-tbar + u - t*u) + m_v2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       m_B*m_v*(m_v2*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u) +
                                       q2*(-1 + 2.0*u + 4.0*t2*(2.0 + u) - 4.0*t*(2.0 + u))) +
                                       m_B2*sigma2*(-(2.0*m_v2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)) +
                                       4.0*q2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)) +
                                       m_B*m_v*(-9.0 + t*(32 - 60*u) + 18.0*u + 4.0*t2*(-8.0 + 15.0*u)) +
                                       12.0*m_B2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u))) -
                                       4.0*m_B3*sigma3*(m_B*(-3.0 + t*(10.0 - 22*u) + 6.0*u + 2.0*t2*(-5.0 + 11.0*u)) +
                                       m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       2.0*m_B*sigma*(m_B2*m_v*(3.0 - 6.0*u - 4.0*t2*(-5.0 + 6.0*u) + 4.0*t*(-5.0 + 6.0*u)) +
                                       m_B3*(2.0 - 4.0*u - 4.0*t2*(-7.0 + 9.0*u) + 4.0*t*(-7.0 + 9.0*u)) +
                                       2.0*m_B*(m_v2*t*tbar + q2*(1 - 2.0*u + 2.0*t2*(-3.0 + u) - 2.0*t*(-3.0 + u))) +
                                       2.0*(m_v*q2*(1 + t - 2.0*u + 6.0*t*u - t2*(1 + 6.0*u)) +
                                       m_v3*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))))) +
                                       4.0*u2*(-(4.0*m_v*(m_v2 - q2)*t*tbar) + 8.0*m_B3*t*tbar*(-1 + u) +
                                       m_B2*m_v*(1 - 8.0*t2*(-1 + u) + 8.0*t*(-1 + u) - 2.0*u) +
                                       4.0*m_B3*sigma3*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       4.0*m_B2*sigma2*(m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       2.0*m_B*(-1 + 2.0*u + t*(3.0 - 7.0*u) + t2*(-3.0 + 7.0*u))) -
                                       2.0*m_B*(4.0*q2*t*(-tbar + u - t*u) + m_v2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))) +
                                       m_B*sigma*(m_B*m_v*(5.0 + 16.0*t2 - 16.0*t)*(-1 + 2.0*u) + 8.0*q2*t*(-tbar + u - t*u) +
                                       4.0*m_B2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)) +
                                       2.0*m_v2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))))*pow(omega,2.0) +
                                       m_B*(4.0*m_v4*t*tbar + 8.0*m_B4*t*tbar*(-1 + u) + 8.0*m_B4*sigma4*t*tbar*(-1 + u) +
                                       4.0*m_B3*m_v*t*tbar*(-2.0 + u) - 4.0*m_B2*t*tbar*(2.0*q2*(-1 + u) + m_v2*(-1 + 2.0*u)) +
                                       2.0*m_v2*q2*(1 - 2.0*u + t2*(2.0 - 8.0*u) + t*(-2.0 + 8.0*u)) +
                                       m_B*m_v*(-(4.0*m_v2*t*tbar*(-2.0 + u)) + q2*(1 - 2.0*u - 4.0*t2*(1 + u) + 4.0*t*(1 + u))) +
                                       4.0*m_B2*sigma3*(q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B*t*tbar*(8.0*m_B*(-1 + u) + m_v*(-2.0 + u))) -
                                       4.0*m_B*sigma2*(-(12.0*m_B3*t*tbar*(-1 + u)) - 3.0*m_B2*m_v*t*tbar*(-2.0 + u) +
                                       m_v*q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       m_B*(m_v2*t*(-tbar + 2.0*u - 2.0*t*u) + 2.0*q2*(-1 + t2 - t + 2.0*u + 3.0*t2*u - 3.0*t*u))) +
                                       sigma*(-32*m_B4*t*tbar*(-1 + u) - 12.0*m_B3*m_v*t*tbar*(-2.0 + u) +
                                       4.0*m_B2*(q2*(-1 + 4.0*t2 - 4.0*t + 2.0*u) + 2.0*m_v2*t*(-tbar + 2.0*u - 2.0*t*u)) +
                                       m_B*(4.0*m_v3*t*tbar*(-2.0 + u) +
                                       m_v*q2*(-5.0 + 10.0*u + 4.0*t2*(1 + 7.0*u) - 4.0*t*(1 + 7.0*u))) +
                                       2.0*m_v2*q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/2.0;

                    return C_2 * chi_bar_4;
                }

                inline
                double I1_A_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to chi_double_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double sigma2    = pow(sigma, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*(-(m_B*sigmabar*t*tbar*(-(4.0*m_B*sigmabar) + m_v*(-2.0 + u))) -
                                       2.0*omega*u*(-(2.0*t*tbar*(m_v*(-1 + u) + 2.0*m_B*u)) +
                                       2.0*m_B*sigma2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)) +
                                       sigma*(m_B*(2.0 - 4.0*u + t2*(4.0 - 20.0*u) + 4.0*t*(-1 + 5.0*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)))) +
                                       16.0*u2*sigmabar*t*tbar*(-1 + u)*pow(omega,2.0))*pow(sigmabar,-3.0));

                    return C_1 * chi_double_bar_4;
                }

                inline
                double I2_A_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_double_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_2 = pow(m_B,-2.0)*pow(sigmabar,-4.0)*(2.0*omega*u*
                                       (4.0*m_B3*t*tbar*(-2.0 + u) - 4.0*m_B*q2*t*tbar*(-2.0 + u) +
                                       2.0*m_B3*sigma4*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_B2*m_v*(-1 + 2.0*t2*(-1 + u) - 2.0*t*(-1 + u) + 2.0*u) -
                                       4.0*m_B*m_v2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_v3*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v*q2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_B2*sigma3*(m_B*(6.0 - 12.0*u + t2*(4.0 - 20.0*u) + 4.0*t*(-1 + 5.0*u)) +
                                       3.0*m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       sigma*(m_B3*(2.0 - 4.0*u + 4.0*t2*(-5.0 + u) - 4.0*t*(-5.0 + u)) +
                                       m_B2*m_v*(5.0 + t2*(10.0 - 22*u) - 10.0*u + 2.0*t*(-5.0 + 11.0*u)) +
                                       m_v3*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_v*q2*(-3.0 + 6.0*u + t*(2.0 - 14.0*u) + 2.0*t2*(-1 + 7.0*u)) +
                                       2.0*m_B*(m_v2*(-1 + 2.0*t2 - 2.0*t + 2.0*u) +
                                       q2*(-3.0 + 6.0*u + 2.0*t2*(1 + 7.0*u) - 2.0*t*(1 + 7.0*u)))) +
                                       m_B*sigma2*(6.0*m_B2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)) +
                                       m_B*m_v*(-7.0 + t*(14.0 - 38*u) + 14.0*u + 2.0*t2*(-7.0 + 19.0*u)) +
                                       2.0*(q2*(3.0 - 6.0*u + t2*(2.0 - 16.0*u) + 2.0*t*(-1 + 8.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))))) +
                                       4.0*u2*sigmabar*(-(4.0*(3.0*m_B2 + m_B*m_v + m_v2 - 3.0*q2)*t*tbar*(-1 + u)) +
                                       2.0*m_B2*sigma2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)) +
                                       m_B*sigma*(m_v*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u) +
                                       m_B*(2.0 - 4.0*u - 8.0*t2*(-3.0 + 4.0*u) + 8.0*t*(-3.0 + 4.0*u))))*pow(omega,2.0) -
                                       m_B*sigmabar*(4.0*m_B3*sigma3*t*(-tbar + 2.0*u - 2.0*t*u) -
                                       t*tbar*(3.0*m_B2*m_v*(-2.0 + u) + m_v3*(-2.0 + u) + m_v*q2*(2.0 - 3.0*u) + m_B3*(-4.0 + 8.0*u) -
                                       4.0*m_B*(m_v2*(1 + u) + q2*(-1 + 2.0*u))) +
                                       sigma*(6.0*m_B2*m_v*t*tbar*(-2.0 + u) + 12.0*m_B3*t*tbar*(-1 + 2.0*u) +
                                       m_v*q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       m_B*(-(4.0*m_v2*t*tbar*(1 + u)) + q2*(-2.0 + t*(4.0 - 24*u) + 4.0*u + 4.0*t2*(-1 + 6.0*u)))) +
                                       m_B*sigma2*(-(3.0*m_B*t*tbar*(m_v*(-2.0 + u) + m_B*(-4.0 + 8.0*u))) -
                                       2.0*q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0)))));

                    return C_2 * chi_double_bar_4;
                }

                inline
                double I3_A_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4), m_B5 = pow(m_B, 5);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5), sigma6 = pow(sigma, 6);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(pow(m_B,-2.0)*pow(sigmabar,-5.0)*(-(m_B*sigmabar*
                                       (-(8.0*m_B5*t*tbar*(-1 + u)) + 8.0*m_B5*sigma5*t*tbar*(-1 + u) - 4.0*m_B4*m_v*t*tbar*(-2.0 + u) -
                                       m_v*q2*(m_v2 + q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B2*m_v*(-(4.0*m_v2*t*tbar*(-2.0 + u)) + q2*(1 - 8.0*t2 + 8.0*t - 2.0*u)) +
                                       4.0*m_B3*t*tbar*(4.0*q2*(-1 + u) + m_v2*(-2.0 + 3.0*u)) +
                                       4.0*m_B3*sigma4*(q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B*t*tbar*(10.0*m_B*(-1 + u) + m_v*(-2.0 + u))) +
                                       m_B*sigma2*(-80*m_B4*t*tbar*(-1 + u) - 24*m_B3*m_v*t*tbar*(-2.0 + u) +
                                       4.0*q2*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u - 2.0*m_v2*t*tbar*u) +
                                       12.0*m_B2*(q2*(-1 + 4.0*t2 - 4.0*t + 2.0*u) + m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                       m_B*(4.0*m_v3*t*tbar*(-2.0 + u) +
                                       m_v*q2*(-9.0 + 18.0*u + 8.0*t2*(1 + 6.0*u) - 8.0*t*(1 + 6.0*u)))) -
                                       4.0*m_B2*sigma3*(-(20.0*m_B3*t*tbar*(-1 + u)) - 4.0*m_B2*m_v*t*tbar*(-2.0 + u) +
                                       m_v*q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       m_B*(m_v2*t*tbar*(-2.0 + 3.0*u) + q2*(-3.0 + 6.0*u - 4.0*t*(1 + 2.0*u) + t2*(4.0 + 8.0*u))))
                                       + 4.0*m_B*(-(m_v4*t*tbar*u) + m_v2*q2*(1 - 2.0*u + t2*(2.0 - 5.0*u) + t*(-2.0 + 5.0*u)) -
                                       2.0*t*tbar*(-1 + u)*pow(q2,2.0)) -
                                       2.0*sigma*(-(20.0*m_B5*t*tbar*(-1 + u)) - 8.0*m_B4*m_v*t*tbar*(-2.0 + u) +
                                       m_B2*(4.0*m_v3*t*tbar*(-2.0 + u) +
                                       m_v*q2*(-3.0 + 6.0*u + 4.0*t2*(2.0 + 3.0*u) - 4.0*t*(2.0 + 3.0*u))) -
                                       2.0*m_B3*(-(3.0*m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                       q2*(1 - 2.0*u + 4.0*t2*(-3.0 + 2.0*u) + t*(12.0 - 8.0*u))) +
                                       m_v*(1 + (-2.0 - 4.0*t2 + 4.0*t)*u)*pow(q2,2.0) +
                                       2.0*m_B*(-(m_v4*t*tbar*u) + m_v2*q2*(1 - 2.0*u + t2*(2.0 - 3.0*u) + t*(-2.0 + 3.0*u)) +
                                       (1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u))*pow(q2,2.0))))) +
                                       2.0*omega*u*(-(16.0*m_B5*t*tbar*(-1 + u)) -
                                       2.0*m_B2*m_v*(-(4.0*q2*t*tbar*(-2.0 + u)) +
                                       m_v2*(-1 + 8.0*t2*(-1 + u) - 8.0*t*(-1 + u) + 2.0*u)) +
                                       m_B4*m_v*(-1 + 2.0*u + 4.0*t2*(-4.0 + 3.0*u) - 4.0*t*(-4.0 + 3.0*u)) +
                                       4.0*m_B5*sigma6*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       4.0*m_B4*sigma5*(m_B*(-5.0 + t*(14.0 - 34*u) + 10.0*u + 2.0*t2*(-7.0 + 17.0*u)) +
                                       m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) -
                                       4.0*m_B3*(-(8.0*q2*t*tbar*(-1 + u)) + m_v2*(-1 + 2.0*u + t*(4.0 - 9.0*u) + t2*(-4.0 + 9.0*u))) +
                                       2.0*m_B*sigma2*(m_B3*m_v*(-11.0 + t*(64 - 84*u) + 22*u + t2*(-64 + 84*u)) +
                                       10.0*m_B4*(-1 + 2.0*u + 2.0*t2*(-5.0 + 7.0*u) - 2.0*t*(-5.0 + 7.0*u)) +
                                       2.0*q2*(q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) -
                                       6.0*m_B2*(-(12.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-2.0 + 4.0*u + t*(6.0 - 15.0*u) + 3.0*t2*(-2.0 + 5.0*u))) +
                                       m_B*(-(m_v3*(5.0 + 16.0*t2 - 16.0*t)*(-1 + 2.0*u)) +
                                       2.0*m_v*q2*(-1 + 2.0*u - 2.0*t*(4.0 + 3.0*u) + t2*(8.0 + 6.0*u)))) -
                                       2.0*m_B2*sigma3*(2.0*m_B2*m_v*(-7.0 + t*(28 - 48*u) + 14.0*u + 4.0*t2*(-7.0 + 12.0*u)) +
                                       20.0*m_B3*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)) +
                                       m_B*(40*q2*t*tbar*(-1 + u) + m_v2*(8.0 + t2*(20.0 - 54*u) - 16.0*u + t*(-20.0 + 54*u))) -
                                       m_v*(2.0*m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       q2*(1 - 2.0*u - 4.0*t2*(1 + 2.0*u) + t*(4.0 + 8.0*u)))) +
                                       m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*(m_v4 - pow(q2,2.0)) +
                                       4.0*m_B*(m_v4*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                       m_v2*q2*(1 - 2.0*u + t2*(4.0 - 7.0*u) + t*(-4.0 + 7.0*u)) - 4.0*t*tbar*(-1 + u)*pow(q2,2.0)) -
                                       2.0*sigma*(4.0*m_B4*m_v*(-1 + 9.0*t2*(-1 + u) - 9.0*t*(-1 + u) + 2.0*u) +
                                       m_v*(m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_B5*(-2.0 + t*(44 - 52*u) + 4.0*u + t2*(-44 + 52*u)) +
                                       m_B3*(56*q2*t*tbar*(-1 + u) + m_v2*(8.0 + t2*(28 - 66*u) - 16.0*u + t*(-28 + 66*u))) -
                                       m_B2*m_v*(q2*(1 - 20.0*t2 + 20.0*t - 2.0*u) +
                                       4.0*m_v2*(-1 + 2.0*u + t*(5.0 - 7.0*u) + t2*(-5.0 + 7.0*u))) +
                                       2.0*m_B*(3.0*m_v2*q2*t*tbar*(-2.0 + u) + m_v4*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                       (1 - 2.0*u + 2.0*t2*(-3.0 + u) - 2.0*t*(-3.0 + u))*pow(q2,2.0))) +
                                       m_B3*sigma4*(m_B*m_v*(-17.0 + 34*u + 12.0*t2*(-4.0 + 9.0*u) - 12.0*t*(-4.0 + 9.0*u)) -
                                       4.0*(-(4.0*q2*t*tbar*(-1 + u)) + m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       40*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))) -
                                       4.0*u2*sigmabar*pow(omega,2.0)*(-(8.0*m_B4*t*tbar*(-1 + u)) + 8.0*(m_v2 - q2)*q2*t*tbar*(-1 + u) +
                                       m_B3*m_v*(-1 + 8.0*t2*(-1 + u) - 8.0*t*(-1 + u) + 2.0*u) -
                                       m_B*m_v*(m_v2 + q2)*(-1 + 8.0*t2*(-1 + u) - 8.0*t*(-1 + u) + 2.0*u) +
                                       4.0*m_B4*sigma4*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       4.0*m_B3*sigma3*(m_B*(-3.0 + 6.0*u + t*(8.0 - 20.0*u) + 4.0*t2*(-2.0 + 5.0*u)) +
                                       m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) -
                                       4.0*m_B2*(-(4.0*q2*t*tbar*(-1 + u)) + m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) -
                                       2.0*m_B*sigma*(m_v*q2*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u) +
                                       4.0*m_v3*t*(-tbar + u - t*u) + 2.0*m_B3*
                                       (-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_B2*m_v*(-3.0 + 6.0*u + 4.0*t2*(-3.0 + 5.0*u) - 4.0*t*(-3.0 + 5.0*u)) -
                                       2.0*m_B*(q2*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_v2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))) +
                                       m_B2*sigma2*(m_B*m_v*(-9.0 + 18.0*u + 8.0*t2*(-3.0 + 7.0*u) - 8.0*t*(-3.0 + 7.0*u)) +
                                       12.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0) -
                                       4.0*(-(2.0*m_v2*t*tbar*(-1 + u)) + q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))))));

                    return C_3 * chi_double_bar_4;
                }

                inline
                double I3d1_A_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4), m_B5 = pow(m_B, 5);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5), sigma6 = pow(sigma, 6);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);
                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(pow(sigmabar,-5.0)*(-(sigmabar*(-(8.0*m_B5*t*tbar*(-1 + u)) + 8.0*m_B5*sigma5*t*tbar*(-1 + u) -
                                       4.0*m_B4*m_v*t*tbar*(-2.0 + u) - m_v*q2*(m_v2 + q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B2*m_v*(-(4.0*m_v2*t*tbar*(-2.0 + u)) + q2*(1 - 8.0*t2 + 8.0*t - 2.0*u)) +
                                       4.0*m_B3*t*tbar*(4.0*q2*(-1 + u) + m_v2*(-2.0 + 3.0*u)) +
                                       4.0*m_B3*sigma4*(q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B*t*tbar*(10.0*m_B*(-1 + u) + m_v*(-2.0 + u))) +
                                       m_B*sigma2*(-80*m_B4*t*tbar*(-1 + u) - 24*m_B3*m_v*t*tbar*(-2.0 + u) +
                                       4.0*q2*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u - 2.0*m_v2*t*tbar*u) +
                                       12.0*m_B2*(q2*(-1 + 4.0*t2 - 4.0*t + 2.0*u) + m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                       m_B*(4.0*m_v3*t*tbar*(-2.0 + u) +
                                       m_v*q2*(-9.0 + 18.0*u + 8.0*t2*(1 + 6.0*u) - 8.0*t*(1 + 6.0*u)))) -
                                       4.0*m_B2*sigma3*(-(20.0*m_B3*t*tbar*(-1 + u)) - 4.0*m_B2*m_v*t*tbar*(-2.0 + u) +
                                       m_v*q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       m_B*(m_v2*t*tbar*(-2.0 + 3.0*u) + q2*(-3.0 + 6.0*u - 4.0*t*(1 + 2.0*u) + t2*(4.0 + 8.0*u))))
                                       + 4.0*m_B*(-(m_v4*t*tbar*u) + m_v2*q2*(1 - 2.0*u + t2*(2.0 - 5.0*u) + t*(-2.0 + 5.0*u)) -
                                       2.0*t*tbar*(-1 + u)*pow(q2,2.0)) -
                                       2.0*sigma*(-(20.0*m_B5*t*tbar*(-1 + u)) - 8.0*m_B4*m_v*t*tbar*(-2.0 + u) +
                                       m_B2*(4.0*m_v3*t*tbar*(-2.0 + u) +
                                       m_v*q2*(-3.0 + 6.0*u + 4.0*t2*(2.0 + 3.0*u) - 4.0*t*(2.0 + 3.0*u))) -
                                       2.0*m_B3*(-(3.0*m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                       q2*(1 - 2.0*u + 4.0*t2*(-3.0 + 2.0*u) + t*(12.0 - 8.0*u))) +
                                       m_v*(1 + (-2.0 - 4.0*t2 + 4.0*t)*u)*pow(q2,2.0) +
                                       2.0*m_B*(-(m_v4*t*tbar*u) + m_v2*q2*(1 - 2.0*u + t2*(2.0 - 3.0*u) + t*(-2.0 + 3.0*u)) +
                                       (1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u))*pow(q2,2.0))))) +
                                       2.0*omega*u*pow(m_B,-1)*(-(16.0*m_B5*t*tbar*(-1 + u)) -
                                       2.0*m_B2*m_v*(-(4.0*q2*t*tbar*(-2.0 + u)) +
                                       m_v2*(-1 + 8.0*t2*(-1 + u) - 8.0*t*(-1 + u) + 2.0*u)) +
                                       m_B4*m_v*(-1 + 2.0*u + 4.0*t2*(-4.0 + 3.0*u) - 4.0*t*(-4.0 + 3.0*u)) +
                                       4.0*m_B5*sigma6*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       4.0*m_B4*sigma5*(m_B*(-5.0 + t*(14.0 - 34*u) + 10.0*u + 2.0*t2*(-7.0 + 17.0*u)) +
                                       m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) -
                                       4.0*m_B3*(-(8.0*q2*t*tbar*(-1 + u)) + m_v2*(-1 + 2.0*u + t*(4.0 - 9.0*u) + t2*(-4.0 + 9.0*u))) +
                                       2.0*m_B*sigma2*(m_B3*m_v*(-11.0 + t*(64 - 84*u) + 22*u + t2*(-64 + 84*u)) +
                                       10.0*m_B4*(-1 + 2.0*u + 2.0*t2*(-5.0 + 7.0*u) - 2.0*t*(-5.0 + 7.0*u)) +
                                       2.0*q2*(q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) -
                                       6.0*m_B2*(-(12.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-2.0 + 4.0*u + t*(6.0 - 15.0*u) + 3.0*t2*(-2.0 + 5.0*u))) +
                                       m_B*(-(m_v3*(5.0 + 16.0*t2 - 16.0*t)*(-1 + 2.0*u)) +
                                       2.0*m_v*q2*(-1 + 2.0*u - 2.0*t*(4.0 + 3.0*u) + t2*(8.0 + 6.0*u)))) -
                                       2.0*m_B2*sigma3*(2.0*m_B2*m_v*(-7.0 + t*(28 - 48*u) + 14.0*u + 4.0*t2*(-7.0 + 12.0*u)) +
                                       20.0*m_B3*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)) +
                                       m_B*(40*q2*t*tbar*(-1 + u) + m_v2*(8.0 + t2*(20.0 - 54*u) - 16.0*u + t*(-20.0 + 54*u))) -
                                       m_v*(2.0*m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       q2*(1 - 2.0*u - 4.0*t2*(1 + 2.0*u) + t*(4.0 + 8.0*u)))) +
                                       m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*(m_v4 - pow(q2,2.0)) +
                                       4.0*m_B*(m_v4*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                       m_v2*q2*(1 - 2.0*u + t2*(4.0 - 7.0*u) + t*(-4.0 + 7.0*u)) - 4.0*t*tbar*(-1 + u)*pow(q2,2.0)) -
                                       2.0*sigma*(4.0*m_B4*m_v*(-1 + 9.0*t2*(-1 + u) - 9.0*t*(-1 + u) + 2.0*u) +
                                       m_v*(m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_B5*(-2.0 + t*(44 - 52*u) + 4.0*u + t2*(-44 + 52*u)) +
                                       m_B3*(56*q2*t*tbar*(-1 + u) + m_v2*(8.0 + t2*(28 - 66*u) - 16.0*u + t*(-28 + 66*u))) -
                                       m_B2*m_v*(q2*(1 - 20.0*t2 + 20.0*t - 2.0*u) +
                                       4.0*m_v2*(-1 + 2.0*u + t*(5.0 - 7.0*u) + t2*(-5.0 + 7.0*u))) +
                                       2.0*m_B*(3.0*m_v2*q2*t*tbar*(-2.0 + u) + m_v4*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                       (1 - 2.0*u + 2.0*t2*(-3.0 + u) - 2.0*t*(-3.0 + u))*pow(q2,2.0))) +
                                       m_B3*sigma4*(m_B*m_v*(-17.0 + 34*u + 12.0*t2*(-4.0 + 9.0*u) - 12.0*t*(-4.0 + 9.0*u)) -
                                       4.0*(-(4.0*q2*t*tbar*(-1 + u)) + m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       40*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))) -
                                       4.0*u2*sigmabar*pow(m_B,-1)*pow(omega,2.0)*
                                       (-(8.0*m_B4*t*tbar*(-1 + u)) + 8.0*(m_v2 - q2)*q2*t*tbar*(-1 + u) +
                                       m_B3*m_v*(-1 + 8.0*t2*(-1 + u) - 8.0*t*(-1 + u) + 2.0*u) -
                                       m_B*m_v*(m_v2 + q2)*(-1 + 8.0*t2*(-1 + u) - 8.0*t*(-1 + u) + 2.0*u) +
                                       4.0*m_B4*sigma4*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       4.0*m_B3*sigma3*(m_B*(-3.0 + 6.0*u + t*(8.0 - 20.0*u) + 4.0*t2*(-2.0 + 5.0*u)) +
                                       m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) -
                                       4.0*m_B2*(-(4.0*q2*t*tbar*(-1 + u)) + m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) -
                                       2.0*m_B*sigma*(m_v*q2*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u) +
                                       4.0*m_v3*t*(-tbar + u - t*u) + 2.0*m_B3*
                                       (-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_B2*m_v*(-3.0 + 6.0*u + 4.0*t2*(-3.0 + 5.0*u) - 4.0*t*(-3.0 + 5.0*u)) -
                                       2.0*m_B*(q2*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_v2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))) +
                                       m_B2*sigma2*(m_B*m_v*(-9.0 + 18.0*u + 8.0*t2*(-3.0 + 7.0*u) - 8.0*t*(-3.0 + 7.0*u)) +
                                       12.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0) -
                                       4.0*(-(2.0*m_v2*t*tbar*(-1 + u)) + q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))))));

                    const double C_3d1 = 8.0*m_B4*t*tbar*(-1 + u) + 8.0*m_B3*omega*u*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                         2.0*m_B*(-(16.0*omega*q2*(-3.0 + sigma)*t*tbar*(-1 + u)*u) - 4.0*m_v3*t*tbar*(-2.0 + u) +
                                         4.0*m_v2*omega*u*(-((1 + 6.0*t2 - 6.0*t)*(-1 + 2.0*u)) +
                                         sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                         m_v*(q2*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u) +
                                         2.0*sigma*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) +
                                         4.0*u2*(-1 + 2.0*u + 4.0*t2*(1 + u) - 4.0*t*(1 + u) +
                                         2.0*sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)))*pow(omega,2.0)))*pow(sigmabar,-3.0)
                                         + 2.0*pow(m_B,-1)*(-(16.0*m_v4*omega*u*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) +
                                         4.0*m_v2*omega*q2*u*(-5.0 + t*(14.0 - 32*u) + 10.0*u + 2.0*t2*(-7.0 + 16.0*u) +
                                         3.0*sigma*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) +
                                         m_v*q2*(q2*(-1 + 3.0*sigma)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                         4.0*u2*(1 - 12.0*t2*(-1 + u) + 12.0*t*(-1 + u) - 2.0*u +
                                         3.0*sigma*(-1 + 4.0*t2*(-1 + u) - 4.0*t*(-1 + u) + 2.0*u))*pow(omega,2.0)) +
                                         2.0*m_v3*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u +
                                         4.0*u2*(1 - 6.0*t2*sigmabar*(-1 + u) - 2.0*u + 6.0*t*(-sigmabar + u - sigma*u))*pow(omega,2.0)) -
                                         4.0*omega*u*(-1 + 2.0*u + 2.0*t2*(-7.0 + 9.0*u) - 2.0*t*(-7.0 + 9.0*u) +
                                         sigma*(-3.0 + 6.0*u + 6.0*t2*(1 + u) - 6.0*t*(1 + u)))*pow(q2,2.0))*pow(sigmabar,-5.0) -
                                         2.0*omega*(m_v2 - q2)*u*(-64*omega*q2*sigmabar*t*tbar*(-1 + u)*u +
                                         5.0*m_v3*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                         m_v*q2*(-3.0 + 8.0*sigma)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))*pow(m_B,-2.0)*pow(sigmabar,-6.0) -
                                         2.0*m_B2*pow(sigmabar,-2.0)*(-(2.0*m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                         m_v*omega*u*(3.0 - 6.0*u - 4.0*t2*(2.0 + 3.0*u) + 4.0*t*(2.0 + 3.0*u)) +
                                         8.0*u2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(omega,2.0) + 2.0*q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)) -
                                         4.0*pow(sigmabar,-4.0)*(3.0*m_v4*t*tbar*u +
                                         m_v3*omega*u*(1 - 2.0*u - 4.0*t2*(-5.0 + 3.0*u) + 4.0*t*(-5.0 + 3.0*u) +
                                         4.0*sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                         m_v*omega*q2*u*(1 - 20.0*t2*(-1 + u) + 20.0*t*(-1 + u) - 2.0*u -
                                         2.0*sigma*(-1 + 2.0*u - 4.0*t*(1 + 2.0*u) + t2*(4.0 + 8.0*u))) +
                                         m_v2*(q2*(-3.0 + 6.0*u + t*(6.0 - 17.0*u - 4.0*sigma*u) + t2*(-6.0 + 17.0*u + 4.0*sigma*u)) -
                                         4.0*u2*(3.0 - 6.0*u + 4.0*t2*(1 + sigma*(-1 + u) - 4.0*u) - 4.0*t*(1 + sigma*(-1 + u) - 4.0*u))*
                                         pow(omega,2.0)) + q2*(q2*(1 - 2.0*u + t2*(6.0 - 10.0*u) + 2.0*t*(-3.0 + 5.0*u) +
                                         sigma*(2.0 + (-4.0 - 8.0*t2 + 8.0*t)*u)) -
                                         4.0*u2*pow(omega,2.0)*(-1 + 2.0*u - 4.0*t2*(-2.0 + u) + 4.0*t*(-2.0 + u) +
                                         2.0*sigma*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))));

                    return C_3 * chi_bar_4 + C_3d1 * chi_double_bar_4;
                }

                double integrand_A(const std::array<double, 4> & args, const double & q2) const
                {
                    const double sigma = args[0];
                    const double x     = args[1];
                    const double t     = args[2];
                    const double u     = args[3];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = - f_B() * m_B / (8.0  * M_PI * M_PI * f_P() * (m_B * m_B - m_P() * m_P() - q2));
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_P2 = pow(m_P(), 2);
                    const double M4   = pow(M2, 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

                    const double I1 = I1_A_phi_3(sigma, omega, t, u, q2)             + I1_A_phi_4(sigma, omega, t, u, q2)
                                    + I1_A_phi_bar_3(sigma, omega, t, u, q2)         + I1_A_phi_bar_4(sigma, omega, t, u, q2)
                                    + I1_A_phi_double_bar_3(sigma, omega, t, u, q2)  + I1_A_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I1_A_psi_bar_4(sigma, omega, t, u, q2)         + I1_A_chi_bar_4(sigma, omega, t, u, q2)
                                    + I1_A_psi_double_bar_4(sigma, omega, t, u, q2)  + I1_A_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I2 = I2_A_phi_bar_3(sigma, omega, t, u, q2)         + I2_A_phi_bar_4(sigma, omega, t, u, q2)
                                    + I2_A_phi_double_bar_3(sigma, omega, t, u, q2)  + I2_A_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I2_A_psi_bar_4(sigma, omega, t, u, q2)         + I2_A_chi_bar_4(sigma, omega, t, u, q2)
                                    + I2_A_psi_double_bar_4(sigma, omega, t, u, q2)  + I2_A_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3 = I3_A_phi_double_bar_3(sigma, omega, t, u, q2)  + I3_A_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I3_A_psi_double_bar_4(sigma, omega, t, u, q2)  + I3_A_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result = 0.0;
                    result += - I1;
                    result +=   I2 / M2;
                    result += - I3 / (2.0 * M4);
                    result *=   prefactor * exp;

                    return result;
                }

                double surface_A(const std::array<double, 3> & args, const double & sigma, const double & q2) const
                {
                    const double x     = args[0];
                    const double t     = args[1];
                    const double u     = args[2];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = - f_B() * m_B / (8.0  * M_PI * M_PI * f_P() * (m_B * m_B - m_P() * m_P() - q2));
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_B2 = pow(m_B(), 2);
                    const double m_v2 = pow(m_v(), 2);
                    const double m_P2 = pow(m_P(), 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

                    const double sigmabar = 1.0 - sigma, sigmabar2 = pow(sigmabar, 2);
                    const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
                    const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;

                    const double I2   = I2_A_phi_bar_3(sigma, omega, t, u, q2)           + I2_A_phi_bar_4(sigma, omega, t, u, q2)
                                      + I2_A_phi_double_bar_3(sigma, omega, t, u, q2)    + I2_A_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I2_A_psi_bar_4(sigma, omega, t, u, q2)           + I2_A_chi_bar_4(sigma, omega, t, u, q2)
                                      + I2_A_psi_double_bar_4(sigma, omega, t, u, q2)    + I2_A_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3   = I3_A_phi_double_bar_3(sigma, omega, t, u, q2)    + I3_A_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3_A_psi_double_bar_4(sigma, omega, t, u, q2)    + I3_A_chi_double_bar_4(sigma, omega, t, u, q2);

                    const double I3d1 = I3d1_A_phi_double_bar_3(sigma, omega, t, u, q2)  + I3d1_A_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3d1_A_psi_double_bar_4(sigma, omega, t, u, q2)  + I3d1_A_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result = 0.0;
                    result +=  eta * I2 / m_B2;
                    result += -0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
                    result *=  prefactor * exp;

                    return result;
                }

                double integrand_A_m1(const std::array<double, 4> & args, const double & q2) const
                {
                    const double sigma = args[0];
                    const double x     = args[1];
                    const double t     = args[2];
                    const double u     = args[3];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = - f_B() * m_B / (8.0  * M_PI * M_PI * f_P() * (m_B * m_B - m_P() * m_P() - q2));
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_P2 = pow(m_P(), 2);
                    const double M4   = pow(M2, 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

                    const double I1 = I1_A_phi_3(sigma, omega, t, u, q2)             + I1_A_phi_4(sigma, omega, t, u, q2)
                                    + I1_A_phi_bar_3(sigma, omega, t, u, q2)         + I1_A_phi_bar_4(sigma, omega, t, u, q2)
                                    + I1_A_phi_double_bar_3(sigma, omega, t, u, q2)  + I1_A_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I1_A_psi_bar_4(sigma, omega, t, u, q2)         + I1_A_chi_bar_4(sigma, omega, t, u, q2)
                                    + I1_A_psi_double_bar_4(sigma, omega, t, u, q2)  + I1_A_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I2 = I2_A_phi_bar_3(sigma, omega, t, u, q2)         + I2_A_phi_bar_4(sigma, omega, t, u, q2)
                                    + I2_A_phi_double_bar_3(sigma, omega, t, u, q2)  + I2_A_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I2_A_psi_bar_4(sigma, omega, t, u, q2)         + I2_A_chi_bar_4(sigma, omega, t, u, q2)
                                    + I2_A_psi_double_bar_4(sigma, omega, t, u, q2)  + I2_A_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3 = I3_A_phi_double_bar_3(sigma, omega, t, u, q2)  + I3_A_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I3_A_psi_double_bar_4(sigma, omega, t, u, q2)  + I3_A_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result1 =    0.0;
                           result1 += - I1;
                           result1 +=   I2 / M2;
                           result1 += - I3 / (2.0 * M4);
                           result1 *=   exp * s(sigma, q2);

                    double result2 =    0.0;
                           result2 += - I2;
                           result2 +=   I3 / M2;
                           result2 *=   exp;

                    return prefactor * (result1 + result2);
                }

                double surface_A_m1(const std::array<double, 3> & args, const double & sigma, const double & q2) const
                {
                    const double x     = args[0];
                    const double t     = args[1];
                    const double u     = args[2];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = - f_B() * m_B / (8.0  * M_PI * M_PI * f_P() * (m_B * m_B - m_P() * m_P() - q2));
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_B2 = pow(m_B(), 2);
                    const double m_v2 = pow(m_v(), 2);
                    const double m_P2 = pow(m_P(), 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_P2) / M2());

                    const double sigmabar = 1.0 - sigma, sigmabar2 = pow(sigmabar, 2);
                    const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
                    const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;

                    const double I2   = I2_A_phi_bar_3(sigma, omega, t, u, q2)           + I2_A_phi_bar_4(sigma, omega, t, u, q2)
                                      + I2_A_phi_double_bar_3(sigma, omega, t, u, q2)    + I2_A_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I2_A_psi_bar_4(sigma, omega, t, u, q2)           + I2_A_chi_bar_4(sigma, omega, t, u, q2)
                                      + I2_A_psi_double_bar_4(sigma, omega, t, u, q2)    + I2_A_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3   = I3_A_phi_double_bar_3(sigma, omega, t, u, q2)    + I3_A_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3_A_psi_double_bar_4(sigma, omega, t, u, q2)    + I3_A_chi_double_bar_4(sigma, omega, t, u, q2);

                    const double I3d1 = I3d1_A_phi_double_bar_3(sigma, omega, t, u, q2)  + I3d1_A_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3d1_A_psi_double_bar_4(sigma, omega, t, u, q2)  + I3d1_A_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result1 =   0.0;
                           result1 +=  eta * I2 / m_B2;
                           result1 += -0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
                           result1 *=  exp * s(sigma, q2);

                    double result2 =  0.0;
                           result2 += 0.5 * eta * I3 / m_B2;
                           result2 *= exp;

                    return prefactor * (result1 + result2);
                }

                double A(const double & q2) const
                {
                    const double sigma_0 = this->sigma_0(q2, s0_0_A(), s0_1_A());

                    double value_integral = 0.0;
                    double value_surface  = 0.0;

                    const std::function<double (const std::array<double, 4> &)> integrand = std::bind(&LCSR::integrand_A, this, std::placeholders::_1, q2);
                    const std::function<double (const std::array<double, 3> &)> surface = std::bind(&LCSR::surface_A, this, std::placeholders::_1, sigma_0, q2);

                    value_integral = integrate(integrand, { 0.0, 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0, 1.0 }, cubature::Config());
                    value_surface  = integrate(surface, { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, cubature::Config());

                    return value_integral + value_surface;
                }

                double normalized_first_moment_A(const double & q2) const
                {
                    const double sigma_0 = this->sigma_0(q2, s0_0_A(), s0_1_A());

                    double value_integral = 0.0;
                    double value_surface  = 0.0;

                    const std::function<double (const std::array<double, 4> &)> integrand = std::bind(&LCSR::integrand_A, this, std::placeholders::_1, q2);
                    const std::function<double (const std::array<double, 3> &)> surface = std::bind(&LCSR::surface_A, this, std::placeholders::_1, sigma_0, q2);

                    value_integral = integrate(integrand, { 0.0, 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0, 1.0 }, cubature::Config());
                    value_surface  = integrate(surface, { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, cubature::Config());

                    const double denominator = value_integral + value_surface;

                    double value_integral_m1 = 0.0;
                    double value_surface_m1  = 0.0;

                    const std::function<double (const std::array<double, 4> &)> integrand_m1 = std::bind(&LCSR::integrand_A_m1, this, std::placeholders::_1, q2);
                    const std::function<double (const std::array<double, 3> &)> surface_m1 = std::bind(&LCSR::surface_A_m1, this, std::placeholders::_1, sigma_0, q2);

                    value_integral_m1 = integrate(integrand_m1, { 0.0, 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0, 1.0 }, cubature::Config());
                    value_surface_m1  = integrate(surface_m1, { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, cubature::Config());

                    const double numerator = value_integral_m1 + value_surface_m1;

                    return numerator / denominator;
                }

                /* Diagnostics */

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    /* paramaters */
                    results.add({ this->m_v(),                                                     "m_v(mu) in the MSbar scheme"                   });
                    results.add({ this->m_c(),                                                     "m_c(mu) in the MSbar scheme"                   });
                    results.add({ this->f_P(),                                                     "final state decay constant"                    });

                    /* I phi_3(sigma, omega, t, u, q2) */
                    results.add({ this->I1_A_phi_3(0.02, 0.4, 0.3, 0.7, 2.0),                      "I1_A_phi_3(0.02, 0.4, 0.3, 0.7, 2.0)"          });
                    results.add({ this->I1_A_phi_3(0.01, 0.7, 0.9, 0.1, 3.0),                      "I1_A_phi_3(0.01, 0.7, 0.9, 0.1, 3.0)"          });

                    results.add({ this->I1_A_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),                  "I1_A_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I1_A_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),                  "I1_A_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)"      });
                    results.add({ this->I2_A_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),                  "I2_A_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I2_A_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),                  "I2_A_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)"      });

                    results.add({ this->I1_A_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),           "I1_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I1_A_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),           "I1_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I2_A_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),           "I2_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I2_A_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),           "I2_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3_A_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),           "I3_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I3_A_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),           "I3_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3d1_A_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),         "I3d1_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)" });
                    results.add({ this->I3d1_A_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),         "I3d1_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)" });

                    /* I phi_4(sigma, omega, t, u, q2) */
                    results.add({ this->I1_A_phi_4(0.02, 0.4, 0.3, 0.7, 2.0),                      "I1_A_phi_4(0.02, 0.4, 0.3, 0.7, 2.0)"          });
                    results.add({ this->I1_A_phi_4(0.01, 0.7, 0.9, 0.1, 3.0),                      "I1_A_phi_4(0.01, 0.7, 0.9, 0.1, 3.0)"          });

                    results.add({ this->I1_A_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I1_A_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I1_A_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I1_A_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });
                    results.add({ this->I2_A_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I2_A_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I2_A_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I2_A_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });

                    results.add({ this->I1_A_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I1_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I1_A_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I1_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I2_A_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I2_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I2_A_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I2_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3_A_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I3_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I3_A_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I3_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3d1_A_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),         "I3d1_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)" });
                    results.add({ this->I3d1_A_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),         "I3d1_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)" });

                    /* I psi_4(sigma, omega, t, u, q2) */
                    results.add({ this->I1_A_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I1_A_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I1_A_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I1_A_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });
                    results.add({ this->I2_A_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I2_A_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I2_A_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I2_A_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });

                    results.add({ this->I1_A_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I1_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I1_A_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I1_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I2_A_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I2_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I2_A_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I2_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3_A_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I3_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I3_A_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I3_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3d1_A_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),         "I3d1_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)" });
                    results.add({ this->I3d1_A_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),         "I3d1_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)" });

                    /* I chi_4(sigma, omega, t, u, q2) */
                    results.add({ this->I1_A_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I1_A_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I1_A_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I1_A_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });
                    results.add({ this->I2_A_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I2_A_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I2_A_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I2_A_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });

                    results.add({ this->I1_A_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I1_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I1_A_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I1_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I2_A_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I2_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I2_A_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I2_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3_A_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I3_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I3_A_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I3_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3d1_A_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),         "I3d1_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)" });
                    results.add({ this->I3d1_A_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),         "I3d1_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)" });

                    return results;
                }

                virtual complex<double> H_plus(const double & q2) const
                {
                    const double m_B     = this-> m_B();
                    const double Q_c     = 2.0 / 3.0; //Charm quark charge
                    auto         wc      = model->wilson_coefficients_b_to_s(mu_ren(), "mu", false);
                    auto         C_1_EOS = wc.c1();
                    auto         C_2_EOS = wc.c2();
                    auto         C_1_AK  = C_2_EOS - 1.0 / 6.0 * C_1_EOS;

                    return - (2.0 * C_1_AK) * q2 * Q_c * this->A(q2) / (2.0 * m_B * m_B);
                }

                virtual complex<double> Hhat_plus(const double & q2) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_A(const double & q2) const
                {
                    return this->normalized_first_moment_A(q2);
                }

                static NonlocalFormFactorPtr<nc::PToP> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nc::PToP>(new LCSR(p, o));
                }
        };




        /*
         * Parametrize the entire formfactor, i.e., both leading and all sub-leading powers as described in [GvDV:2020].
         */
        template <typename Process_>
        class GvDV2020 :
            public NonlocalFormFactor<nc::PToP>
        {
            public:
                std::shared_ptr<Model> model;

                //Polynomial expansion parameters
                UsedParameter re_alpha_0_plus;
                UsedParameter im_alpha_0_plus;
                UsedParameter re_alpha_1_plus;
                UsedParameter im_alpha_1_plus;
                UsedParameter re_alpha_2_plus;
                UsedParameter im_alpha_2_plus;

                //Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;

                // final state meson parameters
                UsedParameter m_P;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;

                GvDV2020(const Parameters & p, const Options & o) :
                    re_alpha_0_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_0^plus}@GvDV2020"], *this),
                    im_alpha_0_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_0^plus}@GvDV2020"], *this),
                    re_alpha_1_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_1^plus}@GvDV2020"], *this),
                    im_alpha_1_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_1^plus}@GvDV2020"], *this),
                    re_alpha_2_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_2^plus}@GvDV2020"], *this),
                    im_alpha_2_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_2^plus}@GvDV2020"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_d"], *this),

                    m_P(p["mass::K_d"], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),
                    chiOPE(p["b->sccbar::chiOPE@GvDV2020"], *this)
                {
                }

                ~GvDV2020() = default;

                inline complex<double> phi(const double & q2, const unsigned phiParam[4]) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF                        a    b    c    d
                    // 0(P->P) aka plus          3    3    2    2
                    // perp(P->V) = par(P->V)    3    1    3    0
                    // 0(P->V) aka long          3    1    2    2

                    const double m_P2  = pow(m_P, 2);
                    const double m_B2  = pow(m_B, 2),  m_B4 =  pow(m_B, 4);
                    const double m_D02 = pow(m_D0, 2), m_D04 = pow(m_D0, 4);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nc_utils::z(q2, 4.0 * pow(m_D0, 2), s_0);
                    const double Q2   = this->t_s();
                    const double chi = this->chiOPE();

                    const double a = phiParam[0], b = phiParam[1], c = phiParam[2], d = phiParam[3];

                    const complex<double> Nlambda = 4. * M_PI * pow(m_B2, 0.5*(a-b+c+d) - 1.) * pow( 2.*(4.*m_D02-s_0)/3./chi, 0.5); //(C6)
                    const complex<double> phi1 = -pow(2.*pow((4.*m_D02-Q2)*(4.*m_D02-s_0), 0.5) + 8.*m_D02 - Q2 - s_0, 0.5) /
                                                (2.*pow((4.*m_D02-Q2)*(4.*m_D02-s_0), 0.5) + 8.*m_D02 + Q2*(z-1.) - s_0*(z+1.)); //(C7)
                    const complex<double> phi2 = pow(m_B4*pow(z-1., 4.) - 2.*m_B2*pow(z-1., 2)*(-16*m_D02*z + m_P2*pow(z-1., 2) + s_0*pow(z+1., 2)) +
                                                pow(16*m_D02*z + m_P2*pow(z-1., 2) - s_0*pow(z+1., 2), 2), 0.5);//(C8)
                    const complex<double> phi3 = pow(8*m_D02 + 4*pow(4*m_D04-s_0*m_D02, 0.5) - s_0, 0.5)/(-8*m_D02 - 4*pow(4*m_D04-s_0*m_D02, 0.5) + s_0*(z+1.)); //(C9)
                    const complex<double> phi4 = pow(s_0 * pow(z+1., 2.) - 16.*z*m_D02, -0.5); //(C10)

                    return Nlambda * pow(1.+z, 0.5) * pow(1.-z, a-b+c+d-1.5) * pow(phi1, a) * pow(phi2, 0.5*b) * pow(phi3, c) * pow(phi4, d); //(C5)
                }

                inline complex<double> H_residue_jpsi(const unsigned phiParam[4], const complex<double> & alpha_0, const complex<double> & alpha_1,
                                                      const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto zBP     = eos::nc_utils::z(pow(m_B+m_P, 2), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(m_Jpsi2,         s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(m_psi2S2,        s_p, s_0);

                    return eos::nc_utils::PGvDV2020(z_Jpsi, zBP, alpha_0, alpha_1, alpha_2) / phi(m_Jpsi2, phiParam) *
                            (1 - norm(z_Jpsi)) * (1. - z_Jpsi*std::conj(z_psi2S)) / (z_Jpsi - z_psi2S);
                };

                inline complex<double> H_residue_psi2s(const unsigned phiParam[4], const complex<double> & alpha_0, const complex<double> & alpha_1,
                                                       const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto zBP     = eos::nc_utils::z(pow(m_B+m_P, 2), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(m_Jpsi2,         s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(m_psi2S2,        s_p, s_0);

                    return eos::nc_utils::PGvDV2020(z_psi2S, zBP, alpha_0, alpha_1, alpha_2) / phi(m_psi2S2, phiParam) *
                            (1 - norm(z_psi2S)) * (1. - z_psi2S*std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi);
                };


                virtual complex<double> H_plus(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,              s_p, s_0);
                    const auto zBP     = eos::nc_utils::z(pow(m_B+m_P, 2), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(pow(m_Jpsi, 2),  s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(pow(m_psi2S, 2), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nc_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[4] = {3, 3, 2, 2};

                    return eos::nc_utils::PGvDV2020(z, zBP, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                }


                virtual complex<double> Hhat_plus(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2, s_p, s_0);
                    const auto zBP     = eos::nc_utils::z(pow(m_B+m_P, 2), s_p, s_0);

                    return eos::nc_utils::PGvDV2020(z, zBP, alpha_0, alpha_1, alpha_2);
                }


                virtual complex<double> H_plus_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const unsigned phiParam[4] = {3, 3, 2, 2};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_plus_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const unsigned phiParam[4] = {3, 3, 2, 2};

                    return H_residue_psi2s(phiParam, alpha_0, alpha_1, alpha_2);
                };


////// TODOs
                virtual complex<double> normalized_moment_A(const double & q2) const
                {
                    return 0.0;
                }
///////

                static NonlocalFormFactorPtr<nc::PToP> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nc::PToP>(new GvDV2020<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    const unsigned phiParam[4] = {3, 3, 2, 2}; //plus polarization

                    results.add({ real(1./this->phi(0.0, phiParam)), "Re{1/phi_+(q2 = 0.0)}" });
                    results.add({ imag(1./this->phi(0.0, phiParam)), "Im{1/phi_+(q2 = 0.0)}" });
                    results.add({ real(this->phi(16.0, phiParam)), "Re{phi_+(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phiParam)), "Im{phi_+(q2 = 16.0)}" });

                    const double s_0   = this->t_0();
                    const auto z1 = eos::nc_utils::z(1.0, 4.0 * pow(m_D0, 2), s_0);

                    results.add({ real(eos::nc_utils::PGvDV2020(z1, complex<double>(0.6,0.8), 2.0, 3.0, 4.0)),
                                "Re{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, 2.0, 3.0, 4.0)}" });
                    results.add({ imag(eos::nc_utils::PGvDV2020(z1, complex<double>(0.6,0.8), 2.0, 3.0, 4.0)),
                                "Im{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, 2.0, 3.0, 4.0)}" });

                    return results;
                }
        };



        /*
         * Parametrize the entire formfactor, i.e., both leading and all sub-leading powers as described in [GRvDV:2021].
         */
        template <typename Process_>
        class GRvDV2021 :
            public NonlocalFormFactor<nc::PToP>
        {
            public:
                std::shared_ptr<Model> model;

                //Polynomial expansion parameters
                UsedParameter re_alpha_0_plus;
                UsedParameter im_alpha_0_plus;
                UsedParameter re_alpha_1_plus;
                UsedParameter im_alpha_1_plus;
                UsedParameter re_alpha_2_plus;
                UsedParameter im_alpha_2_plus;

                //Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;
                UsedParameter m_Bsst;

                // final state meson parameters
                UsedParameter m_P;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;

                GRvDV2021(const Parameters & p, const Options & o) :
                    re_alpha_0_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_0^plus}@GRvDV2021"], *this),
                    im_alpha_0_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_0^plus}@GRvDV2021"], *this),
                    re_alpha_1_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_1^plus}@GRvDV2021"], *this),
                    im_alpha_1_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_1^plus}@GRvDV2021"], *this),
                    re_alpha_2_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_2^plus}@GRvDV2021"], *this),
                    im_alpha_2_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_2^plus}@GRvDV2021"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_d"], *this),
                    m_Bsst(p["mass::B_s^*"], *this),

                    m_P(p["mass::K_d"], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),
                    chiOPE(p["b->sccbar::chiOPE@GRvDV2021"], *this)
                {
                }

                ~GRvDV2021() = default;

                inline complex<double> phi(const double & q2, const unsigned phiParam[5]) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF                        a    b    c    d    e
                    // 0(P->P) aka plus          5    3    2    2    2
                    // perp(P->V) = par(P->V)    5    1    3    0    2
                    // 0(P->V) aka long          5    1    2    2    2

                    const double m_P2  = pow(m_P, 2);
                    const double m_Bsst2 = pow(m_Bsst, 2);
                    const double m_B2  = pow(m_B, 2),  m_B4 =  pow(m_B, 4);
                    const double m_D02 = pow(m_D0, 2), m_D04 = pow(m_D0, 4);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nc_utils::z(q2, 4.0 * pow(m_D0, 2), s_0);
                    const double Q2   = this->t_s();
                    const double chi = this->chiOPE();

                    const double a = phiParam[0], b = phiParam[1], c = phiParam[2], d = phiParam[3], e = phiParam[4];

                    const complex<double> Nlambda = 4. * M_PI * pow(m_B2, 0.5*(a-b+c+d-e) - 1.) * pow( 2.*(4.*m_D02-s_0)/3./chi, 0.5);
                    const complex<double> phi1 = -pow(2.*pow((4.*m_D02-Q2)*(4.*m_D02-s_0), 0.5) + 8.*m_D02 - Q2 - s_0, 0.5) /
                                                (2.*pow((4.*m_D02-Q2)*(4.*m_D02-s_0), 0.5) + 8.*m_D02 + Q2*(z-1.) - s_0*(z+1.));
                    const complex<double> phi2 = pow(m_B4*pow(z-1., 4.) - 2.*m_B2*pow(z-1., 2)*(-16*m_D02*z + m_P2*pow(z-1., 2) + s_0*pow(z+1., 2)) +
                                                pow(16*m_D02*z + m_P2*pow(z-1., 2) - s_0*pow(z+1., 2), 2), 0.5);
                    const complex<double> phi3 = pow(8*m_D02 + 4*pow(4*m_D04-s_0*m_D02, 0.5) - s_0, 0.5)/(-8*m_D02 - 4*pow(4*m_D04-s_0*m_D02, 0.5) + s_0*(z+1.));
                    const complex<double> phi4 = pow(s_0 * pow(z+1., 2.) - 16.*z*m_D02, -0.5);
                    const complex<double> phi5 = pow(s_0 * pow(z+1., 2.) - 16.*z*m_D02 - m_Bsst2*pow(-z+1., 2.), 0.5);

                    return Nlambda * pow(1.+z, 0.5) * pow(1.-z, a-b+c+d-e-1.5) * pow(phi1, a) * pow(phi2, 0.5*b) * pow(phi3, c) * pow(phi4, d) * pow(phi5, e);
                }

                inline complex<double> H_residue_jpsi(const unsigned phiParam[5], const complex<double> & alpha_0, const complex<double> & alpha_1,
                                                      const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z_Jpsi  = eos::nc_utils::z(m_Jpsi2,         s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(m_psi2S2,        s_p, s_0);

                    return eos::nc_utils::P(z_Jpsi, alpha_0, alpha_1, alpha_2) / phi(m_Jpsi2, phiParam) *
                            (1 - norm(z_Jpsi)) * (1. - z_Jpsi*std::conj(z_psi2S)) / (z_Jpsi - z_psi2S);
                };

                inline complex<double> H_residue_psi2s(const unsigned phiParam[5], const complex<double> & alpha_0, const complex<double> & alpha_1,
                                                       const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z_Jpsi  = eos::nc_utils::z(m_Jpsi2,         s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(m_psi2S2,        s_p, s_0);

                    return eos::nc_utils::P(z_psi2S, alpha_0, alpha_1, alpha_2) / phi(m_psi2S2, phiParam) *
                            (1 - norm(z_psi2S)) * (1. - z_psi2S*std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi);
                };


                virtual complex<double> H_plus(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,              s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(pow(m_Jpsi, 2),  s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(pow(m_psi2S, 2), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nc_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[5] = {5, 3, 2, 2, 2};

                    return eos::nc_utils::P(z, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                }


                virtual complex<double> Hhat_plus(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2, s_p, s_0);

                    return eos::nc_utils::P(z, alpha_0, alpha_1, alpha_2);
                }


                virtual complex<double> H_plus_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const unsigned phiParam[5] = {5, 3, 2, 2, 2};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_plus_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const unsigned phiParam[5] = {5, 3, 2, 2, 2};

                    return H_residue_psi2s(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> normalized_moment_A(const double & q2) const
                {
                    return 0.0;
                }

                static NonlocalFormFactorPtr<nc::PToP> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nc::PToP>(new GRvDV2021<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    const unsigned phiParam[5] = {5, 3, 2, 2, 2}; //plus polarization

                    results.add({ real(this->phi(16.0, phiParam)), "Re{phi_+(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phiParam)), "Im{phi_+(q2 = 16.0)}" });

                    const double s_0   = this->t_0();
                    const auto z1 = eos::nc_utils::z(1.0, 4.0 * pow(m_D0, 2), s_0);

                    results.add({ real(eos::nc_utils::P(z1, 2.0, 3.0, 4.0)), "Re{P(q2 = 1.0, 2.0, 3.0, 4.0)}" });
                    results.add({ imag(eos::nc_utils::P(z1, 2.0, 3.0, 4.0)), "Im{P(q2 = 1.0, 2.0, 3.0, 4.0)}" });
                    results.add({ real(eos::nc_utils::P(z1, complex<double>(2.0,5.0), complex<double>(3.0,6.0), complex<double>(4.0,7.0))),
                                "Re{P(q2 = 1.0, (2.0,5.0), (3.0,6.0), (4.0,7.0))}" });
                    results.add({ imag(eos::nc_utils::P(z1, complex<double>(2.0,5.0), complex<double>(3.0,6.0), complex<double>(4.0,7.0))),
                                "Im{P(q2 = 1.0, (2.0,5.0), (3.0,6.0), (4.0,7.0))}" });

                    return results;
                }
        };

    }

    NonlocalFormFactorPtr<nc::PToP>
    NonlocalFormFactor<nc::PToP>::make(const QualifiedName & name, const Parameters & p, const Options & o)
    {
        typedef QualifiedName KeyType;
        typedef std::function<NonlocalFormFactorPtr<nc::PToP> (const Parameters &, const Options &)> ValueType;
        std::map<KeyType, ValueType> entries
        {
            // trivial
            std::make_pair("B->K::naive",         &nc_p_to_p::Naive::make),
            // LCSR
            std::make_pair("B->K::LCSR",          &nc_p_to_p::LCSR::make),
            // parametrizations
            std::make_pair("B->K::GvDV2020",      &nc_p_to_p::GvDV2020<nc::BToK>::make),
            std::make_pair("B->K::GRvDV2021",     &nc_p_to_p::GRvDV2021<nc::BToK>::make),
        };

        auto i = entries.find(name);
        if (entries.end() == i)
        {
            return NonlocalFormFactorPtr<nc::PToP>(nullptr);
        }

        return i->second(p, o);
    }

    template <typename Process_>
    struct Implementation<NonlocalFormFactorObservable<Process_, nc::PToP>>
    {
        NameOption opt_formfactor;
        NonlocalFormFactorPtr<nc::PToP> nc;

        Implementation(const Parameters & p, const Options & o) :
            opt_formfactor(o, "formfactor", qnp::Name("GvDV2020")),
            nc(NonlocalFormFactor<nc::PToP>::make(QualifiedName(qnp::Prefix(Process_::label), opt_formfactor.value()), p, o))
        {
        }

        ~Implementation() = default;
    };

    template <typename Process_>
    NonlocalFormFactorObservable<Process_, nc::PToP>::NonlocalFormFactorObservable(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<NonlocalFormFactorObservable<Process_, nc::PToP>>(new Implementation<NonlocalFormFactorObservable<Process_, nc::PToP>>(p, o))
    {
    }

    template <typename Process_>
    NonlocalFormFactorObservable<Process_, nc::PToP>::~NonlocalFormFactorObservable() = default;

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToP>::re_H_plus(const double & q2) const
    {
        return real(this->_imp->nc->H_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToP>::im_H_plus(const double & q2) const
    {
        return imag(this->_imp->nc->H_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToP>::abs_H_plus(const double & q2) const
    {
        return abs(this->_imp->nc->H_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToP>::re_Hhat_plus(const double & q2) const
    {
        return real(this->_imp->nc->Hhat_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToP>::im_Hhat_plus(const double & q2) const
    {
        return imag(this->_imp->nc->Hhat_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToP>::abs_Hhat_plus(const double & q2) const
    {
        return abs(this->_imp->nc->Hhat_plus(q2));
    }


    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToP>::re_normalized_moment_A(const double & q2) const
    {
        return real(this->_imp->nc->normalized_moment_A(q2));
    }

    template class NonlocalFormFactorObservable<nc::BToK, nc::PToP>;
}
