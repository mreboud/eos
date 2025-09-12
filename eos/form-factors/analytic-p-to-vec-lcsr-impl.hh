/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_P_TO_VEC_LCSR_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_P_TO_VEC_LCSR_IMPL_HH 1

#include <eos/form-factors/analytic-p-to-vec-lcsr.hh>
#include <eos/form-factors/vec-lcdas.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/log.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/stringify.hh>

#include <functional>

namespace eos
{
    using namespace std::literals::string_literals;

    template <typename Process_>
    struct Implementation<AnalyticFormFactorPToVLCSR<Process_>>
    {
        std::shared_ptr<Model> model;

        // B-meson and final state meson parameters
        UsedParameter m_B;
        UsedParameter f_B;
        UsedParameter m_V;

        // sum rule parameters
        UsedParameter s0_0_A1;
        UsedParameter s0_1_A1;
        UsedParameter s0_0_A2;
        UsedParameter s0_1_A2;
        UsedParameter s0_0_A30;
        UsedParameter s0_1_A30;
        UsedParameter s0_0_V;
        UsedParameter s0_1_V;
        UsedParameter s0_0_T1;
        UsedParameter s0_1_T1;
        UsedParameter s0_0_T23A;
        UsedParameter s0_1_T23A;
        UsedParameter s0_0_T23B;
        UsedParameter s0_1_T23B;
        UsedParameter M2;

        // renormalization scale
        UsedParameter mu;

        std::shared_ptr<VectorLCDAs> lcdas;

        // // switches to enable/disable certain contributions
        TwistOption opt_2pt;
        std::array<double, 4> switch_2pt;
        // SwitchOption opt_3pt;
        // double switch_2pt_g;
        // double switch_3pt;

        // config for integration
        cubature::Config cub_conf;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make("SM", p, o)),
            m_B(p["mass::" + stringify(Process_::name_B)], u),
            f_B(p["decay-constant::" + stringify(Process_::name_B)], u),
            m_V(p["mass::" + stringify(Process_::name_V)], u),
            s0_0_A1(p[stringify(Process_::label) + "::s_0^A1,0@Kstar-LCSR"], u),
            s0_1_A1(p[stringify(Process_::label) + "::s_0^A1,1@Kstar-LCSR"], u),
            s0_0_A2(p[stringify(Process_::label) + "::s_0^A2,0@Kstar-LCSR"], u),
            s0_1_A2(p[stringify(Process_::label) + "::s_0^A2,1@Kstar-LCSR"], u),
            s0_0_A30(p[stringify(Process_::label) + "::s_0^A30,0@Kstar-LCSR"], u),
            s0_1_A30(p[stringify(Process_::label) + "::s_0^A30,1@Kstar-LCSR"], u),
            s0_0_V(p[stringify(Process_::label) + "::s_0^V,0@Kstar-LCSR"], u),
            s0_1_V(p[stringify(Process_::label) + "::s_0^V,1@Kstar-LCSR"], u),
            s0_0_T1(p[stringify(Process_::label) + "::s_0^T1,0@Kstar-LCSR"], u),
            s0_1_T1(p[stringify(Process_::label) + "::s_0^T1,1@Kstar-LCSR"], u),
            s0_0_T23A(p[stringify(Process_::label) + "::s_0^T23A,0@Kstar-LCSR"], u),
            s0_1_T23A(p[stringify(Process_::label) + "::s_0^T23A,1@Kstar-LCSR"], u),
            s0_0_T23B(p[stringify(Process_::label) + "::s_0^T23B,0@Kstar-LCSR"], u),
            s0_1_T23B(p[stringify(Process_::label) + "::s_0^T23B,1@Kstar-LCSR"], u),
            M2(p[stringify(Process_::label) + "::M^2@Kstar-LCSR"], u),
            mu(p[stringify(Process_::label) + "::mu@Kstar-LCSR"], u),
            lcdas(VectorLCDAs::make(Process_::lcsr_V, p, o)),
            opt_2pt(o, options, "2pt-twist"_ok),
            cub_conf(cubature::Config().epsrel(1e-3))
        {
            u.uses(*lcdas);

            switch_2pt[0] = (opt_2pt.value() && Twist::two);
            switch_2pt[1] = (opt_2pt.value() && Twist::three);
            switch_2pt[2] = (opt_2pt.value() && Twist::four);
            switch_2pt[3] = (opt_2pt.value() && Twist::five);
        }

        ~Implementation() = default;

        /* quark masses for the propagating quark */
        double m_b() const
        {
            return model->m_b_msbar(mu());
        }

        /* forwarding the LCDAs */
        // {{{
        inline
        double f_V_perp() const
        {
            return lcdas->fperp(mu());
        }
        inline
        double f_V_para() const
        {
            return lcdas->fpara();
        }

        inline
        double phi2perp(const double & u) const
        {
            return switch_2pt[0] * lcdas->phi2perp(u, mu());
        }
        inline
        double Barphi2perp(const double & u) const
        {
            return switch_2pt[0] * lcdas->Barphi2perp(u, mu());
        }
        inline
        double BarBarphi2perp(const double & u) const
        {
            return switch_2pt[0] * lcdas->BarBarphi2perp(u, mu());
        }
        inline
        double phi2para(const double & u) const
        {
            return switch_2pt[0] * lcdas->phi2para(u, mu());
        }
        inline
        double Barphi2para(const double & u) const
        {
            return switch_2pt[0] * lcdas->Barphi2para(u, mu());
        }
        inline
        double BarBarphi2para(const double & u) const
        {
            return switch_2pt[0] * lcdas->BarBarphi2para(u, mu());
        }
        inline
        double phi3perp(const double & u) const
        {
            return switch_2pt[1] * lcdas->phi3perp(u, mu());
        }
        inline
        double Barphi3perp(const double & u) const
        {
            return switch_2pt[1] * lcdas->Barphi3perp(u, mu());
        }
        inline
        double BarBarphi3perp(const double & u) const
        {
            return switch_2pt[1] * lcdas->BarBarphi3perp(u, mu());
        }
        inline
        double phi3para(const double & u) const
        {
            return switch_2pt[1] * lcdas->phi3para(u, mu());
        }
        inline
        double Barphi3para(const double & u) const
        {
            return switch_2pt[1] * lcdas->Barphi3para(u, mu());
        }
        inline
        double BarBarphi3para(const double & u) const
        {
            return switch_2pt[1] * lcdas->BarBarphi3para(u, mu());
        }
        inline
        double psi3perp(const double & u) const
        {
            return switch_2pt[1] * lcdas->psi3perp(u, mu());
        }
        inline
        double psi3para(const double & u) const
        {
            return switch_2pt[1] * lcdas->psi3para(u, mu());
        }
        inline
        double phi4perp(const double & u) const
        {
            return switch_2pt[2] * lcdas->phi4perp(u, mu());
        }
        inline
        double phi4perpprime(const double & u) const
        {
            return switch_2pt[2] * lcdas->phi4perpprime(u, mu());;
        }
        inline
        double phi4para(const double & u) const
        {
            return switch_2pt[2] * lcdas->phi4para(u, mu());
        }
        inline
        double phi4paraprime(const double & u) const
        {
            return switch_2pt[2] * lcdas->phi4paraprime(u, mu());;
        }
        inline
        double Barphi4para(const double & u) const
        {
            std::function<double (const double &)> integrand = [this](const double & v) -> double
            {
                return switch_2pt[2] * phi4para(v);
            };
            return integrate(integrand, 0.0, u, cub_conf);
        }
        inline
        double psi4perp(const double & u) const
        {
            return switch_2pt[2] * lcdas->psi4perp(u, mu());
        }
        inline
        double Barpsi4perp(const double & u) const
        {
            std::function<double (const double &)> integrand = [this](const double & v) -> double
            {
                return switch_2pt[2] * psi4perp(v);
            };
            return integrate(integrand, 0.0, u, cub_conf);
        }
        inline
        double BarBarpsi4perp(const double & u) const
        {
            // int_[0,u] int_[0,v] f(w) dw dv = u * int_[0,u] (u - v) * f(v) dv
            std::function<double (const double &)> integrand = [this, u](const double & v) -> double
            {
                return switch_2pt[2] * (u - v) * psi4perp(v);
            };
            return integrate(integrand, 0.0, u, cub_conf);
        }
        inline
        double psi4para(const double & u) const
        {
            return switch_2pt[2] * lcdas->psi4para(u, mu());
        }
        inline
        double Barpsi4para(const double & u) const
        {
            std::function<double (const double &)> integrand = [this](const double & v) -> double
            {
                return switch_2pt[2] * psi4para(v);
            };
            return integrate(integrand, 0.0, u, cub_conf);
        }
        inline
        double BarBarpsi4para(const double & u) const
        {
            std::function<double (const double &)> integrand = [this, u](const double & v) -> double
            {
                return switch_2pt[2] * (u - v) * psi4para(v);
            };
            return integrate(integrand, 0.0, u, cub_conf);
        }
        inline
        double phi5perp(const double & u) const
        {
            return switch_2pt[3] * lcdas->phi5perp(u, mu());
        }
        inline
        double phi5perpprime(const double &) const
        {
            return switch_2pt[3] * 0.0;
        }
        inline
        double Barphi5perp(const double & u) const
        {
            std::function<double (const double &)> integrand = [this](const double & v) -> double
            {
                return switch_2pt[3] * phi5perp(v);
            };
            return integrate(integrand, 0.0, u, cub_conf);
        }
        inline
        double psi5perp(const double & u) const
        {
            return switch_2pt[3] * lcdas->psi5perp(u, mu());
        }
        inline
        double psi5perpprime(const double &) const
        {
            return switch_2pt[3] * 0.0;
        }
        inline
        double psi5perpsecond(const double &) const
        {
            return switch_2pt[3] * 0.0;
        }
        // }}}


        /* auxilliary functions */
        // {{{
        double s(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b;

            return (m_b2 + u * (1.0 - u) * m_V2 - (1.0 - u) * q2) / u;
        }

        double u0(const double & s0, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b;
            const double x = m_V2 + q2 - s0;

            const double result = (x + sqrt(4.0 * m_V2 * (m_b2 - q2) + x * x)) / (2.0 * m_V2);

            if ((result < 0) || (result > 1))
            {
                throw InternalError("The threshold value u0 = " + stringify(result) + " lies outside the [0, 1] interval.");
            }

            return result;
        }

        // eta = du / ds
        double eta(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;

            return - u2 / (m_b2 - q2 + u2 * m_V2);
        }

        double etaprime(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;

            return 2.0 * u * (q2 - m_b2) / (m_b2 - q2 + u2 * m_V2) / (m_b2 - q2 + u2 * m_V2);
        }

        double etasecond(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;

            return 2.0 * (m_b2 - q2 - 3.0 * u2 * m_V2) * (q2 - m_b2) / power_of<3>(m_b2 - q2 + u2 * m_V2);
        }

        const std::array<double, 4> I_perp(const double & u, const double &) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            std::array<double, 4> result{{0, 0, 0, 0}};

            const double phi4perp_u = phi4perp(u);

            // I_perp_0
            result[0] = phi2perp(u) / u;
            // I_perp_1
            result[1] = -0.25 * (m_V2 * phi4perp_u) / u2 + (f_V_para * m_b * m_V * psi3perp(u)) / (2.0 * f_V_perp * u2);
            // I_perp_2
            result[2] = -0.25 * (m_b2 * m_V2 * phi4perp_u) / power_of<3>(u);
            // I_perp_3
            result[3] = -0.125 * (f_V_para * power_of<3>(m_b * m_V) * psi5perp(u)) / (f_V_perp * power_of<4>(u));

            return result;
        }

        const std::array<std::array<double, 4>, 4> D_I_perp(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            const std::array<double, 4> I_perp = this->I_perp(u, q2);
            // result[n][m] = D_eta^m[I_n]
            std::array<std::array<double, 4>, 4> result{{{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}};

            const double psi5perp_u = psi5perp(u), psi5perpprime_u = psi5perpprime(u);

            // I_perp_0
            result[0][0] = I_perp[0];
            // I_perp_1
            result[1][0] = I_perp[1];
            // I_perp_2
            result[2][0] = I_perp[2];
            result[2][1] = I_perp[2] * etaprime(u, q2) + (3.0 * m_b2 * m_V2 * eta(u, q2) * phi4perp(u)) / (4.0 * power_of<4>(u)) - (m_b2 * m_V2 * eta(u, q2) * phi4perpprime(u)) / (4.0 * power_of<3>(u));
            // I_perp_3
            result[3][0] = I_perp[3];
            result[3][1] = I_perp[3] * etaprime(u, q2) + (f_V_para * power_of<3>(m_b * m_V) * eta(u, q2) * psi5perp_u) / (2.0 * f_V_perp * power_of<5>(u)) - (f_V_para * power_of<3>(m_b * m_V) * eta(u, q2) * psi5perpprime_u) / (8.0 * f_V_perp * power_of<4>(u));
            result[3][2] = result[3][1] * etaprime(u, q2) + eta(u, q2) * (I_perp[3] * etasecond(u, q2) + (-5 * f_V_para * power_of<3>(m_b * m_V) * eta(u, q2) * psi5perp_u) / (2.0 * f_V_perp * power_of<6>(u)) + (f_V_para * power_of<3>(m_b * m_V) * psi5perp_u * etaprime(u, q2)) / (f_V_perp * power_of<5>(u)) + (f_V_para * power_of<3>(m_b * m_V) * eta(u, q2) * psi5perpprime_u) / (f_V_perp * power_of<5>(u)) - (f_V_para * power_of<3>(m_b * m_V) * etaprime(u, q2) * psi5perpprime_u) / (4.0 * f_V_perp * power_of<4>(u)) - (f_V_para * power_of<3>(m_b * m_V) * eta(u, q2) * psi5perpsecond(u)) / (8.0 * f_V_perp * power_of<4>(u)));

            return result;
        }


        const std::array<double, 4> I_eta(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            std::array<double, 4> result{{0, 0, 0, 0}};

            const double phi2perp_u = phi2perp(u), phi4perp_u = phi4perp(u);
            const double Barphi2perp_u = Barphi2perp(u), Barpsi4perp_u = Barpsi4perp(u);
            const double BarBarphi2perp_u = BarBarphi2perp(u), BarBarphi3para_u = BarBarphi3para(u), BarBarpsi4perp_u = BarBarpsi4perp(u);

            // I_eta_0
            result[0] = (-4.0 * f_V_perp * m_V2 * BarBarphi2perp_u + 8.0 * f_V_perp * m_V2 * BarBarphi3para_u - 4.0 * f_V_perp * m_V2 * BarBarpsi4perp_u + 4.0 * f_V_perp * m_V2 * u * Barphi2perp_u - 4.0 * f_V_perp * m_V2 * u * Barpsi4perp_u + 4.0 * f_V_perp * m_b2 * phi2perp_u - 4.0 * f_V_perp * q2 * phi2perp_u + 4.0 * f_V_perp * m_V2 * u2 * phi2perp_u + 8.0 * f_V_para * m_b * m_V * u * phi3perp(u) + f_V_perp * m_V2 * phi4perp_u + 4.0 * f_V_perp * m_V2 * u * psi3para(u)) / (8.0 * f_V_perp * u2);
            // I_eta_1
            result[1] = -((f_V_para * m_b * power_of<3>(m_V) * BarBarphi2para(u)) / (f_V_perp * u2)) + (m_V2 * (m_b2 - q2 + m_V2 * u2) * BarBarphi2perp_u) / (2.0 * power_of<3>(u)) - (m_V2 * (m_b2 - q2 + m_V2 * u2) * BarBarphi3para_u) / power_of<3>(u) + (2.0 * f_V_para * m_b * power_of<3>(m_V) * BarBarphi3perp(u)) / (f_V_perp * u2) - (f_V_para * m_b * power_of<3>(m_V) * BarBarpsi4para(u)) / (f_V_perp * u2) + (m_V2 * (m_b2 - q2 + m_V2 * u2) * BarBarpsi4perp_u) / (2.0 * power_of<3>(u)) + (m_b2 * m_V2 * Barphi2perp_u) / u2 - (m_b2 * m_V2 * Barpsi4perp_u) / u2 + (m_V2 * (m_b2 + q2 - m_V2 * u2) * phi4perp_u) / (8.0 * power_of<3>(u));
            // I_eta_2
            result[2] = -0.125 * (m_b2 * m_V2 * (m_b2 - q2 + m_V2 * u2) * phi4perp_u) / power_of<4>(u) - (f_V_para * power_of<3>(m_b * m_V) * phi5perp(u)) / (4.0 * f_V_perp * power_of<3>(u));

            return result;
        }

        const std::array<std::array<double, 4>, 4> D_I_eta(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            const std::array<double, 4> I_eta = this->I_eta(u, q2);
            // result[n][m] = D_eta^m[I_n]
            std::array<std::array<double, 4>, 4> result{{{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}};

            // I_eta_0
            result[0][0] = I_eta[0];
            // I_eta_1
            result[1][0] = I_eta[1];
            // I_eta_2
            result[2][0] = I_eta[2];
            result[2][1] = I_eta[2] * etaprime(u, q2) + eta(u, q2) * (-0.25 * (m_b2 * power_of<4>(m_V) * phi4perp(u)) / power_of<3>(u) + (m_b2 * m_V2 * (m_b2 - q2 + m_V2 * u2) * phi4perp(u)) / (2.0 * power_of<5>(u)) + (3.0 * f_V_para * power_of<3>(m_b * m_V) * phi5perp(u)) / (4.0 * f_V_perp * power_of<4>(u)) - (m_b2 * m_V2 * (m_b2 - q2 + m_V2 * u2) * phi4perpprime(u)) / (8.0 * power_of<4>(u)) - (f_V_para * power_of<3>(m_b * m_V) * phi5perpprime(u)) / (4.0 * f_V_perp * power_of<3>(u)));

            return result;
        }


        const std::array<double, 4> I_q(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            std::array<double, 4> result{{0, 0, 0, 0}};

            const double phi4perp_u = phi4perp(u);
            const double BarBarphi2perp_u = BarBarphi2perp(u), BarBarphi3para_u = BarBarphi3para(u), BarBarpsi4perp_u = BarBarpsi4perp(u);

            // I_q_0
            result[0] = phi2perp(u) / u;
            // I_q_1
            result[1] = (m_V2 * (-2.0 + 3.0 * u) * BarBarphi2perp_u) / power_of<3>(u) - (2.0 * m_V2 * (-2.0 + 3.0 * u) * BarBarphi3para_u) / power_of<3>(u) + (m_V2 * (-2.0 + 3.0 * u) * BarBarpsi4perp_u) / power_of<3>(u) + (m_V2 * (-1.0 + u) * Barphi2perp(u)) / u2 - (2.0 * f_V_para * m_b * m_V * Barphi3perp(u)) / (f_V_perp * u2) - (m_V2 * (-1.0 + u) * Barpsi4perp(u)) / u2 + (2.0 * f_V_para * m_b * m_V * phi2para(u)) / (f_V_perp * u2) - (m_V2 * phi4perp_u) / (4.0 * u2) - (m_V2 * (-1.0 + u) * psi3para(u)) / u2;
            // I_q_2
            result[2] = (2.0 * f_V_para * m_b * power_of<3>(m_V) * (-1.0 + u) * BarBarphi2para(u)) / (f_V_perp * power_of<3>(u)) + (m_V2 * (q2 * (-1.0 + u) - m_V2 * (-1.0 + u) * u2 + m_b2 * (1 + u)) * BarBarphi2perp_u) / power_of<4>(u) - (2.0 * m_V2 * (q2 * (-1.0 + u) - m_V2 * (-1.0 + u) * u2 + m_b2 * (1 + u)) * BarBarphi3para_u) / power_of<4>(u) - (4.0 * f_V_para * m_b * power_of<3>(m_V) * (-1.0 + u) * BarBarphi3perp(u)) / (f_V_perp * power_of<3>(u)) + (2.0 * f_V_para * m_b * power_of<3>(m_V) * (-1.0 + u) * BarBarpsi4para(u)) / (f_V_perp * power_of<3>(u)) + (m_V2 * (q2 * (-1.0 + u) - m_V2 * (-1.0 + u) * u2 + m_b2 * (1 + u)) * BarBarpsi4perp_u) / power_of<4>(u) - (m_b2 * m_V2 * phi4perp_u) / (4.0 * power_of<3>(u));
            // I_q_3
            result[3] = -0.5 * (f_V_para * power_of<3>(m_b * m_V) * Barphi4para(u)) / (f_V_perp * power_of<4>(u)) + (f_V_para * power_of<3>(m_b * m_V) * Barphi5perp(u)) / (2.0 * f_V_perp * power_of<4>(u));

            return result;
        }

        const std::array<std::array<double, 4>, 4> D_I_q(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            const std::array<double, 4> I_q = this->I_q(u, q2);
            // result[n][m] = D_eta^m[I_n]
            std::array<std::array<double, 4>, 4> result{{{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}};

            // I_q_0
            result[0][0] = I_q[0];
            // I_q_1
            result[1][0] = I_q[1];
            // I_q_2
            result[2][0] = I_q[2];
            result[2][1] = I_q[2] * etaprime(u, q2) + eta(u, q2) * ((-6.0 * f_V_para * m_b * power_of<3>(m_V) * (-1.0 + u) * BarBarphi2para(u)) / (f_V_perp * power_of<4>(u)) + (2.0 * f_V_para * m_b * power_of<3>(m_V) * BarBarphi2para(u)) / (f_V_perp * power_of<3>(u)) + (m_V2 * (m_b2 + q2 - 2.0 * m_V2 * (-1.0 + u) * u - m_V2 * u2) * BarBarphi2perp(u)) / power_of<4>(u) - (4.0 * m_V2 * (q2 * (-1.0 + u) - m_V2 * (-1.0 + u) * u2 + m_b2 * (1 + u)) * BarBarphi2perp(u)) / power_of<5>(u) - (2.0 * m_V2 * (m_b2 + q2 - 2.0 * m_V2 * (-1.0 + u) * u - m_V2 * u2) * BarBarphi3para(u)) / power_of<4>(u) + (8.0 * m_V2 * (q2 * (-1.0 + u) - m_V2 * (-1.0 + u) * u2 + m_b2 * (1 + u)) * BarBarphi3para(u)) / power_of<5>(u) + (12 * f_V_para * m_b * power_of<3>(m_V) * (-1.0 + u) * BarBarphi3perp(u)) / (f_V_perp * power_of<4>(u)) - (4.0 * f_V_para * m_b * power_of<3>(m_V) * BarBarphi3perp(u)) / (f_V_perp * power_of<3>(u)) - (6 * f_V_para * m_b * power_of<3>(m_V) * (-1.0 + u) * BarBarpsi4para(u)) / (f_V_perp * power_of<4>(u)) + (2.0 * f_V_para * m_b * power_of<3>(m_V) * BarBarpsi4para(u)) / (f_V_perp * power_of<3>(u)) + (m_V2 * (m_b2 + q2 - 2.0 * m_V2 * (-1.0 + u) * u - m_V2 * u2) * BarBarpsi4perp(u)) / power_of<4>(u) - (4.0 * m_V2 * (q2 * (-1.0 + u) - m_V2 * (-1.0 + u) * u2 + m_b2 * (1 + u)) * BarBarpsi4perp(u)) / power_of<5>(u) + (3.0 * m_b2 * m_V2 * phi4perp(u)) / (4.0 * power_of<4>(u)) + (2.0 * f_V_para * m_b * power_of<3>(m_V) * (-1.0 + u) * Barphi2para(u)) / (f_V_perp * power_of<3>(u)) + (m_V2 * (q2 * (-1.0 + u) - m_V2 * (-1.0 + u) * u2 + m_b2 * (1 + u)) * Barphi2perp(u)) / power_of<4>(u) - (2.0 * m_V2 * (q2 * (-1.0 + u) - m_V2 * (-1.0 + u) * u2 + m_b2 * (1 + u)) * Barphi3para(u)) / power_of<4>(u) - (4.0 * f_V_para * m_b * power_of<3>(m_V) * (-1.0 + u) * Barphi3perp(u)) / (f_V_perp * power_of<3>(u)) + (2.0 * f_V_para * m_b * power_of<3>(m_V) * (-1.0 + u) * Barpsi4para(u)) / (f_V_perp * power_of<3>(u)) + (m_V2 * (q2 * (-1.0 + u) - m_V2 * (-1.0 + u) * u2 + m_b2 * (1 + u)) * Barpsi4perp(u)) / power_of<4>(u) - (m_b2 * m_V2 * phi4perpprime(u)) / (4.0 * power_of<3>(u)));
            // I_q_3
            result[3][0] = I_q[3];
            result[3][1] = I_q[3] * etaprime(u, q2) + eta(u, q2) * ((2.0 * f_V_para * power_of<3>(m_b * m_V) * Barphi4para(u)) / (f_V_perp * power_of<5>(u)) - (2.0 * f_V_para * power_of<3>(m_b * m_V) * Barphi5perp(u)) / (f_V_perp * power_of<5>(u)) - (f_V_para * power_of<3>(m_b * m_V) * phi4para(u)) / (2.0 * f_V_perp * power_of<4>(u)) + (f_V_para * power_of<3>(m_b * m_V) * phi5perp(u)) / (2.0 * f_V_perp * power_of<4>(u)));
            result[3][2] = result[3][1] * etaprime(u, q2) + eta(u, q2) * (I_q[3] * etasecond(u, q2) + 2.0 * ((2.0 * f_V_para * power_of<3>(m_b * m_V) * Barphi4para(u)) / (f_V_perp * power_of<5>(u)) - (2.0 * f_V_para * power_of<3>(m_b * m_V) * Barphi5perp(u)) / (f_V_perp * power_of<5>(u)) - (f_V_para * power_of<3>(m_b * m_V) * phi4para(u)) / (2.0 * f_V_perp * power_of<4>(u)) + (f_V_para * power_of<3>(m_b * m_V) * phi5perp(u)) / (2.0 * f_V_perp * power_of<4>(u))) * etaprime(u, q2) + eta(u, q2) * ((-10 * f_V_para * power_of<3>(m_b * m_V) * Barphi4para(u)) / (f_V_perp * power_of<6>(u)) + (10 * f_V_para * power_of<3>(m_b * m_V) * Barphi5perp(u)) / (f_V_perp * power_of<6>(u)) + (4.0 * f_V_para * power_of<3>(m_b * m_V) * phi4para(u)) / (f_V_perp * power_of<5>(u)) - (4.0 * f_V_para * power_of<3>(m_b * m_V) * phi5perp(u)) / (f_V_perp * power_of<5>(u)) - (f_V_para * power_of<3>(m_b * m_V) * phi4paraprime(u)) / (2.0 * f_V_perp * power_of<4>(u)) + (f_V_para * power_of<3>(m_b * m_V) * phi5perpprime(u)) / (2.0 * f_V_perp * power_of<4>(u))));

            return result;
        }

        const std::array<double, 4> I_p(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            std::array<double, 4> result{{0, 0, 0, 0}};

            const double phi4perp_u = phi4perp(u);
            const double BarBarphi2perp_u = BarBarphi2perp(u), BarBarphi3para_u = BarBarphi3para(u), BarBarpsi4perp_u = BarBarpsi4perp(u);

            // I_p_0
            result[0] = -(phi2perp(u) / u);
            // I_p_1
            result[1] = (-3.0 * m_V2 * BarBarphi2perp_u) / u2 + (6 * m_V2 * BarBarphi3para_u) / u2 - (3.0 * m_V2 * BarBarpsi4perp_u) / u2 - (m_V2 * Barphi2perp(u)) / u + (2.0 * f_V_para * m_b * m_V * Barphi3perp(u)) / (f_V_perp * u2) + (m_V2 * Barpsi4perp(u)) / u - (2.0 * f_V_para * m_b * m_V * phi2para(u)) / (f_V_perp * u2) + (m_V2 * phi4perp_u) / (4.0 * u2) + (m_V2 * psi3para(u)) / u;
            // I_p_2
            result[2] = (-2.0 * f_V_para * m_b * power_of<3>(m_V) * BarBarphi2para(u)) / (f_V_perp * u2) + (m_V2 * (-m_b2 - q2 + m_V2 * u2) * BarBarphi2perp_u) / power_of<3>(u) + (2.0 * m_V2 * (m_b2 + q2 - m_V2 * u2) * BarBarphi3para_u) / power_of<3>(u) + (4.0 * f_V_para * m_b * power_of<3>(m_V) * BarBarphi3perp(u)) / (f_V_perp * u2) - (2.0 * f_V_para * m_b * power_of<3>(m_V) * BarBarpsi4para(u)) / (f_V_perp * u2) + (m_V2 * (-m_b2 - q2 + m_V2 * u2) * BarBarpsi4perp_u) / power_of<3>(u) + (m_b2 * m_V2 * phi4perp_u) / (4.0 * power_of<3>(u));
            // I_p_3
            result[3] = (f_V_para * power_of<3>(m_b * m_V) * Barphi4para(u)) / (2.0 * f_V_perp * power_of<4>(u)) - (f_V_para * power_of<3>(m_b * m_V) * Barphi5perp(u)) / (2.0 * f_V_perp * power_of<4>(u));

            return result;
        }

        const std::array<std::array<double, 4>, 4> D_I_p(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            const std::array<double, 4> I_p = this->I_p(u, q2);
            // result[n][m] = D_eta^m[I_n]
            std::array<std::array<double, 4>, 4> result{{{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}};

            // I_p_0
            result[0][0] = I_p[0];
            // I_p_1
            result[1][0] = I_p[1];
            // I_p_2
            result[2][0] = I_p[2];
            result[2][1] = I_p[2] * etaprime(u, q2) + eta(u, q2) * ((4.0 * f_V_para * m_b * power_of<3>(m_V) * BarBarphi2para(u)) / (f_V_perp * power_of<3>(u)) + (2.0 * power_of<4>(m_V) * BarBarphi2perp(u)) / u2 - (3.0 * m_V2 * (-m_b2 - q2 + m_V2 * u2) * BarBarphi2perp(u)) / power_of<4>(u) - (4.0 * power_of<4>(m_V) * BarBarphi3para(u)) / u2 - (6 * m_V2 * (m_b2 + q2 - m_V2 * u2) * BarBarphi3para(u)) / power_of<4>(u) - (8.0 * f_V_para * m_b * power_of<3>(m_V) * BarBarphi3perp(u)) / (f_V_perp * power_of<3>(u)) + (4.0 * f_V_para * m_b * power_of<3>(m_V) * BarBarpsi4para(u)) / (f_V_perp * power_of<3>(u)) + (2.0 * power_of<4>(m_V) * BarBarpsi4perp(u)) / u2 - (3.0 * m_V2 * (-m_b2 - q2 + m_V2 * u2) * BarBarpsi4perp(u)) / power_of<4>(u) - (3.0 * m_b2 * m_V2 * phi4perp(u)) / (4.0 * power_of<4>(u)) - (2.0 * f_V_para * m_b * power_of<3>(m_V) * Barphi2para(u)) / (f_V_perp * u2) + (m_V2 * (-m_b2 - q2 + m_V2 * u2) * Barphi2perp(u)) / power_of<3>(u) + (2.0 * m_V2 * (m_b2 + q2 - m_V2 * u2) * Barphi3para(u)) / power_of<3>(u) + (4.0 * f_V_para * m_b * power_of<3>(m_V) * Barphi3perp(u)) / (f_V_perp * u2) - (2.0 * f_V_para * m_b * power_of<3>(m_V) * Barpsi4para(u)) / (f_V_perp * u2) + (m_V2 * (-m_b2 - q2 + m_V2 * u2) * Barpsi4perp(u)) / power_of<3>(u) + (m_b2 * m_V2 * phi4perpprime(u)) / (4.0 * power_of<3>(u)));
            // I_p_3
            result[3][0] = I_p[3];
            result[3][1] = I_p[3] * etaprime(u, q2) + eta(u, q2) * ((-2.0 * f_V_para * power_of<3>(m_b * m_V) * Barphi4para(u)) / (f_V_perp * power_of<5>(u)) + (2.0 * f_V_para * power_of<3>(m_b * m_V) * Barphi5perp(u)) / (f_V_perp * power_of<5>(u)) + (f_V_para * power_of<3>(m_b * m_V) * phi4para(u)) / (2.0 * f_V_perp * power_of<4>(u)) - (f_V_para * power_of<3>(m_b * m_V) * phi5perp(u)) / (2.0 * f_V_perp * power_of<4>(u)));
            result[3][2] = result[3][1] * etaprime(u, q2) + eta(u, q2) * (I_p[3] * etasecond(u, q2) + 2.0 * ((-2.0 * f_V_para * power_of<3>(m_b * m_V) * Barphi4para(u)) / (f_V_perp * power_of<5>(u)) + (2.0 * f_V_para * power_of<3>(m_b * m_V) * Barphi5perp(u)) / (f_V_perp * power_of<5>(u)) + (f_V_para * power_of<3>(m_b * m_V) * phi4para(u)) / (2.0 * f_V_perp * power_of<4>(u)) - (f_V_para * power_of<3>(m_b * m_V) * phi5perp(u)) / (2.0 * f_V_perp * power_of<4>(u))) * etaprime(u, q2) + eta(u, q2) * ((10 * f_V_para * power_of<3>(m_b * m_V) * Barphi4para(u)) / (f_V_perp * power_of<6>(u)) - (10 * f_V_para * power_of<3>(m_b * m_V) * Barphi5perp(u)) / (f_V_perp * power_of<6>(u)) - (4.0 * f_V_para * power_of<3>(m_b * m_V) * phi4para(u)) / (f_V_perp * power_of<5>(u)) + (4.0 * f_V_para * power_of<3>(m_b * m_V) * phi5perp(u)) / (f_V_perp * power_of<5>(u)) + (f_V_para * power_of<3>(m_b * m_V) * phi4paraprime(u)) / (2.0 * f_V_perp * power_of<4>(u)) - (f_V_para * power_of<3>(m_b * m_V) * phi5perpprime(u)) / (2.0 * f_V_perp * power_of<4>(u))));

            return result;
        }

        const std::array<double, 4> I_T_perp(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            std::array<double, 4> result{{0, 0, 0, 0}};

            const double psi3perp_u = psi3perp(u), phi5perp_u = phi5perp(u), psi5perp_u = psi5perp(u);
            const double Barphi4para_u = Barphi4para(u), Barphi5perp_u = Barphi5perp(u);

            // I_T_perp_0
            result[0] = (-4.0 * f_V_para * m_V * Barphi3perp(u) + 4.0 * f_V_para * m_V * phi2para(u) + 4.0 * f_V_perp * m_b * phi2perp(u) + 4.0 * f_V_para * m_V * u * phi3perp(u) + f_V_para * m_V * psi3perp_u) / (4.0 * f_V_perp * u);
            // I_T_perp_1
            result[1] = (m_b * m_V2 * BarBarphi2perp(u)) / u2 - (2.0 * m_b * m_V2 * BarBarphi3para(u)) / u2 + (m_b * m_V2 * BarBarpsi4perp(u)) / u2 + (m_b * m_V2 * Barphi2perp(u)) / u - (f_V_para * power_of<3>(m_V) * Barphi4para_u) / (4.0 * f_V_perp * u2) + (f_V_para * power_of<3>(m_V) * Barphi5perp_u) / (4.0 * f_V_perp * u2) - (m_b * m_V2 * Barpsi4perp(u)) / u - (f_V_para * power_of<3>(m_V) * phi5perp_u) / (4.0 * f_V_perp * u) + (f_V_para * m_V * (m_b2 + q2 - m_V2 * u2) * psi3perp_u) / (4.0 * f_V_perp * u2);
            // I_T_perp_2
            result[2] = -0.25 * (f_V_para * m_b2 * power_of<3>(m_V) * Barphi4para_u) / (f_V_perp * power_of<3>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * Barphi5perp_u) / (4.0 * f_V_perp * power_of<3>(u)) - (power_of<3>(m_b) * m_V2 * phi4perp(u)) / (4.0 * power_of<3>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * phi5perp_u) / (4.0 * f_V_perp * u2) + (f_V_para * (-(power_of<3>(m_V) * q2) + power_of<5>(m_V) * u2) * psi5perp_u) / (16.0 * f_V_perp * power_of<3>(u));
            // I_T_perp_3
            result[3] = -1.0 / 16.0 * (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 + q2 - m_V2 * u2) * psi5perp_u) / (f_V_perp * power_of<4>(u));

            return result;
        }

        const std::array<std::array<double, 4>, 4> D_I_T_perp(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            const std::array<double, 4> I_T_perp = this->I_T_perp(u, q2);
            // result[n][m] = D_eta^m[I_n]
            std::array<std::array<double, 4>, 4> result{{{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}};

            // I_T_perp_0
            result[0][0] = I_T_perp[0];
            // I_T_perp_1
            result[1][0] = I_T_perp[1];
            // I_T_perp_2
            result[2][0] = I_T_perp[2];
            result[2][1] = I_T_perp[2] * etaprime(u, q2) + eta(u, q2) * ((3.0 * f_V_para * m_b2 * power_of<3>(m_V) * Barphi4para(u)) / (4.0 * f_V_perp * power_of<4>(u)) - (3.0 * f_V_para * m_b2 * power_of<3>(m_V) * Barphi5perp(u)) / (4.0 * f_V_perp * power_of<4>(u)) + (3.0 * power_of<3>(m_b) * m_V2 * phi4perp(u)) / (4.0 * power_of<4>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * phi5perp(u)) / (2.0 * f_V_perp * power_of<3>(u)) + (f_V_para * power_of<5>(m_V) * psi5perp(u)) / (8.0 * f_V_perp * u2) - (3.0 * f_V_para * (-(power_of<3>(m_V) * q2) + power_of<5>(m_V) * u2) * psi5perp(u)) / (16.0 * f_V_perp * power_of<4>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * phi4para(u)) / (4.0 * f_V_perp * power_of<3>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * phi5perp(u)) / (4.0 * f_V_perp * power_of<3>(u)) - (power_of<3>(m_b) * m_V2 * phi4perpprime(u)) / (4.0 * power_of<3>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * phi5perpprime(u)) / (4.0 * f_V_perp * u2) + (f_V_para * (-(power_of<3>(m_V) * q2) + power_of<5>(m_V) * u2) * psi5perpprime(u)) / (16.0 * f_V_perp * power_of<3>(u)));
            // I_T_perp_3
            result[3][0] = I_T_perp[3];
            result[3][1] = I_T_perp[3] * etaprime(u, q2) + (f_V_para * m_b2 * power_of<5>(m_V) * eta(u, q2) * psi5perp(u)) / (8.0 * f_V_perp * power_of<3>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 + q2 - m_V2 * u2) * eta(u, q2) * psi5perp(u)) / (4.0 * f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 + q2 - m_V2 * u2) * eta(u, q2) * psi5perpprime(u)) / (16.0 * f_V_perp * power_of<4>(u));
            result[3][2] = result[3][1] * etaprime(u, q2) + eta(u, q2) * (I_T_perp[3] * etasecond(u, q2) + (-7.0 * f_V_para * m_b2 * power_of<5>(m_V) * eta(u, q2) * psi5perp(u)) / (8.0 * f_V_perp * power_of<4>(u)) - (5 * f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 + q2 - m_V2 * u2) * eta(u, q2) * psi5perp(u)) / (4.0 * f_V_perp * power_of<6>(u)) + (f_V_para * m_b2 * power_of<5>(m_V) * psi5perp(u) * etaprime(u, q2)) / (4.0 * f_V_perp * power_of<3>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 + q2 - m_V2 * u2) * psi5perp(u) * etaprime(u, q2)) / (2.0 * f_V_perp * power_of<5>(u)) + (f_V_para * m_b2 * power_of<5>(m_V) * eta(u, q2) * psi5perpprime(u)) / (4.0 * f_V_perp * power_of<3>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 + q2 - m_V2 * u2) * eta(u, q2) * psi5perpprime(u)) / (2.0 * f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 + q2 - m_V2 * u2) * etaprime(u, q2) * psi5perpprime(u)) / (8.0 * f_V_perp * power_of<4>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 + q2 - m_V2 * u2) * eta(u, q2) * psi5perpsecond(u)) / (16.0 * f_V_perp * power_of<4>(u)));

            return result;
        }

        const std::array<double, 4> I_T5_eta(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            std::array<double, 4> result{{0, 0, 0, 0}};

            const double phi2para_u = phi2para(u), phi2perp_u = phi2perp(u), phi3perp_u = phi3perp(u), phi4perp_u = phi4perp(u), phi5perp_u = phi5perp(u), psi5perp_u = psi5perp(u);
            const double Barphi2perp_u = Barphi2perp(u), Barphi3perp_u = Barphi3perp(u), Barphi4para_u = Barphi4para(u), Barpsi4perp_u = Barpsi4perp(u), Barphi5perp_u = Barphi5perp(u);
            const double BarBarphi2perp_u = BarBarphi2perp(u), BarBarphi3para_u = BarBarphi3para(u), BarBarpsi4perp_u = BarBarpsi4perp(u);

            // I_T5_eta_0
            result[0] = (-4.0 * f_V_perp * m_b * m_V2 * BarBarphi2perp_u + 8.0 * f_V_perp * m_b * m_V2 * BarBarphi3para_u - 4.0 * f_V_perp * m_b * m_V2 * BarBarpsi4perp_u - 4.0 * f_V_perp * m_b * m_V2 * u * Barphi2perp_u - 4.0 * f_V_para * m_b2 * m_V * Barphi3perp_u + 4.0 * f_V_para * m_V * q2 * Barphi3perp_u + 4.0 * f_V_para * power_of<3>(m_V) * u2 * Barphi3perp_u + f_V_para * power_of<3>(m_V) * Barphi4para_u - f_V_para * power_of<3>(m_V) * Barphi5perp_u + 4.0 * f_V_perp * m_b * m_V2 * u * Barpsi4perp_u + 4.0 * f_V_para * m_b2 * m_V * phi2para_u - 4.0 * f_V_para * m_V * q2 * phi2para_u - 4.0 * f_V_para * power_of<3>(m_V) * u2 * phi2para_u + 4.0 * f_V_perp * power_of<3>(m_b) * phi2perp_u - 4.0 * f_V_perp * m_b * q2 * phi2perp_u - 4.0 * f_V_perp * m_b * m_V2 * u2 * phi2perp_u + 4.0 * f_V_para * m_b2 * m_V * u * phi3perp_u + 4.0 * f_V_para * m_V * q2 * u * phi3perp_u - 4.0 * f_V_para * power_of<3>(m_V) * power_of<3>(u) * phi3perp_u + f_V_para * power_of<3>(m_V) * u * phi5perp_u) / (8.0 * f_V_perp * u2);
            // I_T5_eta_1
            result[1] = (m_b * m_V2 * (m_b2 - q2 - m_V2 * u2) * BarBarphi2perp_u) / (2.0 * power_of<3>(u)) + (m_b * m_V2 * (-m_b2 + q2 + m_V2 * u2) * BarBarphi3para_u) / power_of<3>(u) + (m_b * m_V2 * (m_b2 - q2 - m_V2 * u2) * BarBarpsi4perp_u) / (2.0 * power_of<3>(u)) + (m_b * m_V2 * (m_b2 + q2 - m_V2 * u2) * Barphi2perp_u) / (2.0 * u2) + (f_V_para * power_of<3>(m_V) * (m_b2 + q2 + m_V2 * u2) * Barphi4para_u) / (8.0 * f_V_perp * power_of<3>(u)) - (f_V_para * power_of<3>(m_V) * (m_b2 + q2 + m_V2 * u2) * Barphi5perp_u) / (8.0 * f_V_perp * power_of<3>(u)) - (m_b * m_V2 * (m_b2 + q2 - m_V2 * u2) * Barpsi4perp_u) / (2.0 * u2) + (power_of<3>(m_b) * m_V2 * phi4perp_u) / (4.0 * power_of<3>(u)) + (f_V_para * power_of<3>(m_V) * (m_b2 - q2 + m_V2 * u2) * phi5perp_u) / (8.0 * f_V_perp * u2) + (f_V_para * m_V * (q2 * q2 + (m_b2 - m_V2 * u2) * (m_b2 - m_V2 * u2) - 2.0 * q2 * (m_b2 + m_V2 * u2)) * psi3perp(u)) / (8.0 * f_V_perp * power_of<3>(u)) - (f_V_para * power_of<3>(m_V) * (q2 + m_V2 * u2) * psi5perp_u) / (16.0 * f_V_perp * power_of<3>(u));
            // I_T5_eta_2
            result[2] = (f_V_para * m_b2 * power_of<3>(m_V) * (-m_b2 + q2 + m_V2 * u2) * Barphi4para_u) / (8.0 * f_V_perp * power_of<4>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * u2) * Barphi5perp_u) / (8.0 * f_V_perp * power_of<4>(u)) + (power_of<3>(m_b) * m_V2 * (-m_b2 + q2 + m_V2 * u2) * phi4perp_u) / (8.0 * power_of<4>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 + q2 - m_V2 * u2) * phi5perp_u) / (8.0 * f_V_perp * power_of<3>(u)) - (f_V_para * power_of<3>(m_V) * (-3.0 * power_of<4>(m_b) + q2 * q2 + 2.0 * m_b2 * m_V2 * u2 + power_of<4>(m_V) * power_of<4>(u) + 2.0 * q2 * (m_b2 - m_V2 * u2)) * psi5perp_u) / (32.0 * f_V_perp * power_of<4>(u));
            // I_T5_eta_3
            result[3] = -1.0 / 32.0 * (f_V_para * m_b2 * power_of<3>(m_V) * (q2 * q2 + (m_b2 - m_V2 * u2) * (m_b2 - m_V2 * u2) - 2.0 * q2 * (m_b2 + m_V2 * u2)) * psi5perp_u) / (f_V_perp * power_of<5>(u));

            return result;
        }

        const std::array<std::array<double, 4>, 4> D_I_T5_eta(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            const std::array<double, 4> I_T5_eta = this->I_T5_eta(u, q2);
            // result[n][m] = D_eta^m[I_n]
            std::array<std::array<double, 4>, 4> result{{{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}};

            // I_T5_eta_0
            result[0][0] = I_T5_eta[0];
            // I_T5_eta_1
            result[1][0] = I_T5_eta[1];
            // I_T5_eta_2
            result[2][0] = I_T5_eta[2];
            result[2][1] = I_T5_eta[2] * etaprime(u, q2) + eta(u, q2) * ((f_V_para * m_b2 * power_of<5>(m_V) * Barphi4para(u)) / (4.0 * f_V_perp * power_of<3>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (-m_b2 + q2 + m_V2 * u2) * Barphi4para(u)) / (2.0 * f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<5>(m_V) * Barphi5perp(u)) / (4.0 * f_V_perp * power_of<3>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * u2) * Barphi5perp(u)) / (2.0 * f_V_perp * power_of<5>(u)) + (power_of<3>(m_b) * power_of<4>(m_V) * phi4perp(u)) / (4.0 * power_of<3>(u)) - (power_of<3>(m_b) * m_V2 * (-m_b2 + q2 + m_V2 * u2) * phi4perp(u)) / (2.0 * power_of<5>(u)) + (f_V_para * m_b2 * power_of<5>(m_V) * phi5perp(u)) / (4.0 * f_V_perp * u2) + (3.0 * f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 + q2 - m_V2 * u2) * phi5perp(u)) / (8.0 * f_V_perp * power_of<4>(u)) - (f_V_para * power_of<3>(m_V) * (4.0 * m_b2 * m_V2 * u - 4.0 * m_V2 * q2 * u + 4.0 * power_of<4>(m_V) * power_of<3>(u)) * psi5perp(u)) / (32.0 * f_V_perp * power_of<4>(u)) + (f_V_para * power_of<3>(m_V) * (-3.0 * power_of<4>(m_b) + q2 * q2 + 2.0 * m_b2 * m_V2 * u2 + power_of<4>(m_V) * power_of<4>(u) + 2.0 * q2 * (m_b2 - m_V2 * u2)) * psi5perp(u)) / (8.0 * f_V_perp * power_of<5>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (-m_b2 + q2 + m_V2 * u2) * phi4para(u)) / (8.0 * f_V_perp * power_of<4>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * u2) * phi5perp(u)) / (8.0 * f_V_perp * power_of<4>(u)) + (power_of<3>(m_b) * m_V2 * (-m_b2 + q2 + m_V2 * u2) * phi4perpprime(u)) / (8.0 * power_of<4>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 + q2 - m_V2 * u2) * phi5perpprime(u)) / (8.0 * f_V_perp * power_of<3>(u)) - (f_V_para * power_of<3>(m_V) * (-3.0 * power_of<4>(m_b) + q2 * q2 + 2.0 * m_b2 * m_V2 * u2 + power_of<4>(m_V) * power_of<4>(u) + 2.0 * q2 * (m_b2 - m_V2 * u2)) * psi5perpprime(u)) / (32.0 * f_V_perp * power_of<4>(u)));
            // I_T5_eta_3
            result[3][0] = I_T5_eta[3];
            result[3][1] = I_T5_eta[3] * etaprime(u, q2) - 1.0 / 32.0 * (f_V_para * m_b2 * power_of<3>(m_V) * (-4.0 * m_V2 * q2 * u - 4.0 * m_V2 * u * (m_b2 - m_V2 * u2)) * eta(u, q2) * psi5perp(u)) / (f_V_perp * power_of<5>(u)) + (5 * f_V_para * m_b2 * power_of<3>(m_V) * (q2 * q2 + (m_b2 - m_V2 * u2) * (m_b2 - m_V2 * u2) - 2.0 * q2 * (m_b2 + m_V2 * u2)) * eta(u, q2) * psi5perp(u)) / (32.0 * f_V_perp * power_of<6>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (q2 * q2 + (m_b2 - m_V2 * u2) * (m_b2 - m_V2 * u2) - 2.0 * q2 * (m_b2 + m_V2 * u2)) * eta(u, q2) * psi5perpprime(u)) / (32.0 * f_V_perp * power_of<5>(u));
            result[3][2] = result[3][1] * etaprime(u, q2) + eta(u, q2) * (I_T5_eta[3] * etasecond(u, q2) - 1.0 / 32.0 * (f_V_para * m_b2 * power_of<3>(m_V) * (-4.0 * m_V2 * q2 + 8.0 * power_of<4>(m_V) * u2 - 4.0 * m_V2 * (m_b2 - m_V2 * u2)) * eta(u, q2) * psi5perp(u)) / (f_V_perp * power_of<5>(u)) + (5 * f_V_para * m_b2 * power_of<3>(m_V) * (-4.0 * m_V2 * q2 * u - 4.0 * m_V2 * u * (m_b2 - m_V2 * u2)) * eta(u, q2) * psi5perp(u)) / (16.0 * f_V_perp * power_of<6>(u)) - (15 * f_V_para * m_b2 * power_of<3>(m_V) * (q2 * q2 + (m_b2 - m_V2 * u2) * (m_b2 - m_V2 * u2) - 2.0 * q2 * (m_b2 + m_V2 * u2)) * eta(u, q2) * psi5perp(u)) / (16.0 * f_V_perp * power_of<7>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (-4.0 * m_V2 * q2 * u - 4.0 * m_V2 * u * (m_b2 - m_V2 * u2)) * psi5perp(u) * etaprime(u, q2)) / (16.0 * f_V_perp * power_of<5>(u)) + (5 * f_V_para * m_b2 * power_of<3>(m_V) * (q2 * q2 + (m_b2 - m_V2 * u2) * (m_b2 - m_V2 * u2) - 2.0 * q2 * (m_b2 + m_V2 * u2)) * psi5perp(u) * etaprime(u, q2)) / (16.0 * f_V_perp * power_of<6>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (-4.0 * m_V2 * q2 * u - 4.0 * m_V2 * u * (m_b2 - m_V2 * u2)) * eta(u, q2) * psi5perpprime(u)) / (16.0 * f_V_perp * power_of<5>(u)) + (5 * f_V_para * m_b2 * power_of<3>(m_V) * (q2 * q2 + (m_b2 - m_V2 * u2) * (m_b2 - m_V2 * u2) - 2.0 * q2 * (m_b2 + m_V2 * u2)) * eta(u, q2) * psi5perpprime(u)) / (16.0 * f_V_perp * power_of<6>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (q2 * q2 + (m_b2 - m_V2 * u2) * (m_b2 - m_V2 * u2) - 2.0 * q2 * (m_b2 + m_V2 * u2)) * etaprime(u, q2) * psi5perpprime(u)) / (16.0 * f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (q2 * q2 + (m_b2 - m_V2 * u2) * (m_b2 - m_V2 * u2) - 2.0 * q2 * (m_b2 + m_V2 * u2)) * eta(u, q2) * psi5perpsecond(u)) / (32.0 * f_V_perp * power_of<5>(u)));

            return result;
        }

        const std::array<double, 4> I_T5_q(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            std::array<double, 4> result{{0, 0, 0, 0}};

            const double phi2para_u = phi2para(u), psi3perp_u = psi3perp(u), phi5perp_u = phi5perp(u), psi5perp_u = psi5perp(u);
            const double Barphi3perp_u = Barphi3perp(u), Barphi4para_u = Barphi4para(u), Barphi5perp_u = Barphi5perp(u);
            const double BarBarphi2perp_u = BarBarphi2perp(u), BarBarphi3para_u = BarBarphi3para(u), BarBarpsi4perp_u = BarBarpsi4perp(u);

            // I_T5_q_0
            result[0] = (-4.0 * f_V_para * m_V * (-1.0 + u) * Barphi3perp_u + 4.0 * f_V_para * m_V * (-1.0 + u) * phi2para_u + 4.0 * f_V_perp * m_b * u * phi2perp(u) + f_V_para * m_V * u * (4.0 * (-1.0 + u) * phi3perp(u) + psi3perp_u)) / (4.0 * f_V_perp * u2);
            // I_T5_q_1
            result[1] = (m_b * m_V2 * (-2.0 + u) * BarBarphi2perp_u) / power_of<3>(u) - (2.0 * m_b * m_V2 * (-2.0 + u) * BarBarphi3para_u) / power_of<3>(u) + (m_b * m_V2 * (-2.0 + u) * BarBarpsi4perp_u) / power_of<3>(u) + (m_b * m_V2 * (-1.0 + u) * Barphi2perp(u)) / u2 + (f_V_para * m_V * (-m_b2 + q2 * (1 - 2.0 * u) + m_V2 * u2) * Barphi3perp_u) / (f_V_perp * power_of<3>(u)) - (f_V_para * power_of<3>(m_V) * (-2.0 + u) * Barphi4para_u) / (4.0 * f_V_perp * power_of<3>(u)) + (f_V_para * power_of<3>(m_V) * (-2.0 + u) * Barphi5perp_u) / (4.0 * f_V_perp * power_of<3>(u)) - (m_b * m_V2 * (-1.0 + u) * Barpsi4perp(u)) / u2 + (f_V_para * m_V * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * phi2para_u) / (f_V_perp * power_of<3>(u)) - (f_V_para * power_of<3>(m_V) * (-1.0 + u) * phi5perp_u) / (4.0 * f_V_perp * u2) - (f_V_para * m_V * (-m_b2 + q2 + m_V2 * (-2.0 + u) * u) * psi3perp_u) / (4.0 * f_V_perp * u2);
            // I_T5_q_2
            result[2] = (m_b * m_V2 * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * BarBarphi2perp_u) / power_of<4>(u) - (2.0 * m_b * m_V2 * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * BarBarphi3para_u) / power_of<4>(u) + (m_b * m_V2 * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * BarBarpsi4perp_u) / power_of<4>(u) + (f_V_para * power_of<3>(m_V) * (q2 * (1 - 2.0 * u) - m_b2 * (-2.0 + u) + m_V2 * u2) * Barphi4para_u) / (4.0 * f_V_perp * power_of<4>(u)) - (f_V_para * power_of<3>(m_V) * (q2 * (1 - 2.0 * u) - m_b2 * (-2.0 + u) + m_V2 * u2) * Barphi5perp_u) / (4.0 * f_V_perp * power_of<4>(u)) - (power_of<3>(m_b) * m_V2 * phi4perp(u)) / (4.0 * power_of<3>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (-1.0 + u) * phi5perp_u) / (4.0 * f_V_perp * power_of<3>(u)) + (f_V_para * power_of<3>(m_V) * (q2 + m_V2 * (-2.0 + u) * u) * psi5perp_u) / (16.0 * f_V_perp * power_of<3>(u));
            // I_T5_q_3
            result[3] = -0.25 * (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * Barphi4para_u) / (f_V_perp * power_of<5>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * Barphi5perp_u) / (4.0 * f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * (-2.0 + u) * u) * psi5perp_u) / (16.0 * f_V_perp * power_of<4>(u));

            return result;
        }

        const std::array<std::array<double, 4>, 4> D_I_T5_q(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            const std::array<double, 4> I_T5_q = this->I_T5_q(u, q2);
            // result[n][m] = D_eta^m[I_n]
            std::array<std::array<double, 4>, 4> result{{{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}};

            // I_T5_q_0
            result[0][0] = I_T5_q[0];
            // I_T5_q_1
            result[1][0] = I_T5_q[1];
            // I_T5_q_2
            result[2][0] = I_T5_q[2];
            result[2][1] = I_T5_q[2] * etaprime(u, q2) + eta(u, q2) * ((m_b * m_V2 * (2.0 * q2 - 2.0 * m_V2 * u) * BarBarphi2perp(u)) / power_of<4>(u) - (4.0 * m_b * m_V2 * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * BarBarphi2perp(u)) / power_of<5>(u) - (2.0 * m_b * m_V2 * (2.0 * q2 - 2.0 * m_V2 * u) * BarBarphi3para(u)) / power_of<4>(u) + (8.0 * m_b * m_V2 * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * BarBarphi3para(u)) / power_of<5>(u) + (m_b * m_V2 * (2.0 * q2 - 2.0 * m_V2 * u) * BarBarpsi4perp(u)) / power_of<4>(u) - (4.0 * m_b * m_V2 * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * BarBarpsi4perp(u)) / power_of<5>(u) + (f_V_para * power_of<3>(m_V) * (-m_b2 - 2.0 * q2 + 2.0 * m_V2 * u) * Barphi4para(u)) / (4.0 * f_V_perp * power_of<4>(u)) - (f_V_para * power_of<3>(m_V) * (q2 * (1 - 2.0 * u) - m_b2 * (-2.0 + u) + m_V2 * u2) * Barphi4para(u)) / (f_V_perp * power_of<5>(u)) - (f_V_para * power_of<3>(m_V) * (-m_b2 - 2.0 * q2 + 2.0 * m_V2 * u) * Barphi5perp(u)) / (4.0 * f_V_perp * power_of<4>(u)) + (f_V_para * power_of<3>(m_V) * (q2 * (1 - 2.0 * u) - m_b2 * (-2.0 + u) + m_V2 * u2) * Barphi5perp(u)) / (f_V_perp * power_of<5>(u)) + (3.0 * power_of<3>(m_b) * m_V2 * phi4perp(u)) / (4.0 * power_of<4>(u)) + (3.0 * f_V_para * m_b2 * power_of<3>(m_V) * (-1.0 + u) * phi5perp(u)) / (4.0 * f_V_perp * power_of<4>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * phi5perp(u)) / (4.0 * f_V_perp * power_of<3>(u)) + (f_V_para * power_of<3>(m_V) * (m_V2 * (-2.0 + u) + m_V2 * u) * psi5perp(u)) / (16.0 * f_V_perp * power_of<3>(u)) - (3.0 * f_V_para * power_of<3>(m_V) * (q2 + m_V2 * (-2.0 + u) * u) * psi5perp(u)) / (16.0 * f_V_perp * power_of<4>(u)) + (m_b * m_V2 * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * Barphi2perp(u)) / power_of<4>(u) - (2.0 * m_b * m_V2 * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * Barphi3para(u)) / power_of<4>(u) + (m_b * m_V2 * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * Barpsi4perp(u)) / power_of<4>(u) + (f_V_para * power_of<3>(m_V) * (q2 * (1 - 2.0 * u) - m_b2 * (-2.0 + u) + m_V2 * u2) * phi4para(u)) / (4.0 * f_V_perp * power_of<4>(u)) - (f_V_para * power_of<3>(m_V) * (q2 * (1 - 2.0 * u) - m_b2 * (-2.0 + u) + m_V2 * u2) * phi5perp(u)) / (4.0 * f_V_perp * power_of<4>(u)) - (power_of<3>(m_b) * m_V2 * phi4perpprime(u)) / (4.0 * power_of<3>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (-1.0 + u) * phi5perpprime(u)) / (4.0 * f_V_perp * power_of<3>(u)) + (f_V_para * power_of<3>(m_V) * (q2 + m_V2 * (-2.0 + u) * u) * psi5perpprime(u)) / (16.0 * f_V_perp * power_of<3>(u)));
            // I_T5_q_3
            result[3][0] = I_T5_q[3];
            result[3][1] = I_T5_q[3] * etaprime(u, q2) + eta(u, q2) * (-0.25 * (f_V_para * m_b2 * power_of<3>(m_V) * (2.0 * q2 - 2.0 * m_V2 * u) * Barphi4para(u)) / (f_V_perp * power_of<5>(u)) + (5 * f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * Barphi4para(u)) / (4.0 * f_V_perp * power_of<6>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (2.0 * q2 - 2.0 * m_V2 * u) * Barphi5perp(u)) / (4.0 * f_V_perp * power_of<5>(u)) - (5 * f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * Barphi5perp(u)) / (4.0 * f_V_perp * power_of<6>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (-(m_V2 * (-2.0 + u)) - m_V2 * u) * psi5perp(u)) / (16.0 * f_V_perp * power_of<4>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * (-2.0 + u) * u) * psi5perp(u)) / (4.0 * f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * phi4para(u)) / (4.0 * f_V_perp * power_of<5>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * phi5perp(u)) / (4.0 * f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * (-2.0 + u) * u) * psi5perpprime(u)) / (16.0 * f_V_perp * power_of<4>(u)));
            result[3][2] = result[3][1] * etaprime(u, q2) + eta(u, q2) * (I_T5_q[3] * etasecond(u, q2) + 2.0 * etaprime(u, q2) * (-0.25 * (f_V_para * m_b2 * power_of<3>(m_V) * (2.0 * q2 - 2.0 * m_V2 * u) * Barphi4para(u)) / (f_V_perp * power_of<5>(u)) + (5 * f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * Barphi4para(u)) / (4.0 * f_V_perp * power_of<6>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (2.0 * q2 - 2.0 * m_V2 * u) * Barphi5perp(u)) / (4.0 * f_V_perp * power_of<5>(u)) - (5 * f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * Barphi5perp(u)) / (4.0 * f_V_perp * power_of<6>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (-(m_V2 * (-2.0 + u)) - m_V2 * u) * psi5perp(u)) / (16.0 * f_V_perp * power_of<4>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * (-2.0 + u) * u) * psi5perp(u)) / (4.0 * f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * phi4para(u)) / (4.0 * f_V_perp * power_of<5>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * phi5perp(u)) / (4.0 * f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * (-2.0 + u) * u) * psi5perpprime(u)) / (16.0 * f_V_perp * power_of<4>(u))) + eta(u, q2) * ((f_V_para * m_b2 * power_of<5>(m_V) * Barphi4para(u)) / (2.0 * f_V_perp * power_of<5>(u)) + (5 * f_V_para * m_b2 * power_of<3>(m_V) * (2.0 * q2 - 2.0 * m_V2 * u) * Barphi4para(u)) / (2.0 * f_V_perp * power_of<6>(u)) - (15 * f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * Barphi4para(u)) / (2.0 * f_V_perp * power_of<7>(u)) - (f_V_para * m_b2 * power_of<5>(m_V) * Barphi5perp(u)) / (2.0 * f_V_perp * power_of<5>(u)) - (5 * f_V_para * m_b2 * power_of<3>(m_V) * (2.0 * q2 - 2.0 * m_V2 * u) * Barphi5perp(u)) / (2.0 * f_V_perp * power_of<6>(u)) + (15 * f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * Barphi5perp(u)) / (2.0 * f_V_perp * power_of<7>(u)) + (f_V_para * m_b2 * power_of<5>(m_V) * psi5perp(u)) / (8.0 * f_V_perp * power_of<4>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (-(m_V2 * (-2.0 + u)) - m_V2 * u) * psi5perp(u)) / (2.0 * f_V_perp * power_of<5>(u)) - (5 * f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * (-2.0 + u) * u) * psi5perp(u)) / (4.0 * f_V_perp * power_of<6>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (2.0 * q2 - 2.0 * m_V2 * u) * phi4para(u)) / (2.0 * f_V_perp * power_of<5>(u)) + (5 * f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * phi4para(u)) / (2.0 * f_V_perp * power_of<6>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (2.0 * q2 - 2.0 * m_V2 * u) * phi5perp(u)) / (2.0 * f_V_perp * power_of<5>(u)) - (5 * f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * phi5perp(u)) / (2.0 * f_V_perp * power_of<6>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (-(m_V2 * (-2.0 + u)) - m_V2 * u) * psi5perpprime(u)) / (8.0 * f_V_perp * power_of<4>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * (-2.0 + u) * u) * psi5perpprime(u)) / (2.0 * f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * phi4paraprime(u)) / (4.0 * f_V_perp * power_of<5>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - m_V2 * u2 + q2 * (-1 + 2.0 * u)) * phi5perpprime(u)) / (4.0 * f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * (-2.0 + u) * u) * psi5perpsecond(u)) / (16.0 * f_V_perp * power_of<4>(u))));

            return result;
        }

        const std::array<double, 4> I_T5_p(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            std::array<double, 4> result{{0, 0, 0, 0}};

            const double phi2para_u = phi2para(u), psi3perp_u = psi3perp(u), phi5perp_u = phi5perp(u), psi5perp_u = psi5perp(u);
            const double Barphi3perp_u = Barphi3perp(u), Barphi4para_u = Barphi4para(u), Barphi5perp_u = Barphi5perp(u);
            const double BarBarphi2perp_u = BarBarphi2perp(u), BarBarphi3para_u = BarBarphi3para(u), BarBarpsi4perp_u = BarBarpsi4perp(u);

            // I_T5_p_0
            result[0] = -0.25 * (-4.0 * f_V_para * m_V * Barphi3perp_u + 4.0 * f_V_para * m_V * phi2para_u + 4.0 * f_V_perp * m_b * phi2perp(u) + 4.0 * f_V_para * m_V * u * phi3perp(u) + f_V_para * m_V * psi3perp_u) / (f_V_perp * u);
            // I_T5_p_1
            result[1] = -((m_b * m_V2 * BarBarphi2perp_u) / u2) + (2.0 * m_b * m_V2 * BarBarphi3para_u) / u2 - (m_b * m_V2 * BarBarpsi4perp_u) / u2 - (m_b * m_V2 * Barphi2perp(u)) / u + (2.0 * f_V_para * m_V * q2 * Barphi3perp_u) / (f_V_perp * u2) + (f_V_para * power_of<3>(m_V) * Barphi4para_u) / (4.0 * f_V_perp * u2) - (f_V_para * power_of<3>(m_V) * Barphi5perp_u) / (4.0 * f_V_perp * u2) + (m_b * m_V2 * Barpsi4perp(u)) / u - (2.0 * f_V_para * m_V * q2 * phi2para_u) / (f_V_perp * u2) + (f_V_para * power_of<3>(m_V) * phi5perp_u) / (4.0 * f_V_perp * u) + (f_V_para * m_V * (-m_b2 + q2 + m_V2 * u2) * psi3perp_u) / (4.0 * f_V_perp * u2);
            // I_T5_p_2
            result[2] = (-2.0 * m_b * m_V2 * q2 * BarBarphi2perp_u) / power_of<3>(u) + (4.0 * m_b * m_V2 * q2 * BarBarphi3para_u) / power_of<3>(u) - (2.0 * m_b * m_V2 * q2 * BarBarpsi4perp_u) / power_of<3>(u) + (f_V_para * power_of<3>(m_V) * (m_b2 + 2.0 * q2) * Barphi4para_u) / (4.0 * f_V_perp * power_of<3>(u)) - (f_V_para * power_of<3>(m_V) * (m_b2 + 2.0 * q2) * Barphi5perp_u) / (4.0 * f_V_perp * power_of<3>(u)) + (power_of<3>(m_b) * m_V2 * phi4perp(u)) / (4.0 * power_of<3>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * phi5perp_u) / (4.0 * f_V_perp * u2) - (f_V_para * power_of<3>(m_V) * (q2 + m_V2 * u2) * psi5perp_u) / (16.0 * f_V_perp * power_of<3>(u));
            // I_T5_p_3
            result[3] = (f_V_para * m_b2 * power_of<3>(m_V) * q2 * Barphi4para_u) / (2.0 * f_V_perp * power_of<4>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * q2 * Barphi5perp_u) / (2.0 * f_V_perp * power_of<4>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * u2) * psi5perp_u) / (16.0 * f_V_perp * power_of<4>(u));

            return result;
        }

        const std::array<std::array<double, 4>, 4> D_I_T5_p(const double & u, const double & q2) const
        {
            const double m_b = this->m_b();
            const double m_V2 = m_V * m_V, m_b2 = m_b * m_b, u2 = u * u;
            const double f_V_perp = this->f_V_perp(), f_V_para = this->f_V_para();
            const std::array<double, 4> I_T5_p = this->I_T5_p(u, q2);
            // result[n][m] = D_eta^m[I_n]
            std::array<std::array<double, 4>, 4> result{{{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}}};

            // I_T5_p_0
            result[0][0] = I_T5_p[0];
            // I_T5_p_1
            result[1][0] = I_T5_p[1];
            // I_T5_p_2
            result[2][0] = I_T5_p[2];
            result[2][1] = I_T5_p[2] * etaprime(u, q2) + eta(u, q2) * ((6 * m_b * m_V2 * q2 * BarBarphi2perp(u)) / power_of<4>(u) - (12 * m_b * m_V2 * q2 * BarBarphi3para(u)) / power_of<4>(u) + (6 * m_b * m_V2 * q2 * BarBarpsi4perp(u)) / power_of<4>(u) - (3.0 * f_V_para * power_of<3>(m_V) * (m_b2 + 2.0 * q2) * Barphi4para(u)) / (4.0 * f_V_perp * power_of<4>(u)) + (3.0 * f_V_para * power_of<3>(m_V) * (m_b2 + 2.0 * q2) * Barphi5perp(u)) / (4.0 * f_V_perp * power_of<4>(u)) - (3.0 * power_of<3>(m_b) * m_V2 * phi4perp(u)) / (4.0 * power_of<4>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * phi5perp(u)) / (2.0 * f_V_perp * power_of<3>(u)) - (f_V_para * power_of<5>(m_V) * psi5perp(u)) / (8.0 * f_V_perp * u2) + (3.0 * f_V_para * power_of<3>(m_V) * (q2 + m_V2 * u2) * psi5perp(u)) / (16.0 * f_V_perp * power_of<4>(u)) - (2.0 * m_b * m_V2 * q2 * Barphi2perp(u)) / power_of<3>(u) + (4.0 * m_b * m_V2 * q2 * Barphi3para(u)) / power_of<3>(u) - (2.0 * m_b * m_V2 * q2 * Barpsi4perp(u)) / power_of<3>(u) + (f_V_para * power_of<3>(m_V) * (m_b2 + 2.0 * q2) * phi4para(u)) / (4.0 * f_V_perp * power_of<3>(u)) - (f_V_para * power_of<3>(m_V) * (m_b2 + 2.0 * q2) * phi5perp(u)) / (4.0 * f_V_perp * power_of<3>(u)) + (power_of<3>(m_b) * m_V2 * phi4perpprime(u)) / (4.0 * power_of<3>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * phi5perpprime(u)) / (4.0 * f_V_perp * u2) - (f_V_para * power_of<3>(m_V) * (q2 + m_V2 * u2) * psi5perpprime(u)) / (16.0 * f_V_perp * power_of<3>(u)));
            // I_T5_p_3
            result[3][0] = I_T5_p[3];
            result[3][1] = I_T5_p[3] * etaprime(u, q2) + eta(u, q2) * ((-2.0 * f_V_para * m_b2 * power_of<3>(m_V) * q2 * Barphi4para(u)) / (f_V_perp * power_of<5>(u)) + (2.0 * f_V_para * m_b2 * power_of<3>(m_V) * q2 * Barphi5perp(u)) / (f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<5>(m_V) * psi5perp(u)) / (8.0 * f_V_perp * power_of<3>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * u2) * psi5perp(u)) / (4.0 * f_V_perp * power_of<5>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * q2 * phi4para(u)) / (2.0 * f_V_perp * power_of<4>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * q2 * phi5perp(u)) / (2.0 * f_V_perp * power_of<4>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * u2) * psi5perpprime(u)) / (16.0 * f_V_perp * power_of<4>(u)));
            result[3][2] = result[3][1] * etaprime(u, q2) + eta(u, q2) * (I_T5_p[3] * etasecond(u, q2) + 2.0 * etaprime(u, q2) * ((-2.0 * f_V_para * m_b2 * power_of<3>(m_V) * q2 * Barphi4para(u)) / (f_V_perp * power_of<5>(u)) + (2.0 * f_V_para * m_b2 * power_of<3>(m_V) * q2 * Barphi5perp(u)) / (f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<5>(m_V) * psi5perp(u)) / (8.0 * f_V_perp * power_of<3>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * u2) * psi5perp(u)) / (4.0 * f_V_perp * power_of<5>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * q2 * phi4para(u)) / (2.0 * f_V_perp * power_of<4>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * q2 * phi5perp(u)) / (2.0 * f_V_perp * power_of<4>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * u2) * psi5perpprime(u)) / (16.0 * f_V_perp * power_of<4>(u))) + eta(u, q2) * ((10 * f_V_para * m_b2 * power_of<3>(m_V) * q2 * Barphi4para(u)) / (f_V_perp * power_of<6>(u)) - (10 * f_V_para * m_b2 * power_of<3>(m_V) * q2 * Barphi5perp(u)) / (f_V_perp * power_of<6>(u)) + (7 * f_V_para * m_b2 * power_of<5>(m_V) * psi5perp(u)) / (8.0 * f_V_perp * power_of<4>(u)) + (5 * f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * u2) * psi5perp(u)) / (4.0 * f_V_perp * power_of<6>(u)) - (4.0 * f_V_para * m_b2 * power_of<3>(m_V) * q2 * phi4para(u)) / (f_V_perp * power_of<5>(u)) + (4.0 * f_V_para * m_b2 * power_of<3>(m_V) * q2 * phi5perp(u)) / (f_V_perp * power_of<5>(u)) - (f_V_para * m_b2 * power_of<5>(m_V) * psi5perpprime(u)) / (4.0 * f_V_perp * power_of<3>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * u2) * psi5perpprime(u)) / (2.0 * f_V_perp * power_of<5>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * q2 * phi4paraprime(u)) / (2.0 * f_V_perp * power_of<4>(u)) - (f_V_para * m_b2 * power_of<3>(m_V) * q2 * phi5perpprime(u)) / (2.0 * f_V_perp * power_of<4>(u)) + (f_V_para * m_b2 * power_of<3>(m_V) * (m_b2 - q2 - m_V2 * u2) * psi5perpsecond(u)) / (16.0 * f_V_perp * power_of<4>(u))));

            return result;
        }
        // }}}

        /* A1 : form factor and daughter sum rule */
        // {{{
        double a_1(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_A1() + q2 * s0_1_A1();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_eta_q2 = I_eta(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = exp(-s(u, q2) / M2()) * I_eta_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_eta(u0, q2);

            double P_eta = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_eta -= exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_eta += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_eta *= f_V_perp() * m_b();

            return 1.0 / (m_B2 * f_B() * (m_B() + m_V())) * exp(m_B2 / M2()) * P_eta;
        }

        double a_1_dsr(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_A1() + q2 * s0_1_A1();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_eta_q2 = I_eta(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = (n * M2() - s(u, q2)) * exp(-s(u, q2) / M2()) * I_eta_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_eta(u0, q2);

            double P_eta_derivative = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_eta_derivative -= (n * M2() - s0) * exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_eta_derivative += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_eta_derivative *= f_V_perp() * m_b();

            return -1.0 / (m_B2 * f_B() * (m_B() + m_V())) * exp(m_B2 / M2()) * P_eta_derivative / a_1(q2);
        }
        // }}}

        /* A2 : form factor and daughter sum rule  */
        // {{{
        double a_2(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_A2() + q2 * s0_1_A2();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_p_q2 = I_p(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = exp(-s(u, q2) / M2()) * I_p_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_p(u0, q2);

            double P_p = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_p -= exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_p += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_p *= f_V_perp() * m_b();

            return -0.5 * (m_B() + m_V()) / m_B2 / f_B() * exp(m_B2 / M2()) * P_p;
        }

        double a_2_dsr(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_A2() + q2 * s0_1_A2();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_p_q2 = I_p(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = (n * M2() - s(u, q2)) * exp(-s(u, q2) / M2()) * I_p_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_p(u0, q2);

            double P_p_derivative = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_p_derivative -= (n * M2() - s0) * exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_p_derivative += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_p_derivative *= f_V_perp() * m_b();

            return 0.5 * (m_B() + m_V()) / m_B2 / f_B() * exp(m_B2 / M2()) * P_p_derivative / a_2(q2);
        }
        // }}}

        /* A30 : form factor and daughter sum rule  */
        // {{{
        double a_30(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_A30() + q2 * s0_1_A30();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_q_q2 = I_q(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = exp(-s(u, q2) / M2()) * I_q_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_q(u0, q2);

            double P_q = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_q -= exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_q += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_q *= f_V_perp() * m_b();

            return - (m_B() + m_V()) / m_B2 / f_B() * exp(m_B2 / M2()) * P_q;
        }

        double a_30_dsr(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_A30() + q2 * s0_1_A30();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_q_q2 = I_q(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = (n * M2() - s(u, q2)) * exp(-s(u, q2) / M2()) * I_q_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_q(u0, q2);

            double P_q_derivative = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_q_derivative -= (n * M2() - s0) * exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_q_derivative += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_q_derivative *= f_V_perp() * m_b();

            return (m_B() + m_V()) / m_B2 / f_B() * exp(m_B2 / M2()) * P_q_derivative / a_30(q2);
        }
        // }}}


        /* V : form factor and daughter sum rule */
        // {{{
        double v(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_V() + q2 * s0_1_V();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_perp_q2 = I_perp(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = exp(-s(u, q2) / M2()) * I_perp_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_perp(u0, q2);

            double P_perp = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_perp -= exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_perp += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_perp *= f_V_perp() * m_b();

            return 0.5 * (m_B() + m_V()) / m_B2 / f_B() * exp(m_B2 / M2()) * P_perp;
        }

        double v_dsr(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_V() + q2 * s0_1_V();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_perp_q2 = I_perp(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = (n * M2() - s(u, q2)) * exp(-s(u, q2) / M2()) * I_perp_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_perp(u0, q2);

            double P_perp_derivative = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_perp_derivative -= (n * M2() - s0) * exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_perp_derivative += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_perp_derivative *= f_V_perp() * m_b();

            return -0.5 * (m_B() + m_V()) / m_B2 / f_B() * exp(m_B2 / M2()) * P_perp_derivative / v(q2);
        }
        // }}}

        /* T1 : form factor and daughter sum rule */
        // {{{
        double t_1(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_T1() + q2 * s0_1_T1();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_T_perp_q2 = I_T_perp(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = exp(-s(u, q2) / M2()) * I_T_perp_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_T_perp(u0, q2);

            double P_T_perp = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_T_perp -= exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_T_perp += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_T_perp *= f_V_perp() * m_b();

            return 0.5 / m_B2 / f_B() * exp(m_B2 / M2()) * P_T_perp;
        }

        double t_1_dsr(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_T1() + q2 * s0_1_T1();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_T_perp_q2 = I_T_perp(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = (n * M2() - s(u, q2)) * exp(-s(u, q2) / M2()) * I_T_perp_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_T_perp(u0, q2);

            double P_T_perp_derivative = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_T_perp_derivative -= (n * M2() - s0) * exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_T_perp_derivative += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_T_perp_derivative *= f_V_perp() * m_b();

            return -0.5 / m_B2 / f_B() * exp(m_B2 / M2()) * P_T_perp_derivative / t_1(q2);
        }
        // }}}

        /* T23A : form factor and daughter sum rule */
        // {{{
        double t_23A(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_T23A() + q2 * s0_1_T23A();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_T5_p_q2 = I_T5_p(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = exp(-s(u, q2) / M2()) * I_T5_p_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_T5_p(u0, q2);

            double P_T5_p = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_T5_p -= exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_T5_p += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_T5_p *= f_V_perp() * m_b();

            return -0.5 / m_B2 / f_B() * exp(m_B2 / M2()) * P_T5_p;
        }

        double t_23A_dsr(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_T23A() + q2 * s0_1_T23A();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_T5_p_q2 = I_T5_p(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = (n * M2() - s(u, q2)) * exp(-s(u, q2) / M2()) * I_T5_p_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_T5_p(u0, q2);

            double P_T5_p_derivative = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_T5_p_derivative -= (n * M2() - s0) * exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_T5_p_derivative += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_T5_p_derivative *= f_V_perp() * m_b();

            return 0.5 / m_B2 / f_B() * exp(m_B2 / M2()) * P_T5_p_derivative / t_23A(q2);
        }
        // }}}

        /* T23B : form factor and daughter sum rule */
        // {{{
        double t_23B(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_T23B() + q2 * s0_1_T23B();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_T5_q_q2 = I_T5_q(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = exp(-s(u, q2) / M2()) * I_T5_q_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_T5_q(u0, q2);

            double P_T5_q = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_T5_q -= exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_T5_q += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_T5_q *= f_V_perp() * m_b();

            return 0.5 / m_B2 / f_B() * exp(m_B2 / M2()) * P_T5_q;
        }

        double t_23B_dsr(const double & q2) const
        {
            const double m_B2 = m_B() * m_B();
            const double s0 = s0_0_T23B() + q2 * s0_1_T23B();
            const double u0 = this->u0(s0, q2);

            std::function<std::array<double, 4> (const double &)>  integrand = [this, q2](const double & u) -> std::array<double, 4>
            {
                std::array<double, 4> res;
                std::array<double, 4> I_T5_q_q2 = I_T5_q(u, q2);
                for (unsigned n = 0; n < 4; n++)
                {
                    res[n] = (n * M2() - s(u, q2)) * exp(-s(u, q2) / M2()) * I_T5_q_q2[n];
                }

                return res;
            };
            // 4D cubature integration
            std::array<double, 4> integral = integrate(integrand, u0, 1.0, cub_conf);

            const std::array<std::array<double, 4>, 4> surface = D_I_T5_q(u0, q2);

            double P_T5_q_derivative = 0.0, power_of_M2 = 1.0;

            for (unsigned n = 0; n < 4; n++)
            {
                for (unsigned i = 0; i < 3 - n; i++)
                {
                    P_T5_q_derivative -= (n * M2() - s0) * exp(-s0 / M2()) * eta(u0, q2) * surface[n + i + 1][i] / power_of_M2;
                }
                P_T5_q_derivative += integral[n] / power_of_M2;
                power_of_M2 *= M2();
            }

            P_T5_q_derivative *= f_V_perp() * m_b();

            return -0.5 / m_B2 / f_B() * exp(m_B2 / M2()) * P_T5_q_derivative / t_23B(q2);
        }
        // }}}

        /* Diagnostics */

        Diagnostics diagnostics() const
        {
            Diagnostics results;

            // Helper functions
            results.add({m_b(),                      "m_b(mu)"                                    });
            results.add({this->model->m_s_msbar(mu), "m_s(mu)"                                    });
            results.add({this->model->m_ud_msbar(mu),"m_ud(mu)"                                   });

            results.add({eta(0.3, 1.0),              "eta(u = 0.3, q2 = 1.0)"                     });
            results.add({etaprime(0.3, 1.0),         "eta'(u = 0.3, q2 = 1.0)"                    });
            results.add({etasecond(0.3, 1.0),        "eta''(u = 0.3, q2 = 1.0)"                   });

            results.add({exp(-s(0.3, 1.0) / M2()),   "exp(-s(0.3, 1.0) / M2())"                   });

            // LCDA wrappers
            results.add({this->lcdas->a1perp(mu),    "a1perp(mu)"                                 });
            results.add({phi2perp(0.3),              "phi2perp(u = 0.3)"                          });
            results.add({phi2para(0.3),              "phi2para(u = 0.3)"                          });
            results.add({phi3perp(0.3),              "phi3perp(u = 0.3)"                          });
            results.add({phi3para(0.3),              "phi3para(u = 0.3)"                          });
            results.add({phi4perp(0.3),              "phi4perp(u = 0.3)"                          });
            results.add({phi4para(0.3),              "phi4para(u = 0.3)"                          });
            results.add({psi4perp(0.3),              "psi4perp(u = 0.3)"                          });
            results.add({psi4para(0.3),              "psi4para(u = 0.3)"                          });
            results.add({Barphi2perp(0.3),           "Barphi2perp(u = 0.3)"                       });
            results.add({Barphi2para(0.3),           "Barphi2para(u = 0.3)"                       });
            results.add({BarBarphi2perp(0.3),        "BarBarphi2perp(u = 0.3)"                    });
            results.add({BarBarphi2para(0.3),        "BarBarphi2para(u = 0.3)"                    });
            results.add({Barphi3perp(0.3),           "Barphi3perp(u = 0.3)"                       });
            results.add({Barphi3para(0.3),           "Barphi3para(u = 0.3)"                       });
            results.add({BarBarphi3perp(0.3),        "BarBarphi3perp(u = 0.3)"                    });
            results.add({BarBarphi3para(0.3),        "BarBarphi3para(u = 0.3)"                    });
            results.add({Barphi4para(0.3),           "Barphi4para(u = 0.3)"                       });
            results.add({Barpsi4perp(0.3),           "Barpsi4perp(u = 0.3)"                       });
            results.add({Barpsi4para(0.3),           "Barpsi4para(u = 0.3)"                       });
            results.add({BarBarpsi4perp(0.3),        "BarBarpsi4perp(u = 0.3)"                    });
            results.add({BarBarpsi4para(0.3),        "BarBarpsi4para(u = 0.3)"                    });
            results.add({Barphi5perp(0.3),           "Barphi5perp(u = 0.3)"                       });

            // I functions
            results.add({I_perp(0.3, 1.0)[0],        "I_perp[0](u = 0.3, q2 = 1.0)"               });
            results.add({I_perp(0.3, 1.0)[1],        "I_perp[1](u = 0.3, q2 = 1.0)"               });
            results.add({I_T5_q(0.3, 1.0)[0],        "I_T5_q[0](u = 0.3, q2 = 1.0)"               });
            results.add({I_T5_q(0.3, 1.0)[1],        "I_T5_q[1](u = 0.3, q2 = 1.0)"               });
            results.add({I_T5_q(0.3, 1.0)[2],        "I_T5_q[2](u = 0.3, q2 = 1.0)"               });
            results.add({I_T5_q(0.3, 1.0)[3],        "I_T5_q[3](u = 0.3, q2 = 1.0)"               });
            results.add({D_I_T5_q(0.3, 1.0)[2][1],   "D_I_T5_q[2][1](u = 0.3, q2 = 1.0)"          });
            results.add({D_I_T5_q(0.3, 1.0)[3][1],   "D_I_T5_q[3][1](u = 0.3, q2 = 1.0)"          });
            results.add({D_I_T5_q(0.3, 1.0)[3][2],   "D_I_T5_q[3][2](u = 0.3, q2 = 1.0)"          });

            results.add({a_1(1.0),                   "a_1(q2 = 1.0)"                              });
            results.add({a_2(1.0),                   "a_2(q2 = 1.0)"                              });
            results.add({a_30(1.0),                  "a_30(q2 = 1.0)"                             });
            results.add({v(1.0),                     "v(q2 = 1.0)"                                });
            results.add({t_1(1.0),                   "t_1(q2 = 1.0)"                              });
            results.add({t_23A(1.0),                 "t_23A(q2 = 1.0)"                            });
            results.add({t_23B(1.0),                 "t_23B(q2 = 1.0)"                            });

            results.add({a_1_dsr(1.0),               "a_1_dsr(q2 = 1.0)"                          });
            results.add({a_2_dsr(1.0),               "a_2_dsr(q2 = 1.0)"                          });
            results.add({a_30_dsr(1.0),              "a_30_dsr(q2 = 1.0)"                         });
            results.add({v_dsr(1.0),                 "v_dsr(q2 = 1.0)"                            });
            results.add({t_1_dsr(1.0),               "t_1_dsr(q2 = 1.0)"                          });
            results.add({t_23A_dsr(1.0),             "t_23A_dsr(q2 = 1.0)"                        });
            results.add({t_23B_dsr(1.0),             "t_23B_dsr(q2 = 1.0)"                        });
            return results;
        }

    };

    template <typename Process_>
    const std::vector<OptionSpecification>
    Implementation<AnalyticFormFactorPToVLCSR<Process_>>::options
    {
        {"2pt-twist"_ok, "2|3|4|5"s, "2|3|4|5"s}
    };

    template <typename Process_>
    AnalyticFormFactorPToVLCSR<Process_>::AnalyticFormFactorPToVLCSR(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<AnalyticFormFactorPToVLCSR<Process_>>(new Implementation<AnalyticFormFactorPToVLCSR<Process_>>(p, o, *this))
    {
    }

    template <typename Process_>
    AnalyticFormFactorPToVLCSR<Process_>::~AnalyticFormFactorPToVLCSR()
    {
    }

    template <typename Process_>
    FormFactors<PToV> *
    AnalyticFormFactorPToVLCSR<Process_>::make(const Parameters & p, const Options & o)
    {
        return new AnalyticFormFactorPToVLCSR<Process_>(p, o);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::a_0(const double & q2) const
    {
        const double m_B = this->_imp->m_B(), m_V = this->_imp->m_V();
        const double m_B2 = m_B * m_B, m_V2 = m_V * m_V;

        return ((m_B + m_V) * this->_imp->a_1(q2) - (m_B2 + q2 - m_V2) / (m_B + m_V) * this->_imp->a_2(q2) - q2 / (m_B + m_V) * this->_imp->a_30(q2)) / (2.0 * m_V);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::a_1(const double & q2) const
    {
        return this->_imp->a_1(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::a_2(const double & q2) const
    {
        return this->_imp->a_2(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::a_12(const double & q2) const
    {
        const double m_B = this->_imp->m_B(), m_V = this->_imp->m_V();
        const double m_B2 = m_B * m_B, m_V2 = m_V * m_V;

        const double c_1 = (m_B + m_V) * (m_B2 - m_V2 - q2) / (16.0 * m_B * m_V2);
        const double c_2 = eos::lambda(m_B2, m_V2, q2) / (16.0 * m_B * m_V2 * (m_B + m_V));

        return c_1 * this->_imp->a_1(q2) - c_2 * this->_imp->a_2(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::v(const double & q2) const
    {
        return this->_imp->v(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::t_1(const double & q2) const
    {
        return this->_imp->t_1(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::t_2(const double & q2) const
    {
        const double m_B = this->_imp->m_B(), m_V = this->_imp->m_V();
        const double m_B2 = m_B * m_B, m_V2 = m_V * m_V;

        const double c_23a = (m_B2 - m_V2 + q2) / (m_B2 - m_V2);
        const double c_23b = -2.0 * q2 / (m_B2 - m_V2);

        return c_23a * this->_imp->t_23A(q2) + c_23b * this->_imp->t_23B(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::t_3(const double & q2) const
    {
        return -1.0 * this->_imp->t_23A(q2) + 2.0 * this->_imp->t_23B(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::t_23(const double & q2) const
    {
        const double m_B = this->_imp->m_B(), m_V = this->_imp->m_V();
        const double m_B2 = m_B * m_B, m_V2 = m_V * m_V;

        const double c_2 = (m_B + m_V) / (8.0 * m_B * m_V2) * (m_B2 + 3.0 * m_V2 - q2);
        const double c_3 = eos::lambda(m_B2, m_V2, q2) / (8.0 * m_B * m_V2 * (m_B - m_V));

        return c_2 * t_2(q2) - c_3 * t_3(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::a_1_dsr(const double & q2) const
    {
        return this->_imp->a_1_dsr(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::a_2_dsr(const double & q2) const
    {
        return this->_imp->a_2_dsr(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::a_30_dsr(const double & q2) const
    {
        return this->_imp->a_30_dsr(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::v_dsr(const double & q2) const
    {
        return this->_imp->v_dsr(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::t_1_dsr(const double & q2) const
    {
        return this->_imp->t_1_dsr(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::t_23A_dsr(const double & q2) const
    {
        return this->_imp->t_23A_dsr(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::t_23B_dsr(const double & q2) const
    {
        return this->_imp->t_23B_dsr(q2);
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::f_perp(const double &) const
    {
        return 0.0;
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::f_para(const double &) const
    {
        return 0.0;
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::f_long(const double &) const
    {
        return 0.0;
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::f_perp_T(const double &) const
    {
        return 0.0;
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::f_para_T(const double &) const
    {
        return 0.0;
    }

    template <typename Process_>
    double
    AnalyticFormFactorPToVLCSR<Process_>::f_long_T(const double &) const
    {
        return 0.0;
    }

    template <typename Process_>
    Diagnostics
    AnalyticFormFactorPToVLCSR<Process_>::diagnostics() const
    {
        return this->_imp->diagnostics();
    }

    template <typename Process_>
    const std::set<ReferenceName>
    AnalyticFormFactorPToVLCSR<Process_>::references
    {
        "BSZ:2015A"_rn,
    };

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    AnalyticFormFactorPToVLCSR<Process_>::begin_options()
    {
        return Implementation<AnalyticFormFactorPToVLCSR<Process_>>::options.cbegin();
    }

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    AnalyticFormFactorPToVLCSR<Process_>::end_options()
    {
        return Implementation<AnalyticFormFactorPToVLCSR<Process_>>::options.cend();
    }
}

#endif
