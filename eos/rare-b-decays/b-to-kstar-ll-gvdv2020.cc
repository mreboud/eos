/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016, 2017 Danny van Dyk
 * Copyright (c) 2016 Christoph Bobeth
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

#include <eos/rare-b-decays/b-to-kstar-ll-gvdv2020.hh>
#include <eos/utils/kinematic.hh>

using namespace std;

namespace eos
{
    using namespace std::placeholders;
    using std::norm;

    BToKstarDileptonAmplitudes<tag::GvDV2020>::BToKstarDileptonAmplitudes(const Parameters & p, const Options & o) :
        AmplitudeGenerator(p, o),
        hbar(p["hbar"], *this),
        m_b_MSbar(p["mass::b(MSbar)"], *this),
        m_s_MSbar(p["mass::s(2GeV)"], *this),
        mu(p["mu"], *this),
        alpha_e(p["QED::alpha_e(m_b)"], *this),
        g_fermi(p["G_Fermi"], *this),
        tau(p["life_time::B_" + o.get("q", "d")], *this),
        q(o, "q", { "d", "u" }, "d"),
        formfactor(o, "formfactor", { "GvDV2020ModelA", "GvDV2020ModelB", "naive" }, "GvDV2020ModelA"),
        nonlocal_formfactor(NonlocalFormFactor<nc::PToV>::make(formfactor.value(), p, o))
    {
    }

    WilsonCoefficients<BToS>
    BToKstarDileptonAmplitudes<tag::GvDV2020>::wilson_coefficients() const
    {
        return model->wilson_coefficients_b_to_s(mu(), lepton_flavour, cp_conjugate);
    }

    /* Amplitudes */
    /* calculate contribution to transversity amplitudes
     * -> from semi-leptonic operators O_9,10,S,T + primed and
     *    el-mg. dipole operator O_7, which factorize naively
     * -> not included 4-quark matrix elements to C_7,9 -> C_7,9^eff,
     *    which are captured by the nonlocal formfactor
     * -> use full QCD form factors (FF)
     */
    BToKstarDilepton::Amplitudes
    BToKstarDileptonAmplitudes<tag::GvDV2020>::amplitudes(const double & s) const
    {
        BToKstarDilepton::Amplitudes result;

        auto wc = this->wilson_coefficients();

        // classic form factors
        const double
                ff_V  = form_factors->v(s),
                ff_A0 = form_factors->a_0(s),
                ff_A1 = form_factors->a_1(s),
                ff_A2 = form_factors->a_2(s),
                ff_T1 = form_factors->t_1(s),
                ff_T2 = form_factors->t_2(s),
                ff_T3 = form_factors->t_3(s);

        // kinematics
        const double
                sqrt_s      = std::sqrt(s),
                m_B         = this->m_B(),
                m_B2        = pow(m_B, 2),
                m_V         = this->m_Kstar(),
                m_V2        = pow(m_V, 2),
                lambda      = eos::lambda(m_B2, m_V2, s),
                sqrt_lambda = std::sqrt(lambda);

        // vectorial form factors, cf. [GvDV2020], eq. (C5)
        const double
                calF_perp = sqrt(2.0) * sqrt_lambda / (m_B * (m_B + m_V)) * ff_V,
                calF_para = sqrt(2.0) * (m_B + m_V) / m_B * ff_A1,
                calF_long = ((m_B2 - m_V2 - s) * pow(m_B + m_V, 2) * ff_A1 - lambda * ff_A2)
                          / (2.0 * m_V * sqrt_s * pow(m_B + m_V, 2)),
                calF_time = ff_A0;

        // tensorial form factors, cf. [GvDV2020], eq. (C6)
        const double
                calF_T_perp = sqrt(2.0) * sqrt_lambda / m_B2 * ff_T1,
                calF_T_para = sqrt(2.0) * (m_B2 - m_V2) / m_B2 * ff_T2,
                calF_T_long = sqrt_s * (m_B - m_V) / (2.0 * m_B * m_V * pow(m_B2 - m_V2, 2))
                            * ((m_B2 - m_V2) * (m_B2 + 3.0 * m_V2 - s) * ff_T2 - lambda * ff_T3);

        const complex<double>
                calH_perp = nonlocal_formfactor->H_perp(s),
                calH_para = nonlocal_formfactor->H_para(s),
                calH_long = nonlocal_formfactor->H_long(s);

        // Wilson coefficients
        const complex<double>
                c910_m_r = (wc.c9() - wc.c9prime()) + (wc.c10() - wc.c10prime()),
                c910_m_l = (wc.c9() - wc.c9prime()) - (wc.c10() - wc.c10prime()),
                c910_p_r = (wc.c9() + wc.c9prime()) + (wc.c10() + wc.c10prime()),
                c910_p_l = (wc.c9() + wc.c9prime()) - (wc.c10() + wc.c10prime()),
                c7_m = (wc.c7() - wc.c7prime()),
                c7_p = (wc.c7() + wc.c7prime());

        // quark masses
        const double
                m_b_msbar = model->m_b_msbar(mu()),
                m_s_msbar = model->m_s_msbar(mu());

        // normalization constant, cf. [GvDV2020], eq. (B11)
        const double calN = g_fermi() * alpha_e * abs(model->ckm_tb() * conj(model->ckm_ts()))
                * sqrt(s * beta_l(s) * sqrt_lambda / (3.0 * 1024 * pow(M_PI, 5) * m_B));

        // vector amplitudes, cf. [GvDV2020], eq. (B7)-(B9)
        result.a_long_right = -calN * (c910_m_r * calF_long + 2.0 * m_B / s * ((m_b_msbar - m_s_msbar) * c7_m * calF_T_long - 16.0 * pow(M_PI, 2) * m_B * calH_long));
        result.a_long_left  = -calN * (c910_m_l * calF_long + 2.0 * m_B / s * ((m_b_msbar - m_s_msbar) * c7_m * calF_T_long - 16.0 * pow(M_PI, 2) * m_B * calH_long));

        result.a_para_right = -calN * (c910_m_r * calF_para + 2.0 * m_B / s * ((m_b_msbar - m_s_msbar) * c7_m * calF_T_para - 16.0 * pow(M_PI, 2) * m_B * calH_para));
        result.a_para_left  = -calN * (c910_m_l * calF_para + 2.0 * m_B / s * ((m_b_msbar - m_s_msbar) * c7_m * calF_T_para - 16.0 * pow(M_PI, 2) * m_B * calH_para));

        result.a_perp_right = +calN * (c910_p_r * calF_perp + 2.0 * m_B / s * ((m_b_msbar + m_s_msbar) * c7_p * calF_T_perp - 16.0 * pow(M_PI, 2) * m_B * calH_perp));
        result.a_perp_left  = +calN * (c910_p_l * calF_perp + 2.0 * m_B / s * ((m_b_msbar + m_s_msbar) * c7_p * calF_T_perp - 16.0 * pow(M_PI, 2) * m_B * calH_perp));

        // scalar amplitude, cf. [GvDV2020], eq. (B10)
        result.a_time = -calN * 2.0 * (wc.c10() - wc.c10prime()) * calF_time;

#if 0
        // timelike amplitude
        amps.a_time = norm_s * sqrt_lambda / sqrt_s
            * (2.0 * (wc.c10() - wc.c10prime())
               + s / m_l / (m_b_MSbar + m_s_MSbar) * (wc.cP() - wc.cPprime())
              ) * ff_A0;

        // scalar amplitude
        amps.a_scal = -2.0 * norm_s * sqrt_lam * (wc.cS() - wc.cSprime()) / (m_b_MSbar + m_s_MSbar) * ff_A0;

        // tensor amplitudes
        const double
            kin_tensor_1 = norm_s / m_Kstar() * ((m_B2 + 3.0 * m_K2 - s) * ff_T2 - lam_s / m2_diff * ff_T3),
            kin_tensor_2 = 2.0 * norm_s * sqrt_lam / sqrt_s * ff_T1,
            kin_tensor_3 = 2.0 * norm_s * m2_diff / sqrt_s * ff_T2;

        // correct the sign of C_T5 from [BHvD2012v4] because of inconsistent use of gamma5 <-> Levi-Civita
        static const double sign = -1;

        amps.a_para_perp = kin_tensor_1 * wc.cT();
        amps.a_time_long = kin_tensor_1 * sign * wc.cT5();

        amps.a_time_perp = kin_tensor_2 * wc.cT();
        amps.a_long_perp = kin_tensor_2 * sign * wc.cT5();

        amps.a_time_para = kin_tensor_3 * sign * wc.cT5();
        amps.a_long_para = kin_tensor_3 * wc.cT();
#endif

        return result;
    }
}
