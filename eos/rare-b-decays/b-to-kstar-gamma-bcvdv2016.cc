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

#include <eos/rare-b-decays/b-to-kstar-gamma-bcvdv2016.hh>
#include <eos/utils/kinematic.hh>

using namespace std;

namespace eos
{
    using namespace std::placeholders;
    using std::norm;

    BToKstarGammaAmplitudes<tag::BCvDV2016>::BToKstarGammaAmplitudes(const Parameters & p, const Options & o) :
        AmplitudeGenerator(p, o),
        m_b_MSbar(p["mass::b(MSbar)"], *this),
        m_s_MSbar(p["mass::s(2GeV)"], *this),
        formfactor(o, "formfactor", { "BCvDV2016ModelA", "BCvDV2016ModelB", "naive" }, "BCvDV2016ModelA"),
        nonlocal_formfactor(NonlocalFormFactor<nff::PToV>::make(formfactor.value(), p, o))
    {
    }

    WilsonCoefficients<BToS>
    BToKstarGammaAmplitudes<tag::BCvDV2016>::wilson_coefficients() const
    {
        return model->wilson_coefficients_b_to_s(mu(),l.value(), cp_conjugate);
    }

    /* Amplitudes */
    /* calculate contribution to transversity amplitudes
     * -> from el-mg. dipole operator O_7, which factorize naively
     * -> not included 4-quark matrix elements to C_7 -> C_7^eff,
     *    which are captured by the nonlocal formfactor
     * -> use full QCD form factors (FF)
     */
    BToKstarGamma::Amplitudes
    BToKstarGammaAmplitudes<tag::BCvDV2016>::amplitudes() const
    {
        BToKstarGamma::Amplitudes result;

        auto wc = this->wilson_coefficients();

        // classic form factors
        // note that T_2(0) = T_1(0)
        const double
                ff_T1 = form_factors->t_1(0.0),
                ff_T2 = ff_T1;

        // kinematics
        const double
                m_B         = this->m_B(),
                m_B2        = pow(m_B, 2),
                m_V         = this->m_Kstar(),
                m_V2        = pow(m_V, 2),
                lambda      = eos::lambda(m_B2, m_V2, 0.0),
                sqrt_lambda = std::sqrt(lambda);

        // tensorial form factors
        const double
                calF_T_perp = sqrt(2.0) * sqrt_lambda / m_B2 * ff_T1,
                calF_T_para = sqrt(2.0) * (m_B2 - m_V2) / m_B2 * ff_T2;

        const complex<double>
                calH_perp = nonlocal_formfactor->H_perp(0.0),
                calH_para = nonlocal_formfactor->H_para(0.0);

        // Wilson coefficients
        const complex<double>
                mc7_m = (model->m_b_msbar(mu) - model->m_s_msbar(mu)) * (wc.c7() - wc.c7prime()),
                mc7_p = (model->m_b_msbar(mu) + model->m_s_msbar(mu)) * (wc.c7() + wc.c7prime());

        // normalization constant
        const complex<double> calN = g_fermi() * model->ckm_tb() * conj(model->ckm_ts())
                * sqrt(alpha_e * m_B * (m_B * m_B - m_Kstar * m_Kstar) / (128 * pow(M_PI, 4)));

        // vector amplitudes
        result.a_para = +calN * (mc7_m * calF_T_para - 16.0 * pow(M_PI, 2) * m_B * calH_para);
        result.a_perp = -calN * (mc7_p * calF_T_perp - 16.0 * pow(M_PI, 2) * m_B * calH_perp);

        return result;
    }
}
