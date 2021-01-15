/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016, 2017, 2019 Danny van Dyk
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

#include <eos/form-factors/mesonic.hh>
#include <eos/rare-b-decays/b-to-kstar-charmonium.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <complex>

namespace eos
{
    using std::abs;
    using std::arg;
    using std::conj;
    using std::norm;
    using std::real;
    using std::sqrt;

    /*!
     * Implementation for the decay @f$\bar{B} \to \bar{K}^* \psi@f$.
     */
    template <>
    struct Implementation<BToKstarCharmonium>
    {
        UsedParameter g_fermi;

        UsedParameter hbar;

        std::shared_ptr<Model> model;

        SwitchOption q;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_Kstar;

        SwitchOption formfactor;

        NonlocalFormFactorPtr<nc::PToV> nonlocal_formfactor;

        SwitchOption psi;

        UsedParameter m_psi;

        UsedParameter f_psi;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        std::function<complex<double> ()> residue_H_long;
        std::function<complex<double> ()> residue_H_perp;
        std::function<complex<double> ()> residue_H_para;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u),
            model(Model::make(o.get("model", "SM"), p, o)),
            q(o, "q", { "d", "u" }, "d"),
            m_B(p["mass::B_" + q.value()], u),
            tau_B(p["life_time::B_" + q.value()], u),
            m_Kstar(p["mass::K_" + q.value() + "^*"], u),
            formfactor(o, "formfactor", { "GvDV2020", "GRvDV2021" }, "GvDV2020"),
            nonlocal_formfactor(NonlocalFormFactor<nc::PToV>::make("B->K^*::" + formfactor.value(), p, o)),
            psi(o, "psi", { "J/psi", "psi(2S)" }, "J/psi"),
            m_psi(p["mass::" + psi.value()], u),
            f_psi(p["decay-constant::" + psi.value()], u),
            form_factors(FormFactorFactory<PToV>::create("B->K^*::" + o.get("form-factors", "BSZ2015"), p))
        {
            if (! nonlocal_formfactor.get())
                throw InternalError("Cannot construct the nonlocal formfactor");

            if ("J/psi" == psi.value())
            {
                residue_H_long = std::bind(&NonlocalFormFactor<nc::PToV>::H_long_residue_jpsi, nonlocal_formfactor);
                residue_H_perp = std::bind(&NonlocalFormFactor<nc::PToV>::H_perp_residue_jpsi, nonlocal_formfactor);
                residue_H_para = std::bind(&NonlocalFormFactor<nc::PToV>::H_para_residue_jpsi, nonlocal_formfactor);
            }
            else
            {
                residue_H_long = std::bind(&NonlocalFormFactor<nc::PToV>::H_long_residue_psi2s, nonlocal_formfactor);
                residue_H_perp = std::bind(&NonlocalFormFactor<nc::PToV>::H_perp_residue_psi2s, nonlocal_formfactor);
                residue_H_para = std::bind(&NonlocalFormFactor<nc::PToV>::H_para_residue_psi2s, nonlocal_formfactor);
            }

            u.uses(*model);
            u.uses(*nonlocal_formfactor);
            u.uses(*form_factors);
        }

        ~Implementation() = default;

        // The amplitudes in the conventions of [BCvDV2016], eq. (B14)
        struct AmplitudesBCvDV2016
        {
            complex<double> A_perp, A_para, A_long;
        };

        // The amplitudes in the conventions of [T2002], eq. (2.38)
        struct AmplitudesExperimental
        {
            complex<double> A_perp, A_para, A_long;
        };

        AmplitudesBCvDV2016 amplitudes_bcvdv2016() const
        {
            const complex<double> res_H_long = this->residue_H_long();
            const complex<double> res_H_perp = this->residue_H_perp();
            const complex<double> res_H_para = this->residue_H_para();

            const double m_B = this->m_B(), m_B2 = pow(m_B, 2);
            const double m_psi = this->m_psi();

            complex<double> A_perp = m_B2 / (f_psi * m_psi) * res_H_perp;
            complex<double> A_para = m_B2 / (f_psi * m_psi) * res_H_para;
            complex<double> A_long = m_B2 / (f_psi * m_psi) * res_H_long;

            return { A_perp, A_para, A_long };
        }

        // Returns amplitudes in convention of e.g. [T:2002A], eq. (2.38),
        // Amplitudes are CP invariant according to [BRY:2006A].
        AmplitudesExperimental amplitudes_experimental() const
        {
            static complex<double> I(0.0, 1.0);

            const auto amps = this->amplitudes_bcvdv2016();

            return {
                    -I * amps.A_perp,
                    -I * amps.A_para,
                    +I * (m_B() + m_Kstar()) / m_B() * amps.A_long
                };
        }

        double branching_ratio() const
        {
            const auto amps = amplitudes_bcvdv2016();
            const auto lambda = eos::lambda(pow(m_B, 2), pow(m_Kstar, 2), pow(m_psi, 2));
            const auto prefactor = pow(g_fermi * abs(model->ckm_cb() * conj(model->ckm_cs())), 2)
                    * tau_B() / hbar() * sqrt(lambda) / (2.0 * M_PI * m_B);
            const auto r = pow(1.0 + m_Kstar / m_B, 2);

            return prefactor * (norm(amps.A_perp) + norm(amps.A_para) + r * norm(amps.A_long));
        }

        complex<double> ratio_perp() const
        {
            const double m_B2     = pow(m_B, 2);
            const double m_Kstar2 = pow(m_Kstar, 2);
            const double m_psi2   = pow(m_psi, 2);

            const double lambda   = eos::lambda(m_B2, m_Kstar2, m_psi2);

            const double F_perp   = sqrt(2.0 * lambda) / (m_B * (m_B + m_Kstar)) * form_factors->v(m_psi2);

            return residue_H_perp() / F_perp;
        }

        complex<double> ratio_para() const
        {
            const double m_psi2   = pow(m_psi, 2);

            const double F_para   = sqrt(2.0) * (m_B + m_Kstar) / m_B * form_factors->a_1(m_psi2);

            return residue_H_para() / F_para;
        }

        complex<double> ratio_long() const
        {
            const double m_B2     = pow(m_B, 2);
            const double m_Kstar2 = pow(m_Kstar, 2);
            const double m_psi2   = pow(m_psi, 2);

            const double lambda   = eos::lambda(m_B2, m_Kstar2, m_psi2);

            const double F_long   = ((m_B2 - m_Kstar2 - m_psi2) * pow(m_B + m_Kstar, 2) * form_factors->a_1(m_psi2) - lambda * form_factors->a_2(m_psi2))
                                  / (2.0 * m_psi * m_Kstar * pow(m_B + m_Kstar, 2));

            return residue_H_long() / F_long;
        }
    };

    BToKstarCharmonium::BToKstarCharmonium(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BToKstarCharmonium>(new Implementation<BToKstarCharmonium>(p, o, *this))
    {
    }

    BToKstarCharmonium::~BToKstarCharmonium() = default;

    double
    BToKstarCharmonium::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    BToKstarCharmonium::perp_polarization() const
    {
        const auto amps = _imp->amplitudes_experimental();

        return norm(amps.A_perp) / (norm(amps.A_perp) + norm(amps.A_para) + norm(amps.A_long));
    }

    double
    BToKstarCharmonium::para_polarization() const
    {
        const auto amps = _imp->amplitudes_experimental();

        return norm(amps.A_para) / (norm(amps.A_perp) + norm(amps.A_para) + norm(amps.A_long));
    }

    double
    BToKstarCharmonium::long_polarization() const
    {
        const auto amps = _imp->amplitudes_experimental();

        return norm(amps.A_long) / (norm(amps.A_perp) + norm(amps.A_para) + norm(amps.A_long));
    }

    double
    BToKstarCharmonium::delta_perp_long() const
    {
        const auto amps = _imp->amplitudes_experimental();

        const auto result = arg(amps.A_perp / amps.A_long);

        // clamp this result between 0 and 2 pi
        if (result < 0)
        {
            return result + 2.0 * M_PI;
        }
        else
        {
            return result;
        }
    }

    double
    BToKstarCharmonium::delta_para_long() const
    {
        const auto amps = _imp->amplitudes_experimental();

        const auto result = arg(amps.A_para / amps.A_long);

        // clamp this result between -2 pi and 0
        if (result > 0)
        {
            return result - 2.0 * M_PI;
        }
        else
        {
            return result;
        }
    }

    double
    BToKstarCharmonium::S_1c_LHCb() const
    {
        return this->long_polarization();
    }

    double
    BToKstarCharmonium::S_1s_LHCb() const
    {
        return 3.0 / 4.0 * (this->perp_polarization() + this->para_polarization());
    }

    double
    BToKstarCharmonium::S_3_LHCb() const
    {
        return 1.0 / 2.0 * (this->perp_polarization() - this->para_polarization());
    }

    double
    BToKstarCharmonium::S_4_LHCb() const
    {
        return +sqrt(this->long_polarization() * this->para_polarization() / 2.0) * cos(this->delta_para_long());
    }

    double
    BToKstarCharmonium::S_8_LHCb() const
    {
        return +sqrt(this->long_polarization() * this->perp_polarization() / 2.0) * sin(-this->delta_perp_long());
    }

    double
    BToKstarCharmonium::S_9_LHCb() const
    {
        return +sqrt(this->para_polarization() * this->perp_polarization() / 2.0) * sin(this->delta_perp_long() - this->delta_para_long());
    }

    double
    BToKstarCharmonium::re_ratio_perp() const
    {
        return real(_imp->ratio_perp());
    }

    double
    BToKstarCharmonium::re_ratio_para() const
    {
        return real(_imp->ratio_para());
    }

    double
    BToKstarCharmonium::re_ratio_long() const
    {
        return real(_imp->ratio_long());
    }

    double
    BToKstarCharmonium::im_ratio_perp() const
    {
        return imag(_imp->ratio_perp());
    }

    double
    BToKstarCharmonium::im_ratio_para() const
    {
        return imag(_imp->ratio_para());
    }

    double
    BToKstarCharmonium::im_ratio_long() const
    {
        return imag(_imp->ratio_long());
    }

    double
    BToKstarCharmonium::abs_ratio_perp() const
    {
        return abs(_imp->ratio_perp());
    }

    double
    BToKstarCharmonium::abs_ratio_para() const
    {
        return abs(_imp->ratio_para());
    }

    double
    BToKstarCharmonium::abs_ratio_long() const
    {
        return abs(_imp->ratio_long());
    }
}
