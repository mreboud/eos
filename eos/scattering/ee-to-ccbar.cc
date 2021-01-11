/*
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

#include <eos/scattering/ee-to-ccbar.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/kmatrix-impl.hh>

#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <vector>
#include <memory>

namespace eos
{

    template <>
    struct Implementation<EEToCCBar>
    {
//        std::shared_ptr<Model> model;

        UsedParameter hbar;
        UsedParameter alpha_em;
        UsedParameter m_e;
        UsedParameter m_D0;
        UsedParameter m_D;

        //Charmonium masses
        UsedParameter m_psi2S;
        UsedParameter m_psi3770;

        // Channel-Resonance couplings
        UsedParameter g0_psi2S_ee;
        UsedParameter g0_psi3770_ee;
        UsedParameter g0_psi2S_D0Dbar0;
        UsedParameter g0_psi3770_D0Dbar0;
        UsedParameter g0_psi2S_DpDm;
        UsedParameter g0_psi3770_DpDm;


        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
//            model(Model::make(o.get("model", "SM"), p, o)),
            hbar(p["hbar"], u),
            alpha_em(p["QED::alpha_e(0)"], u),
            m_e(p["mass::e"], u),
            m_D0(p["mass::D^0"], u),
            m_D(p["mass::D^+"], u),

            m_psi2S(p["mass::psi(2S)"], u),
            m_psi3770(p["mass::psi(3770)"], u),

            g0_psi2S_ee(p["ee->ccbar::g0(psi(2S),ee)"], u),
            g0_psi3770_ee(p["ee->ccbar::g0(psi(3770),ee)"], u),
            g0_psi2S_D0Dbar0(p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"], u),
            g0_psi3770_D0Dbar0(p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"], u),
            g0_psi2S_DpDm(p["ee->ccbar::g0(psi(2S),D^+D^-)"], u),
            g0_psi3770_DpDm(p["ee->ccbar::g0(psi(3770),D^+D^-)"], u)
        {
        }

        const static unsigned nchannels = 3;
        const static unsigned nresonances = 2;

        inline double sigma_eetomumu(const double & s)
        {

            return 4.0 * M_PI * alpha_em*alpha_em / (3.0 * s);
        }

        inline double sigma_eetof(const double & s, const KMatrix<nchannels, nresonances> & K, const unsigned & channel)
        {

            // Conversion factor between GeV^2 and nb
            const double speedoflight = 299792458.; //Exact value
            const double GeVtonb = 10 * pow( 1.e18 * hbar * speedoflight, 2.0);

            // Channel properties
            const double Nf = K._channels[channel]->_N_orbital;
            const double rhof = real(K._channels[channel]->rho(s));

            // Get T-matrix[ee, channel]
            const complex<double> T1f = K.tmatrix_row(0, s)[channel];

            return GeVtonb * 16.*M_PI/s * Nf * rhof * norm(T1f);
        }

        inline double sigma_eetoccbar(const double & s, const KMatrix<nchannels, nresonances> K)
        {

            double total_xsec = 0.0;

            for (int i = 0; i < nchannels; i++)
            {
                total_xsec += sigma_eetof(s, K, i);
            }

            return total_xsec;
        }

        double sigma_eetochannel(const double & E, const unsigned & channel)
        {

            // Build K Matrix
            auto psi2S_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi2S_res", m_psi2S);
            auto psi3770_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi3770_res", m_psi3770);

            std::vector<Parameter> ee_g0s       {{g0_psi2S_ee, g0_psi3770_ee}};
            std::vector<Parameter> D0Dbar0_g0s  {{g0_psi2S_D0Dbar0, g0_psi3770_D0Dbar0}};
            std::vector<Parameter> DpDm_g0s     {{g0_psi2S_DpDm, g0_psi3770_DpDm}};

            auto ee_chan        = std::make_shared<PPchan<nchannels, nresonances>>("ee_chan", m_e, m_e, 3, ee_g0s);
            auto D0Dbar0_chan   = std::make_shared<PPchan<nchannels, nresonances>>("D0Dbar0_chan", m_D0, m_D0, 3, D0Dbar0_g0s);
            auto DpDm_chan      = std::make_shared<PPchan<nchannels, nresonances>>("DpDm_chan", m_D, m_D, 3, DpDm_g0s);

            KMatrix<nchannels, nresonances> K({ee_chan, D0Dbar0_chan, DpDm_chan}, {psi2S_res, psi3770_res}, "KMatrix");

            return sigma_eetof(E*E, K, channel);
        }

        double Rc(const double & E)
        {

            // Build K Matrix
            auto psi2S_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi2S_res", m_psi2S);
            auto psi3770_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi3770_res", m_psi3770);

            std::vector<Parameter> ee_g0s       {{g0_psi2S_ee, g0_psi3770_ee}};
            std::vector<Parameter> D0Dbar0_g0s  {{g0_psi2S_D0Dbar0, g0_psi3770_D0Dbar0}};
            std::vector<Parameter> DpDm_g0s     {{g0_psi2S_DpDm, g0_psi3770_DpDm}};

            auto ee_chan        = std::make_shared<PPchan<nchannels, nresonances>>("ee_chan", m_e, m_e, 3, ee_g0s);
            auto D0Dbar0_chan   = std::make_shared<PPchan<nchannels, nresonances>>("D0Dbar0_chan", m_D0, m_D0, 3, D0Dbar0_g0s);
            auto DpDm_chan      = std::make_shared<PPchan<nchannels, nresonances>>("DpDm_chan", m_D, m_D, 3, DpDm_g0s);

            KMatrix<nchannels, nresonances> K({ee_chan, D0Dbar0_chan, DpDm_chan}, {psi2S_res, psi3770_res}, "KMatrix");

            return sigma_eetoccbar(E*E, K) / sigma_eetomumu(E*E);
        }


    };


    EEToCCBar::EEToCCBar(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<EEToCCBar>(new Implementation<EEToCCBar>(parameters, options, *this))
    {
    }

    EEToCCBar::~EEToCCBar()
    {
    }

    double
    EEToCCBar::sigma_eetoD0Dbar0(const double & E) const
    {
        return _imp->sigma_eetochannel(E, 1);
    }

    double
    EEToCCBar::Rc(const double & E) const
    {
        return _imp->Rc(E);
    }

}
