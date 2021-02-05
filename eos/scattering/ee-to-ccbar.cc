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

    // K matrix channels
    // 0    ee
    // 1    eff
    // 2    D0Dbar0
    // 3    D+D-

    template <>
    struct Implementation<EEToCCBar>
    {
        UsedParameter hbar;
        UsedParameter alpha_em;
        UsedParameter m_e;
        UsedParameter m_D0;
        UsedParameter m_D;
        UsedParameter m_Dst0;
        UsedParameter m_Dst;
        UsedParameter m_Ds;

        //Charmonium masses
        UsedParameter m_psi2S;
        UsedParameter m_psi3770;

        // Channel-Resonance couplings
        UsedParameter g0_psi2S_ee;
        UsedParameter g0_psi3770_ee;
        UsedParameter g0_psi2S_eff;
        UsedParameter g0_psi3770_eff;
        UsedParameter g0_psi2S_D0Dbar0;
        UsedParameter g0_psi3770_D0Dbar0;
        UsedParameter g0_psi2S_DpDm;
        UsedParameter g0_psi3770_DpDm;
        UsedParameter g0_psi2S_D0Dbarst0;
        UsedParameter g0_psi3770_D0Dbarst0;
        UsedParameter g0_psi2S_DpDstm;
        UsedParameter g0_psi3770_DpDstm;
        UsedParameter g0_psi2S_DspDsm;
        UsedParameter g0_psi3770_DspDsm;

        // Non-cc contribution to the Rc ratio
        UsedParameter Rconstant;


        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            hbar(p["hbar"], u),
            alpha_em(p["QED::alpha_e(0)"], u),
            m_e(p["mass::e"], u),
            m_D0(p["mass::D^0"], u),
            m_D(p["mass::D^+"], u),
            m_Dst0(p["mass::D_u^*"], u),
            m_Dst(p["mass::D_d^*"], u),
            m_Ds(p["mass::D_s"], u),

            m_psi2S(p["mass::psi(2S)"], u),
            m_psi3770(p["mass::psi(3770)"], u),

            g0_psi2S_ee(p["ee->ccbar::g0(psi(2S),ee)"], u),
            g0_psi3770_ee(p["ee->ccbar::g0(psi(3770),ee)"], u),
            g0_psi2S_eff(p["ee->ccbar::g0(psi(2S),eff)"], u),
            g0_psi3770_eff(p["ee->ccbar::g0(psi(3770),eff)"], u),
            g0_psi2S_D0Dbar0(p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"], u),
            g0_psi3770_D0Dbar0(p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"], u),
            g0_psi2S_DpDm(p["ee->ccbar::g0(psi(2S),D^+D^-)"], u),
            g0_psi3770_DpDm(p["ee->ccbar::g0(psi(3770),D^+D^-)"], u),
            g0_psi2S_D0Dbarst0(p["ee->ccbar::g0(psi(2S),D^0Dbar^*0)"], u),
            g0_psi3770_D0Dbarst0(p["ee->ccbar::g0(psi(3770),D^0Dbar^*0)"], u),
            g0_psi2S_DpDstm(p["ee->ccbar::g0(psi(2S),D^+D^*-)"], u),
            g0_psi3770_DpDstm(p["ee->ccbar::g0(psi(3770),D^+D^*-)"], u),
            g0_psi2S_DspDsm(p["ee->ccbar::g0(psi(2S),D_s^+D_s^-)"], u),
            g0_psi3770_DspDsm(p["ee->ccbar::g0(psi(3770),D_s^+D_s^-)"], u),

            Rconstant(p["ee->ccbar::Rconstant"], u)
        {
        }

        const static unsigned nchannels = 7;
        const static unsigned nresonances = 2;

        inline double sigma_eetomumu(const double & s)
        {
            // Conversion factor between GeV^2 and nb
            const double speedoflight = 299792458.; //Exact value
            const double GeVtonb = 10 * pow( 1.e18 * hbar * speedoflight, 2.0);

            return GeVtonb * 4.0 * M_PI * alpha_em*alpha_em / (3.0 * s);
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

        inline double sigma_eetoccbar(const double & s, const KMatrix<nchannels, nresonances> & K)
        {

            double total_xsec = 0.0;

            for (unsigned i = 0; i < nchannels; i++)
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

            std::vector<Parameter> ee_g0s         {{g0_psi2S_ee,      g0_psi3770_ee}};
            std::vector<Parameter> eff_g0s        {{g0_psi2S_eff,     g0_psi3770_eff}};
            std::vector<Parameter> D0Dbar0_g0s    {{g0_psi2S_D0Dbar0, g0_psi3770_D0Dbar0}};
            std::vector<Parameter> DpDm_g0s       {{g0_psi2S_DpDm,    g0_psi3770_DpDm}};
            std::vector<Parameter> D0Dbarst0_g0s  {{g0_psi2S_D0Dbarst0, g0_psi3770_D0Dbarst0}};
            std::vector<Parameter> DpDstm_g0s     {{g0_psi2S_DpDstm,    g0_psi3770_DpDstm}};
            std::vector<Parameter> DspDsm_g0s     {{g0_psi2S_DspDsm,    g0_psi3770_DspDsm}};

            auto ee_chan        = std::make_shared<PPchan<nchannels, nresonances>>("ee_chan", m_e, m_e, 3, ee_g0s);
            // Massless effective channel
            auto eff_chan       = std::make_shared<PPchan<nchannels, nresonances>>("eff_chan", m_e, m_e, 3, eff_g0s);
            auto D0Dbar0_chan   = std::make_shared<PPchan<nchannels, nresonances>>("D0Dbar0_chan", m_D0, m_D0, 3, D0Dbar0_g0s);
            auto DpDm_chan      = std::make_shared<PPchan<nchannels, nresonances>>("DpDm_chan", m_D, m_D, 3, DpDm_g0s);
            auto D0Dbarst0_chan = std::make_shared<VPchan<nchannels, nresonances>>("D0Dbarst0_chan", m_D0, m_Dst0, 3, D0Dbarst0_g0s);
            auto DpDstm_chan    = std::make_shared<VPchan<nchannels, nresonances>>("DpDstm_chan", m_D, m_Dst, 3, DpDstm_g0s);
            auto DspDsm_chan    = std::make_shared<PPchan<nchannels, nresonances>>("DspDsm_chan", m_Ds, m_Ds, 3, DspDsm_g0s);

            KMatrix<nchannels, nresonances> K({ee_chan, eff_chan, D0Dbar0_chan, DpDm_chan, D0Dbarst0_chan, DpDstm_chan, DspDsm_chan},
                                              {psi2S_res, psi3770_res},
                                              "KMatrix");

            return sigma_eetof(E*E, K, channel);
        }

        double Rc(const double & E)
        {

            // Build K Matrix
            auto psi2S_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi2S_res", m_psi2S);
            auto psi3770_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi3770_res", m_psi3770);

            std::vector<Parameter> ee_g0s         {{g0_psi2S_ee,      g0_psi3770_ee}};
            std::vector<Parameter> eff_g0s        {{g0_psi2S_eff,     g0_psi3770_eff}};
            std::vector<Parameter> D0Dbar0_g0s    {{g0_psi2S_D0Dbar0, g0_psi3770_D0Dbar0}};
            std::vector<Parameter> DpDm_g0s       {{g0_psi2S_DpDm,    g0_psi3770_DpDm}};
            std::vector<Parameter> D0Dbarst0_g0s  {{g0_psi2S_D0Dbarst0, g0_psi3770_D0Dbarst0}};
            std::vector<Parameter> DpDstm_g0s     {{g0_psi2S_DpDstm,    g0_psi3770_DpDstm}};
            std::vector<Parameter> DspDsm_g0s     {{g0_psi2S_DspDsm,    g0_psi3770_DspDsm}};

            auto ee_chan        = std::make_shared<PPchan<nchannels, nresonances>>("ee_chan", m_e, m_e, 3, ee_g0s);
            // Massless effective channel
            auto eff_chan       = std::make_shared<PPchan<nchannels, nresonances>>("eff_chan", m_e, m_e, 3, eff_g0s);
            auto D0Dbar0_chan   = std::make_shared<PPchan<nchannels, nresonances>>("D0Dbar0_chan", m_D0, m_D0, 3, D0Dbar0_g0s);
            auto DpDm_chan      = std::make_shared<PPchan<nchannels, nresonances>>("DpDm_chan", m_D, m_D, 3, DpDm_g0s);
            auto D0Dbarst0_chan = std::make_shared<VPchan<nchannels, nresonances>>("D0Dbarst0_chan", m_D0, m_Dst0, 3, D0Dbarst0_g0s);
            auto DpDstm_chan    = std::make_shared<VPchan<nchannels, nresonances>>("DpDstm_chan", m_D, m_Dst, 3, DpDstm_g0s);
            auto DspDsm_chan    = std::make_shared<PPchan<nchannels, nresonances>>("DspDsm_chan", m_Ds, m_Ds, 3, DspDsm_g0s);

            KMatrix<nchannels, nresonances> K({ee_chan, eff_chan, D0Dbar0_chan, DpDm_chan, D0Dbarst0_chan, DpDstm_chan, DspDsm_chan},
                                              {psi2S_res, psi3770_res},
                                              "KMatrix");

            return sigma_eetoccbar(E*E, K) / sigma_eetomumu(E*E) + Rconstant; //Add constant
        }

        double psi2S_partial_width(unsigned channel)
        {

            // Build K Matrix
            auto psi2S_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi2S_res", m_psi2S);
            auto psi3770_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi3770_res", m_psi3770);

            std::vector<Parameter> ee_g0s         {{g0_psi2S_ee,      g0_psi3770_ee}};
            std::vector<Parameter> eff_g0s        {{g0_psi2S_eff,     g0_psi3770_eff}};
            std::vector<Parameter> D0Dbar0_g0s    {{g0_psi2S_D0Dbar0, g0_psi3770_D0Dbar0}};
            std::vector<Parameter> DpDm_g0s       {{g0_psi2S_DpDm,    g0_psi3770_DpDm}};
            std::vector<Parameter> D0Dbarst0_g0s  {{g0_psi2S_D0Dbarst0, g0_psi3770_D0Dbarst0}};
            std::vector<Parameter> DpDstm_g0s     {{g0_psi2S_DpDstm,    g0_psi3770_DpDstm}};
            std::vector<Parameter> DspDsm_g0s     {{g0_psi2S_DspDsm,    g0_psi3770_DspDsm}};

            auto ee_chan        = std::make_shared<PPchan<nchannels, nresonances>>("ee_chan", m_e, m_e, 3, ee_g0s);
            // Massless effective channel
            auto eff_chan       = std::make_shared<PPchan<nchannels, nresonances>>("eff_chan", m_e, m_e, 3, eff_g0s);
            auto D0Dbar0_chan   = std::make_shared<PPchan<nchannels, nresonances>>("D0Dbar0_chan", m_D0, m_D0, 3, D0Dbar0_g0s);
            auto DpDm_chan      = std::make_shared<PPchan<nchannels, nresonances>>("DpDm_chan", m_D, m_D, 3, DpDm_g0s);
            auto D0Dbarst0_chan = std::make_shared<VPchan<nchannels, nresonances>>("D0Dbarst0_chan", m_D0, m_Dst0, 3, D0Dbarst0_g0s);
            auto DpDstm_chan    = std::make_shared<VPchan<nchannels, nresonances>>("DpDstm_chan", m_D, m_Dst, 3, DpDstm_g0s);
            auto DspDsm_chan    = std::make_shared<PPchan<nchannels, nresonances>>("DspDsm_chan", m_Ds, m_Ds, 3, DspDsm_g0s);

            KMatrix<nchannels, nresonances> K({ee_chan, eff_chan, D0Dbar0_chan, DpDm_chan, D0Dbarst0_chan, DpDstm_chan, DspDsm_chan},
                                              {psi2S_res, psi3770_res},
                                              "KMatrix");

            return K.partial_width(0, channel);
        }

        double psi2S_total_width()
        {

            // Build K Matrix
            auto psi2S_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi2S_res", m_psi2S);
            auto psi3770_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi3770_res", m_psi3770);

            std::vector<Parameter> ee_g0s         {{g0_psi2S_ee,      g0_psi3770_ee}};
            std::vector<Parameter> eff_g0s        {{g0_psi2S_eff,     g0_psi3770_eff}};
            std::vector<Parameter> D0Dbar0_g0s    {{g0_psi2S_D0Dbar0, g0_psi3770_D0Dbar0}};
            std::vector<Parameter> DpDm_g0s       {{g0_psi2S_DpDm,    g0_psi3770_DpDm}};
            std::vector<Parameter> D0Dbarst0_g0s  {{g0_psi2S_D0Dbarst0, g0_psi3770_D0Dbarst0}};
            std::vector<Parameter> DpDstm_g0s     {{g0_psi2S_DpDstm,    g0_psi3770_DpDstm}};
            std::vector<Parameter> DspDsm_g0s     {{g0_psi2S_DspDsm,    g0_psi3770_DspDsm}};

            auto ee_chan        = std::make_shared<PPchan<nchannels, nresonances>>("ee_chan", m_e, m_e, 3, ee_g0s);
            // Massless effective channel
            auto eff_chan       = std::make_shared<PPchan<nchannels, nresonances>>("eff_chan", m_e, m_e, 3, eff_g0s);
            auto D0Dbar0_chan   = std::make_shared<PPchan<nchannels, nresonances>>("D0Dbar0_chan", m_D0, m_D0, 3, D0Dbar0_g0s);
            auto DpDm_chan      = std::make_shared<PPchan<nchannels, nresonances>>("DpDm_chan", m_D, m_D, 3, DpDm_g0s);
            auto D0Dbarst0_chan = std::make_shared<VPchan<nchannels, nresonances>>("D0Dbarst0_chan", m_D0, m_Dst0, 3, D0Dbarst0_g0s);
            auto DpDstm_chan    = std::make_shared<VPchan<nchannels, nresonances>>("DpDstm_chan", m_D, m_Dst, 3, DpDstm_g0s);
            auto DspDsm_chan    = std::make_shared<PPchan<nchannels, nresonances>>("DspDsm_chan", m_Ds, m_Ds, 3, DspDsm_g0s);

            KMatrix<nchannels, nresonances> K({ee_chan, eff_chan, D0Dbar0_chan, DpDm_chan, D0Dbarst0_chan, DpDstm_chan, DspDsm_chan},
                                              {psi2S_res, psi3770_res},
                                              "KMatrix");

            return K.width(0);
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
    EEToCCBar::psi2S_ee_width() const
    {
        return _imp->psi2S_partial_width(0);
    }

    double
    EEToCCBar::psi2S_eff_width() const
    {
        return _imp->psi2S_partial_width(1);
    }

    double
    EEToCCBar::psi2S_total_width() const
    {
        return _imp->psi2S_total_width();
    }

    double
    EEToCCBar::sigma_eetoee(const double & E) const
    {
        return _imp->sigma_eetochannel(E, 0);
    }

    double
    EEToCCBar::sigma_eetoeff(const double & E) const
    {
        return _imp->sigma_eetochannel(E, 1);
    }

    double
    EEToCCBar::sigma_eetoD0Dbar0(const double & E) const
    {
        return _imp->sigma_eetochannel(E, 2);
    }

    double
    EEToCCBar::sigma_eetoDpDm(const double & E) const
    {
        return _imp->sigma_eetochannel(E, 3);
    }

    double
    EEToCCBar::sigma_eetoD0Dbarst0(const double & E) const
    {
        return _imp->sigma_eetochannel(E, 4);
    }

    double
    EEToCCBar::sigma_eetoDpDstm(const double & E) const
    {
        return _imp->sigma_eetochannel(E, 5);
    }

    double
    EEToCCBar::sigma_eetoDspDsm(const double & E) const
    {
        return _imp->sigma_eetochannel(E, 6);
    }

    double
    EEToCCBar::Rc(const double & E) const
    {
        return _imp->Rc(E);
    }

}
