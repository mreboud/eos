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
        UsedParameter hbar;
        UsedParameter alpha_em;
        UsedParameter m_e;
        UsedParameter m_D0;
        UsedParameter m_D;
        UsedParameter m_Dst0;
        UsedParameter m_Dst;
        UsedParameter m_Ds;
        UsedParameter m_Dsst;

        //Charmonium masses
        UsedParameter m_psi2S;
        UsedParameter m_psi3770;
        UsedParameter m_psi4040;
        UsedParameter m_psi4160;
        UsedParameter m_psi4415;

        // Channel-Resonance couplings
        UsedParameter g0_psi2S_ee;
        UsedParameter g0_psi3770_ee;
        UsedParameter g0_psi4040_ee;
        UsedParameter g0_psi4160_ee;
        UsedParameter g0_psi4415_ee;

        UsedParameter g0_psi2S_eff;
        UsedParameter g0_psi3770_eff;
        UsedParameter g0_psi4040_eff;
        UsedParameter g0_psi4160_eff;
        UsedParameter g0_psi4415_eff;

        UsedParameter g0_psi2S_D0Dbar0;
        UsedParameter g0_psi3770_D0Dbar0;
        UsedParameter g0_psi4040_D0Dbar0;
        UsedParameter g0_psi4160_D0Dbar0;
        UsedParameter g0_psi4415_D0Dbar0;

        UsedParameter g0_psi2S_D0Dbarst0;
        UsedParameter g0_psi3770_D0Dbarst0;
        UsedParameter g0_psi4040_D0Dbarst0;
        UsedParameter g0_psi4160_D0Dbarst0;
        UsedParameter g0_psi4415_D0Dbarst0;

        UsedParameter g0_psi2S_DspDsm;
        UsedParameter g0_psi3770_DspDsm;
        UsedParameter g0_psi4040_DspDsm;
        UsedParameter g0_psi4160_DspDsm;
        UsedParameter g0_psi4415_DspDsm;

        UsedParameter g0_psi2S_Dst0Dbarst0S0;
        UsedParameter g0_psi3770_Dst0Dbarst0S0;
        UsedParameter g0_psi4040_Dst0Dbarst0S0;
        UsedParameter g0_psi4160_Dst0Dbarst0S0;
        UsedParameter g0_psi4415_Dst0Dbarst0S0;

        UsedParameter g0_psi2S_Dst0Dbarst0S2;
        UsedParameter g0_psi3770_Dst0Dbarst0S2;
        UsedParameter g0_psi4040_Dst0Dbarst0S2;
        UsedParameter g0_psi4160_Dst0Dbarst0S2;
        UsedParameter g0_psi4415_Dst0Dbarst0S2;

        UsedParameter g0_psi2S_DspDsstm;
        UsedParameter g0_psi3770_DspDsstm;
        UsedParameter g0_psi4040_DspDsstm;
        UsedParameter g0_psi4160_DspDsstm;
        UsedParameter g0_psi4415_DspDsstm;

        UsedParameter g0_psi2S_DsstpDsstmS0;
        UsedParameter g0_psi3770_DsstpDsstmS0;
        UsedParameter g0_psi4040_DsstpDsstmS0;
        UsedParameter g0_psi4160_DsstpDsstmS0;
        UsedParameter g0_psi4415_DsstpDsstmS0;

        UsedParameter g0_psi2S_DsstpDsstmS2;
        UsedParameter g0_psi3770_DsstpDsstmS2;
        UsedParameter g0_psi4040_DsstpDsstmS2;
        UsedParameter g0_psi4160_DsstpDsstmS2;
        UsedParameter g0_psi4415_DsstpDsstmS2;

        // Non-cc contribution to the Rc ratio
        UsedParameter c_00_00;
        UsedParameter c_00_01;
        UsedParameter c_00_02;
        UsedParameter c_00_04;
        UsedParameter c_00_08;
        UsedParameter c_00_09;
        UsedParameter c_00_10;
        UsedParameter c_00_15;
        UsedParameter c_00_17;
        UsedParameter c_00_18;

        UsedParameter c_01_01;
        UsedParameter c_01_02;
        UsedParameter c_01_04;
        UsedParameter c_01_08;
        UsedParameter c_01_09;
        UsedParameter c_01_10;
        UsedParameter c_01_15;
        UsedParameter c_01_17;
        UsedParameter c_01_18;

        UsedParameter c_02_02;
        UsedParameter c_02_04;
        UsedParameter c_02_08;
        UsedParameter c_02_09;
        UsedParameter c_02_10;
        UsedParameter c_02_15;
        UsedParameter c_02_17;
        UsedParameter c_02_18;

        UsedParameter c_04_04;
        UsedParameter c_04_08;
        UsedParameter c_04_09;
        UsedParameter c_04_10;
        UsedParameter c_04_15;
        UsedParameter c_04_17;
        UsedParameter c_04_18;

        UsedParameter c_08_08;
        UsedParameter c_08_09;
        UsedParameter c_08_10;
        UsedParameter c_08_15;
        UsedParameter c_08_17;
        UsedParameter c_08_18;

        UsedParameter c_09_09;
        UsedParameter c_09_10;
        UsedParameter c_09_15;
        UsedParameter c_09_17;
        UsedParameter c_09_18;

        UsedParameter c_10_10;
        UsedParameter c_10_15;
        UsedParameter c_10_17;
        UsedParameter c_10_18;

        UsedParameter c_15_15;
        UsedParameter c_15_17;
        UsedParameter c_15_18;

        UsedParameter c_17_17;
        UsedParameter c_17_18;

        UsedParameter c_18_18;

        // Constant terms of the K matrix
        UsedParameter Rconstant;

        const static unsigned nchannels = 20;
        const static unsigned nresonances = 5;

        std::shared_ptr<KMatrix<nchannels, nresonances>> K;

        using IntermediateResult = EEToCCBar::IntermediateResult;
        IntermediateResult _intermediate_result;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            hbar(p["hbar"], u),
            alpha_em(p["QED::alpha_e(0)"], u),
            m_e(p["mass::e"], u),
            m_D0(p["mass::D^0"], u),
            m_D(p["mass::D^+"], u),
            m_Dst0(p["mass::D_u^*"], u),
            m_Dst(p["mass::D_d^*"], u),
            m_Ds(p["mass::D_s"], u),
            m_Dsst(p["mass::D_s^*"], u),

            m_psi2S(p["mass::psi(2S)"], u),
            m_psi3770(p["mass::psi(3770)"], u),
            m_psi4040(p["mass::psi(4040)"], u),
            m_psi4160(p["mass::psi(4160)"], u),
            m_psi4415(p["mass::psi(4415)"], u),

            g0_psi2S_ee(p["ee->ccbar::g0(psi(2S),ee)"], u),
            g0_psi3770_ee(p["ee->ccbar::g0(psi(3770),ee)"], u),
            g0_psi4040_ee(p["ee->ccbar::g0(psi(4040),ee)"], u),
            g0_psi4160_ee(p["ee->ccbar::g0(psi(4160),ee)"], u),
            g0_psi4415_ee(p["ee->ccbar::g0(psi(4415),ee)"], u),

            g0_psi2S_eff(p["ee->ccbar::g0(psi(2S),eff)"], u),
            g0_psi3770_eff(p["ee->ccbar::g0(psi(3770),eff)"], u),
            g0_psi4040_eff(p["ee->ccbar::g0(psi(4040),eff)"], u),
            g0_psi4160_eff(p["ee->ccbar::g0(psi(4160),eff)"], u),
            g0_psi4415_eff(p["ee->ccbar::g0(psi(4415),eff)"], u),

            g0_psi2S_D0Dbar0(p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"], u),
            g0_psi3770_D0Dbar0(p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"], u),
            g0_psi4040_D0Dbar0(p["ee->ccbar::g0(psi(4040),D^0Dbar^0)"], u),
            g0_psi4160_D0Dbar0(p["ee->ccbar::g0(psi(4160),D^0Dbar^0)"], u),
            g0_psi4415_D0Dbar0(p["ee->ccbar::g0(psi(4415),D^0Dbar^0)"], u),

            g0_psi2S_D0Dbarst0(p["ee->ccbar::g0(psi(2S),D^0Dbar^*0)"], u),
            g0_psi3770_D0Dbarst0(p["ee->ccbar::g0(psi(3770),D^0Dbar^*0)"], u),
            g0_psi4040_D0Dbarst0(p["ee->ccbar::g0(psi(4040),D^0Dbar^*0)"], u),
            g0_psi4160_D0Dbarst0(p["ee->ccbar::g0(psi(4160),D^0Dbar^*0)"], u),
            g0_psi4415_D0Dbarst0(p["ee->ccbar::g0(psi(4415),D^0Dbar^*0)"], u),

            g0_psi2S_DspDsm(p["ee->ccbar::g0(psi(2S),D_s^+D_s^-)"], u),
            g0_psi3770_DspDsm(p["ee->ccbar::g0(psi(3770),D_s^+D_s^-)"], u),
            g0_psi4040_DspDsm(p["ee->ccbar::g0(psi(4040),D_s^+D_s^-)"], u),
            g0_psi4160_DspDsm(p["ee->ccbar::g0(psi(4160),D_s^+D_s^-)"], u),
            g0_psi4415_DspDsm(p["ee->ccbar::g0(psi(4415),D_s^+D_s^-)"], u),

            g0_psi2S_Dst0Dbarst0S0(p["ee->ccbar::g0(psi(2S),D^*0Dbar^*0S0)"], u),
            g0_psi3770_Dst0Dbarst0S0(p["ee->ccbar::g0(psi(3770),D^*0Dbar^*0S0)"], u),
            g0_psi4040_Dst0Dbarst0S0(p["ee->ccbar::g0(psi(4040),D^*0Dbar^*0S0)"], u),
            g0_psi4160_Dst0Dbarst0S0(p["ee->ccbar::g0(psi(4160),D^*0Dbar^*0S0)"], u),
            g0_psi4415_Dst0Dbarst0S0(p["ee->ccbar::g0(psi(4415),D^*0Dbar^*0S0)"], u),

            g0_psi2S_Dst0Dbarst0S2(p["ee->ccbar::g0(psi(2S),D^*0Dbar^*0S2)"], u),
            g0_psi3770_Dst0Dbarst0S2(p["ee->ccbar::g0(psi(3770),D^*0Dbar^*0S2)"], u),
            g0_psi4040_Dst0Dbarst0S2(p["ee->ccbar::g0(psi(4040),D^*0Dbar^*0S2)"], u),
            g0_psi4160_Dst0Dbarst0S2(p["ee->ccbar::g0(psi(4160),D^*0Dbar^*0S2)"], u),
            g0_psi4415_Dst0Dbarst0S2(p["ee->ccbar::g0(psi(4415),D^*0Dbar^*0S2)"], u),

            g0_psi2S_DspDsstm(p["ee->ccbar::g0(psi(2S),D_s^+D_s^*-)"], u),
            g0_psi3770_DspDsstm(p["ee->ccbar::g0(psi(3770),D_s^+D_s^*-)"], u),
            g0_psi4040_DspDsstm(p["ee->ccbar::g0(psi(4040),D_s^+D_s^*-)"], u),
            g0_psi4160_DspDsstm(p["ee->ccbar::g0(psi(4160),D_s^+D_s^*-)"], u),
            g0_psi4415_DspDsstm(p["ee->ccbar::g0(psi(4415),D_s^+D_s^*-)"], u),

            g0_psi2S_DsstpDsstmS0(p["ee->ccbar::g0(psi(2S),D_s^*+D_s^*-S0)"], u),
            g0_psi3770_DsstpDsstmS0(p["ee->ccbar::g0(psi(3770),D_s^*+D_s^*-S0)"], u),
            g0_psi4040_DsstpDsstmS0(p["ee->ccbar::g0(psi(4040),D_s^*+D_s^*-S0)"], u),
            g0_psi4160_DsstpDsstmS0(p["ee->ccbar::g0(psi(4160),D_s^*+D_s^*-S0)"], u),
            g0_psi4415_DsstpDsstmS0(p["ee->ccbar::g0(psi(4415),D_s^*+D_s^*-S0)"], u),

            g0_psi2S_DsstpDsstmS2(p["ee->ccbar::g0(psi(2S),D_s^*+D_s^*-S2)"], u),
            g0_psi3770_DsstpDsstmS2(p["ee->ccbar::g0(psi(3770),D_s^*+D_s^*-S2)"], u),
            g0_psi4040_DsstpDsstmS2(p["ee->ccbar::g0(psi(4040),D_s^*+D_s^*-S2)"], u),
            g0_psi4160_DsstpDsstmS2(p["ee->ccbar::g0(psi(4160),D_s^*+D_s^*-S2)"], u),
            g0_psi4415_DsstpDsstmS2(p["ee->ccbar::g0(psi(4415),D_s^*+D_s^*-S2)"], u),

            c_00_00(p["ee->ccbar::c(ee,ee)"], u),
            c_00_01(p["ee->ccbar::c(ee,eff)"], u),
            c_00_02(p["ee->ccbar::c(ee,D^0Dbar^0)"], u),
            c_00_04(p["ee->ccbar::c(ee,D^0Dbar^*0)"], u),
            c_00_08(p["ee->ccbar::c(ee,D_s^+D_s^-)"], u),
            c_00_09(p["ee->ccbar::c(ee,D^*0Dbar^*0S0)"], u),
            c_00_10(p["ee->ccbar::c(ee,D^*0Dbar^*0S2)"], u),
            c_00_15(p["ee->ccbar::c(ee,D_s^+D_s^*-)"], u),
            c_00_17(p["ee->ccbar::c(ee,D_s^*+D_s^*-S0)"], u),
            c_00_18(p["ee->ccbar::c(ee,D_s^*+D_s^*-S2)"], u),

            c_01_01(p["ee->ccbar::c(eff,eff)"], u),
            c_01_02(p["ee->ccbar::c(eff,D^0Dbar^0)"], u),
            c_01_04(p["ee->ccbar::c(eff,D^0Dbar^*0)"], u),
            c_01_08(p["ee->ccbar::c(eff,D_s^+D_s^-)"], u),
            c_01_09(p["ee->ccbar::c(eff,D^*0Dbar^*0S0)"], u),
            c_01_10(p["ee->ccbar::c(eff,D^*0Dbar^*0S2)"], u),
            c_01_15(p["ee->ccbar::c(eff,D_s^+D_s^*-)"], u),
            c_01_17(p["ee->ccbar::c(eff,D_s^*+D_s^*-S0)"], u),
            c_01_18(p["ee->ccbar::c(eff,D_s^*+D_s^*-S2)"], u),

            c_02_02(p["ee->ccbar::c(D^0Dbar^0,D^0Dbar^0)"], u),
            c_02_04(p["ee->ccbar::c(D^0Dbar^0,D^0Dbar^*0)"], u),
            c_02_08(p["ee->ccbar::c(D^0Dbar^0,D_s^+D_s^-)"], u),
            c_02_09(p["ee->ccbar::c(D^0Dbar^0,D^*0Dbar^*0S0)"], u),
            c_02_10(p["ee->ccbar::c(D^0Dbar^0,D^*0Dbar^*0S2)"], u),
            c_02_15(p["ee->ccbar::c(D^0Dbar^0,D_s^+D_s^*-)"], u),
            c_02_17(p["ee->ccbar::c(D^0Dbar^0,D_s^*+D_s^*-S0)"], u),
            c_02_18(p["ee->ccbar::c(D^0Dbar^0,D_s^*+D_s^*-S2)"], u),

            c_04_04(p["ee->ccbar::c(D^0Dbar^*0,D^0Dbar^*0)"], u),
            c_04_08(p["ee->ccbar::c(D^0Dbar^*0,D_s^+D_s^-)"], u),
            c_04_09(p["ee->ccbar::c(D^0Dbar^*0,D^*0Dbar^*0S0)"], u),
            c_04_10(p["ee->ccbar::c(D^0Dbar^*0,D^*0Dbar^*0S2)"], u),
            c_04_15(p["ee->ccbar::c(D^0Dbar^*0,D_s^+D_s^*-)"], u),
            c_04_17(p["ee->ccbar::c(D^0Dbar^*0,D_s^*+D_s^*-S0)"], u),
            c_04_18(p["ee->ccbar::c(D^0Dbar^*0,D_s^*+D_s^*-S2)"], u),

            c_08_08(p["ee->ccbar::c(D_s^+D_s^-,D_s^+D_s^-)"], u),
            c_08_09(p["ee->ccbar::c(D_s^+D_s^-,D^*0Dbar^*0S0)"], u),
            c_08_10(p["ee->ccbar::c(D_s^+D_s^-,D^*0Dbar^*0S2)"], u),
            c_08_15(p["ee->ccbar::c(D_s^+D_s^-,D_s^+D_s^*-)"], u),
            c_08_17(p["ee->ccbar::c(D_s^+D_s^-,D_s^*+D_s^*-S0)"], u),
            c_08_18(p["ee->ccbar::c(D_s^+D_s^-,D_s^*+D_s^*-S2)"], u),

            c_09_09(p["ee->ccbar::c(D^*0Dbar^*0S0,D^*0Dbar^*0S0)"], u),
            c_09_10(p["ee->ccbar::c(D^*0Dbar^*0S0,D^*0Dbar^*0S2)"], u),
            c_09_15(p["ee->ccbar::c(D^*0Dbar^*0S0,D_s^+D_s^*-)"], u),
            c_09_17(p["ee->ccbar::c(D^*0Dbar^*0S0,D_s^*+D_s^*-S0)"], u),
            c_09_18(p["ee->ccbar::c(D^*0Dbar^*0S0,D_s^*+D_s^*-S2)"], u),

            c_10_10(p["ee->ccbar::c(D^*0Dbar^*0S2,D^*0Dbar^*0S2)"], u),
            c_10_15(p["ee->ccbar::c(D^*0Dbar^*0S2,D_s^+D_s^*-)"], u),
            c_10_17(p["ee->ccbar::c(D^*0Dbar^*0S2,D_s^*+D_s^*-S0)"], u),
            c_10_18(p["ee->ccbar::c(D^*0Dbar^*0S2,D_s^*+D_s^*-S2)"], u),

            c_15_15(p["ee->ccbar::c(D_s^+D_s^*-,D_s^+D_s^*-)"], u),
            c_15_17(p["ee->ccbar::c(D_s^+D_s^*-,D_s^*+D_s^*-S0)"], u),
            c_15_18(p["ee->ccbar::c(D_s^+D_s^*-,D_s^*+D_s^*-S2)"], u),

            c_17_17(p["ee->ccbar::c(D_s^*+D_s^*-S0,D_s^*+D_s^*-S0)"], u),
            c_17_18(p["ee->ccbar::c(D_s^*+D_s^*-S0,D_s^*+D_s^*-S2)"], u),

            c_18_18(p["ee->ccbar::c(D_s^*+D_s^*-S2,D_s^*+D_s^*-S2)"], u),

            Rconstant(p["ee->ccbar::Rconstant"], u)
        {
            // Build K Matrix IN THE ISOSPIN LIMIT
            auto psi2S_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi2S_res", m_psi2S);
            auto psi3770_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi3770_res", m_psi3770);
            auto psi4040_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi4040_res", m_psi4040);
            auto psi4160_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi4160_res", m_psi4160);
            auto psi4415_res = std::make_shared<charmonium_resonance<nchannels, nresonances>>("psi4415_res", m_psi4415);

            std::vector<std::vector<Parameter>> bkgcst {
                {c_00_00, c_00_01, c_00_02, c_00_02, c_00_04, c_00_04, c_00_04, c_00_04, c_00_08, c_00_09, c_00_10, c_00_10, c_00_10, c_00_10, c_00_10, c_00_15, c_00_15, c_00_17, c_00_18, c_00_18},
                {c_00_01, c_01_01, c_01_02, c_01_02, c_01_04, c_01_04, c_01_04, c_01_04, c_01_08, c_01_09, c_01_10, c_01_10, c_01_10, c_01_10, c_01_10, c_01_15, c_01_15, c_01_17, c_01_18, c_01_18},
                {c_00_02, c_01_02, c_02_02, c_02_02, c_02_04, c_02_04, c_02_04, c_02_04, c_02_08, c_02_09, c_02_10, c_02_10, c_02_10, c_02_10, c_02_10, c_02_15, c_02_15, c_02_17, c_02_18, c_02_18},
                {c_00_02, c_01_02, c_02_02, c_02_02, c_02_04, c_02_04, c_02_04, c_02_04, c_02_08, c_02_09, c_02_10, c_02_10, c_02_10, c_02_10, c_02_10, c_02_15, c_02_15, c_02_17, c_02_18, c_02_18},
                {c_00_04, c_01_04, c_02_04, c_02_04, c_04_04, c_04_04, c_04_04, c_04_04, c_04_08, c_04_09, c_04_10, c_04_10, c_04_10, c_04_10, c_04_10, c_04_15, c_04_15, c_04_17, c_04_18, c_04_18},
                {c_00_04, c_01_04, c_02_04, c_02_04, c_04_04, c_04_04, c_04_04, c_04_04, c_04_08, c_04_09, c_04_10, c_04_10, c_04_10, c_04_10, c_04_10, c_04_15, c_04_15, c_04_17, c_04_18, c_04_18},
                {c_00_04, c_01_04, c_02_04, c_02_04, c_04_04, c_04_04, c_04_04, c_04_04, c_04_08, c_04_09, c_04_10, c_04_10, c_04_10, c_04_10, c_04_10, c_04_15, c_04_15, c_04_17, c_04_18, c_04_18},
                {c_00_04, c_01_04, c_02_04, c_02_04, c_04_04, c_04_04, c_04_04, c_04_04, c_04_08, c_04_09, c_04_10, c_04_10, c_04_10, c_04_10, c_04_10, c_04_15, c_04_15, c_04_17, c_04_18, c_04_18},
                {c_00_08, c_01_08, c_02_08, c_02_08, c_04_08, c_04_08, c_04_08, c_04_08, c_08_08, c_08_09, c_08_10, c_08_10, c_08_10, c_08_10, c_08_10, c_08_15, c_08_15, c_08_17, c_08_18, c_08_18},
                {c_00_09, c_01_09, c_02_09, c_02_09, c_04_09, c_04_09, c_04_09, c_04_09, c_08_09, c_09_09, c_09_10, c_09_10, c_09_10, c_09_10, c_09_10, c_09_15, c_09_15, c_09_17, c_09_18, c_09_18},
                {c_00_10, c_01_10, c_02_10, c_02_10, c_04_10, c_04_10, c_04_10, c_04_10, c_08_10, c_09_10, c_10_10, c_10_10, c_10_10, c_10_10, c_10_10, c_10_15, c_10_15, c_10_17, c_10_18, c_10_18},
                {c_00_10, c_01_10, c_02_10, c_02_10, c_04_10, c_04_10, c_04_10, c_04_10, c_08_10, c_09_10, c_10_10, c_10_10, c_10_10, c_10_10, c_10_10, c_10_15, c_10_15, c_10_17, c_10_18, c_10_18},
                {c_00_10, c_01_10, c_02_10, c_02_10, c_04_10, c_04_10, c_04_10, c_04_10, c_08_10, c_09_10, c_10_10, c_10_10, c_10_10, c_10_10, c_10_10, c_10_15, c_10_15, c_10_17, c_10_18, c_10_18},
                {c_00_10, c_01_10, c_02_10, c_02_10, c_04_10, c_04_10, c_04_10, c_04_10, c_08_10, c_09_10, c_10_10, c_10_10, c_10_10, c_10_10, c_10_10, c_10_15, c_10_15, c_10_17, c_10_18, c_10_18},
                {c_00_10, c_01_10, c_02_10, c_02_10, c_04_10, c_04_10, c_04_10, c_04_10, c_08_10, c_09_10, c_10_10, c_10_10, c_10_10, c_10_10, c_10_10, c_10_15, c_10_15, c_10_17, c_10_18, c_10_18},
                {c_00_15, c_01_15, c_02_15, c_02_15, c_04_15, c_04_15, c_04_15, c_04_15, c_08_15, c_09_15, c_10_15, c_10_15, c_10_15, c_10_15, c_10_15, c_15_15, c_15_15, c_15_17, c_15_18, c_15_18},
                {c_00_15, c_01_15, c_02_15, c_02_15, c_04_15, c_04_15, c_04_15, c_04_15, c_08_15, c_09_15, c_10_15, c_10_15, c_10_15, c_10_15, c_10_15, c_15_15, c_15_15, c_15_17, c_15_18, c_15_18},
                {c_00_17, c_01_17, c_02_17, c_02_17, c_04_17, c_04_17, c_04_17, c_04_17, c_08_17, c_09_17, c_10_17, c_10_17, c_10_17, c_10_17, c_10_17, c_15_17, c_15_17, c_17_17, c_17_18, c_17_18},
                {c_00_18, c_01_18, c_02_18, c_02_18, c_04_18, c_04_18, c_04_18, c_04_18, c_08_18, c_09_18, c_10_18, c_10_18, c_10_18, c_10_18, c_10_18, c_15_18, c_15_18, c_17_18, c_18_18, c_18_18},
                {c_00_18, c_01_18, c_02_18, c_02_18, c_04_18, c_04_18, c_04_18, c_04_18, c_08_18, c_09_18, c_10_18, c_10_18, c_10_18, c_10_18, c_10_18, c_15_18, c_15_18, c_17_18, c_18_18, c_18_18}
            };

            std::vector<Parameter> ee_g0s            {{g0_psi2S_ee,            g0_psi3770_ee,            g0_psi4040_ee,            g0_psi4160_ee,            g0_psi4415_ee            }};
            std::vector<Parameter> eff_g0s           {{g0_psi2S_eff,           g0_psi3770_eff,           g0_psi4040_eff,           g0_psi4160_eff,           g0_psi4415_eff           }};
            std::vector<Parameter> D0Dbar0_g0s       {{g0_psi2S_D0Dbar0,       g0_psi3770_D0Dbar0,       g0_psi4040_D0Dbar0,       g0_psi4160_D0Dbar0,       g0_psi4415_D0Dbar0       }};
            std::vector<Parameter> D0Dbarst0_g0s     {{g0_psi2S_D0Dbarst0,     g0_psi3770_D0Dbarst0,     g0_psi4040_D0Dbarst0,     g0_psi4160_D0Dbarst0,     g0_psi4415_D0Dbarst0     }};
            std::vector<Parameter> DspDsm_g0s        {{g0_psi2S_DspDsm,        g0_psi3770_DspDsm,        g0_psi4040_DspDsm,        g0_psi4160_DspDsm,        g0_psi4415_DspDsm        }};
            std::vector<Parameter> Dst0Dbarst0S0_g0s {{g0_psi2S_Dst0Dbarst0S0, g0_psi3770_Dst0Dbarst0S0, g0_psi4040_Dst0Dbarst0S0, g0_psi4160_Dst0Dbarst0S0, g0_psi4415_Dst0Dbarst0S0 }};
            std::vector<Parameter> Dst0Dbarst0S2_g0s {{g0_psi2S_Dst0Dbarst0S2, g0_psi3770_Dst0Dbarst0S2, g0_psi4040_Dst0Dbarst0S2, g0_psi4160_Dst0Dbarst0S2, g0_psi4415_Dst0Dbarst0S2 }};
            std::vector<Parameter> DspDsstm_g0s      {{g0_psi2S_DspDsstm,      g0_psi3770_DspDsstm,      g0_psi4040_DspDsstm,      g0_psi4160_DspDsstm,      g0_psi4415_DspDsstm      }};
            std::vector<Parameter> DsstpDsstmS0_g0s  {{g0_psi2S_DsstpDsstmS0,  g0_psi3770_DsstpDsstmS0,  g0_psi4040_DsstpDsstmS0,  g0_psi4160_DsstpDsstmS0,  g0_psi4415_DsstpDsstmS0  }};
            std::vector<Parameter> DsstpDsstmS2_g0s  {{g0_psi2S_DsstpDsstmS2,  g0_psi3770_DsstpDsstmS2,  g0_psi4040_DsstpDsstmS2,  g0_psi4160_DsstpDsstmS2,  g0_psi4415_DsstpDsstmS2  }};

            auto ee_chan            = std::make_shared<PPPwavechan<nchannels, nresonances>>("ee_chan",            m_e,    m_e,     3, ee_g0s);
            // Massless effective channel
            auto eff_chan           = std::make_shared<PPPwavechan<nchannels, nresonances>>("eff_chan",           m_e,    m_e,     3, eff_g0s);
            auto D0Dbar0_chan       = std::make_shared<PPPwavechan<nchannels, nresonances>>("D0Dbar0_chan",       m_D0,   m_D0,    3, D0Dbar0_g0s);
            auto DpDm_chan          = std::make_shared<PPPwavechan<nchannels, nresonances>>("DpDm_chan",          m_D,    m_D,     3, D0Dbar0_g0s);
            auto D0Dbarst0_chan     = std::make_shared<VPPwavechan<nchannels, nresonances>>("D0Dbarst0_chan",     m_D0,   m_Dst0,  3, D0Dbarst0_g0s);
            auto Dst0Dbar0_chan     = std::make_shared<VPPwavechan<nchannels, nresonances>>("Dst0Dbar0_chan",     m_Dst0, m_D0,    3, D0Dbarst0_g0s);
            auto DpDstm_chan        = std::make_shared<VPPwavechan<nchannels, nresonances>>("DpDstm_chan",        m_D,    m_Dst,   3, D0Dbarst0_g0s);
            auto DstpDm_chan        = std::make_shared<VPPwavechan<nchannels, nresonances>>("DstpDm_chan",        m_D,    m_Dst,   3, D0Dbarst0_g0s);
            auto DspDsm_chan        = std::make_shared<PPPwavechan<nchannels, nresonances>>("DspDsm_chan",        m_Ds,   m_Ds,    3, DspDsm_g0s);
            auto Dst0Dbarst0P0_chan = std::make_shared<VVPwavechan<nchannels, nresonances>>("Dst0Dbarst0P0_chan", m_Dst0, m_Dst0,  3, Dst0Dbarst0S0_g0s);
            auto Dst0Dbarst0P2_chan = std::make_shared<VVPwavechan<nchannels, nresonances>>("Dst0Dbarst0P2_chan", m_Dst0, m_Dst0,  3, Dst0Dbarst0S2_g0s);
            auto Dst0Dbarst0F2_chan = std::make_shared<VVFwavechan<nchannels, nresonances>>("Dst0Dbarst0F2_chan", m_Dst0, m_Dst0,  7, Dst0Dbarst0S2_g0s);
            auto DstpDstmP0_chan    = std::make_shared<VVPwavechan<nchannels, nresonances>>("DstpDstmP0_chan",    m_Dst,  m_Dst,   3, Dst0Dbarst0S0_g0s);
            auto DstpDstmP2_chan    = std::make_shared<VVPwavechan<nchannels, nresonances>>("DstpDstmP2_chan",    m_Dst,  m_Dst,   3, Dst0Dbarst0S2_g0s);
            auto DstpDstmF2_chan    = std::make_shared<VVFwavechan<nchannels, nresonances>>("DstpDstmF2_chan",    m_Dst,  m_Dst,   7, Dst0Dbarst0S2_g0s);
            auto DspDsstm_chan      = std::make_shared<VPPwavechan<nchannels, nresonances>>("DspDsstm_chan",      m_Ds,   m_Dsst,  3, DspDsstm_g0s);
            auto DsstpDsm_chan      = std::make_shared<VPPwavechan<nchannels, nresonances>>("DsstpDsm_chan",      m_Dsst, m_Ds,    3, DspDsstm_g0s);
            auto DsstpDsstmP0_chan  = std::make_shared<VVPwavechan<nchannels, nresonances>>("DsstpDsstmP0_chan",  m_Dsst, m_Dsst,  3, DsstpDsstmS0_g0s);
            auto DsstpDsstmP2_chan  = std::make_shared<VVPwavechan<nchannels, nresonances>>("DsstpDsstmP2_chan",  m_Dsst, m_Dsst,  3, DsstpDsstmS2_g0s);
            auto DsstpDsstmF2_chan  = std::make_shared<VVFwavechan<nchannels, nresonances>>("DsstpDsstmF2_chan",  m_Dsst, m_Dsst,  7, DsstpDsstmS2_g0s);

            K = std::shared_ptr<KMatrix<nchannels, nresonances>> (
                new KMatrix<nchannels, nresonances>(
                    { ee_chan,
                      eff_chan,
                      D0Dbar0_chan,
                      DpDm_chan,
                      D0Dbarst0_chan,
                      Dst0Dbar0_chan,
                      DpDstm_chan,
                      DstpDm_chan,
                      DspDsm_chan,
                      Dst0Dbarst0P0_chan,
                      Dst0Dbarst0P2_chan,
                      Dst0Dbarst0F2_chan,
                      DstpDstmP0_chan,
                      DstpDstmP2_chan,
                      DstpDstmF2_chan,
                      DspDsstm_chan,
                      DsstpDsm_chan,
                      DsstpDsstmP0_chan,
                      DsstpDsstmP2_chan,
                      DsstpDsstmF2_chan
                    },
                    { psi2S_res,
                      psi3770_res,
                      psi4040_res,
                      psi4160_res,
                      psi4415_res
                    },
                    bkgcst,
                    "KMatrix")
                );
        }

        const IntermediateResult * prepare(const double & E)
        {
            _intermediate_result.tmatrix_row_0 = K->tmatrix_row(0, E*E);

            _intermediate_result.E = E;
            _intermediate_result.s = E*E;

            return &_intermediate_result;
        }


        inline double sigma_eetomumu(const double & E)
        {
            // Conversion factor between GeV^2 and nb
            const double speedoflight = 299792458.; //Exact value
            const double GeVtonb = 10 * pow( 1.e18 * hbar * speedoflight, 2.0);

            return GeVtonb * 4.0 * M_PI * alpha_em*alpha_em / (3.0 * E*E);
        }

        double sigma_eetochannel(const IntermediateResult * intermediate_result, const unsigned & channel)
        {
            // Conversion factor between GeV^2 and nb
            const double speedoflight = 299792458.; //Exact value
            const double GeVtonb = 10 * pow( 1.e18 * hbar * speedoflight, 2.0);

            // Channel properties
            const double Nf = K->_channels[channel]->_N_orbital;
            const double rhof = real(K->_channels[channel]->rho(intermediate_result->s));

            // Get T-matrix[ee, channel]
            const complex<double> T1f = intermediate_result->tmatrix_row_0[channel];

            return GeVtonb * 16. * M_PI / intermediate_result->s * Nf * rhof * norm(T1f);
        }

        double Rc(const IntermediateResult * intermediate_result)
        {
            double total_xsec = 0.0;

            for (unsigned i = 0; i < nchannels; i++)
            {
                total_xsec += sigma_eetochannel(intermediate_result, i);
            }

            return total_xsec / sigma_eetomumu(intermediate_result->E) + Rconstant; //Add constant term
        }

        double psi2S_partial_width(unsigned channel)
        {
            return K->partial_width(0, channel);
        }

        double psi2S_total_width()
        {
            return K->width(0);
        }


    };


    EEToCCBar::EEToCCBar(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<EEToCCBar>(new Implementation<EEToCCBar>(parameters, options, *this))
    {
    }

    EEToCCBar::~EEToCCBar()
    {
    }

    const EEToCCBar::IntermediateResult *
    EEToCCBar::prepare(const double & E) const
    {
        return _imp->prepare(E);
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
    EEToCCBar::sigma_eetoee(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->sigma_eetochannel(ir, 0);
    }

    double
    EEToCCBar::sigma_eetoeff(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->sigma_eetochannel(ir, 1);
    }

    double
    EEToCCBar::sigma_eetoD0Dbar0(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->sigma_eetochannel(ir, 2);
    }

    double
    EEToCCBar::sigma_eetoDpDm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->sigma_eetochannel(ir, 3);
    }

    double
    EEToCCBar::sigma_eetoD0Dbarst0(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->sigma_eetochannel(ir, 4) + _imp->sigma_eetochannel(ir, 5);
    }

    double
    EEToCCBar::sigma_eetoDpDstm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->sigma_eetochannel(ir, 6) + _imp->sigma_eetochannel(ir, 7);
    }

    double
    EEToCCBar::sigma_eetoDspDsm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->sigma_eetochannel(ir, 8);
    }

    double
    EEToCCBar::sigma_eetoDst0Dbarst0(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->sigma_eetochannel(ir, 9) + _imp->sigma_eetochannel(ir, 10);
    }

    double
    EEToCCBar::sigma_eetoDstpDstm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->sigma_eetochannel(ir, 11) + _imp->sigma_eetochannel(ir, 12);
    }

    double
    EEToCCBar::sigma_eetoDspDsstm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->sigma_eetochannel(ir, 13) + _imp->sigma_eetochannel(ir, 14);
    }

    double
    EEToCCBar::sigma_eetoDsstpDsstm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->sigma_eetochannel(ir, 15) + _imp->sigma_eetochannel(ir, 16);
    }

    double
    EEToCCBar::Rc(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->Rc(ir);
    }


    const std::set<ReferenceName>
    EEToCCBar::references
    {
    };

}
