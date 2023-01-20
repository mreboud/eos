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

#include <eos/observable-impl.hh>
#include <eos/scattering/ee-to-ccbar.hh>
#include <eos/utils/concrete_observable.hh>
#include <eos/utils/concrete-cacheable-observable.hh>

namespace eos
{
    // ee -> ccbar
    // {{{
    ObservableGroup
    make_ee_to_ccbar_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $ee \to c\bar{c}$ processes)",
            R"(The options "nchannel" and "nresonance" fix the number of channels and resonances respectively.)",
            {
                make_observable("Jpsi->e^+e^-::decay_width", R"($\Gamma(J/\psi \to ee)$)",
                        Unit::GeV(),
                        &EEToCCBar::Jpsi_ee_width),

                make_observable("Jpsi->eff::decay_width", R"($\Gamma(J/\psi \to \textrm{eff})$)",
                        Unit::GeV(),
                        &EEToCCBar::Jpsi_eff_width),

                make_observable("psi(2S)->e^+e^-::decay_width", R"($\Gamma(\psi(2S) \to ee)$)",
                        Unit::GeV(),
                        &EEToCCBar::psi2S_ee_width),

                make_observable("psi(2S)->eff::decay_width", R"($\Gamma(\psi(2S) \to \textrm{eff})$)",
                        Unit::GeV(),
                        &EEToCCBar::psi2S_eff_width),

                make_observable("psi(4040)->DD::decay_width", R"($\Gamma(\psi(4040) \to D\bar{D})$)",
                        Unit::GeV(),
                        &EEToCCBar::psi4040_DD_width),

                make_observable("psi(4040)->DD^*::decay_width", R"($\Gamma(\psi(4040) \to D\bar{D}^*)$)",
                        Unit::GeV(),
                        &EEToCCBar::psi4040_DDst_width),

                make_observable("psi(4040)->D^*D^*::decay_width", R"($\Gamma(\psi(4040) \to D^*\bar{D}^*)$)",
                        Unit::GeV(),
                        &EEToCCBar::psi4040_DstDst_width),

                make_observable("psi(4160)->DD::decay_width", R"($\Gamma(\psi(4160) \to D\bar{D})$)",
                        Unit::GeV(),
                        &EEToCCBar::psi4160_DD_width),

                make_observable("psi(4160)->DD^*::decay_width", R"($\Gamma(\psi(4160) \to D\bar{D}^*)$)",
                        Unit::GeV(),
                        &EEToCCBar::psi4160_DDst_width),

                make_observable("psi(4160)->D^*D^*::decay_width", R"($\Gamma(\psi(4160) \to D^*\bar{D}^*)$)",
                        Unit::GeV(),
                        &EEToCCBar::psi4160_DstDst_width),

                make_observable("psi(4415)->DD::decay_width", R"($\Gamma(\psi(4415) \to D\bar{D})$)",
                        Unit::GeV(),
                        &EEToCCBar::psi4415_DD_width),

                make_observable("psi(4415)->DD^*::decay_width", R"($\Gamma(\psi(4415) \to D\bar{D}^*)$)",
                        Unit::GeV(),
                        &EEToCCBar::psi4415_DDst_width),

                make_observable("psi(4415)->D^*D^*::decay_width", R"($\Gamma(\psi(4415) \to D^*\bar{D}^*)$)",
                        Unit::GeV(),
                        &EEToCCBar::psi4415_DstDst_width),

                make_expression_observable("psi(4040)::DD_DD^*_ratio",
                        R"(\Gamma(\psi(4040) \to D\bar{D})/\Gamma(\psi(4040) \to D\bar{D}^*))",
                        Unit::None(),
                        R"(
                        <<psi(4040)->DD::decay_width>>
                        /
                        <<psi(4040)->DD^*::decay_width>>
                        )"),

                make_expression_observable("psi(4040)::D^*D^*_DD^*_ratio",
                        R"(\Gamma(\psi(4040) \to D^*\bar{D}^*)/\Gamma(\psi(4040) \to D\bar{D}^*))",
                        Unit::None(),
                        R"(
                        <<psi(4040)->D^*D^*::decay_width>>
                        /
                        <<psi(4040)->DD^*::decay_width>>
                        )"),

                make_expression_observable("psi(4160)::DD_D^*D^*_ratio",
                        R"(\Gamma(\psi(4160) \to D\bar{D})/\Gamma(\psi(4160) \to D^*\bar{D}^*))",
                        Unit::None(),
                        R"(
                        <<psi(4160)->DD::decay_width>>
                        /
                        <<psi(4160)->D^*D^*::decay_width>>
                        )"),

                make_expression_observable("psi(4160)::DD^*_D^*D^*_ratio",
                        R"(\Gamma(\psi(4160) \to D\bar{D}^*)/\Gamma(\psi(4160) \to D^*\bar{D}^*))",
                        Unit::None(),
                        R"(
                        <<psi(4160)->DD^*::decay_width>>
                        /
                        <<psi(4160)->D^*D^*::decay_width>>
                        )"),

                make_expression_observable("psi(4415)::DD_D^*D^*_ratio",
                        R"(\Gamma(\psi(4415) \to D\bar{D})/\Gamma(\psi(4415) \to D^*\bar{D}^*))",
                        Unit::None(),
                        R"(
                        <<psi(4415)->DD::decay_width>>
                        /
                        <<psi(4415)->D^*D^*::decay_width>>
                        )"),

                make_expression_observable("psi(4415)::DD^*_D^*D^*_ratio",
                        R"(\Gamma(\psi(4415) \to D\bar{D}^*)/\Gamma(\psi(4415) \to D^*\bar{D}^*))",
                        Unit::None(),
                        R"(
                        <<psi(4415)->DD^*::decay_width>>
                        /
                        <<psi(4415)->D^*D^*::decay_width>>
                        )"),

                make_observable("Jpsi::total_width", R"($\Gamma_{J/\psi}$)",
                        Unit::GeV(),
                        &EEToCCBar::Jpsi_total_width),

                make_observable("psi(2S)::total_width", R"($\Gamma_{\psi(2S)}$)",
                        Unit::GeV(),
                        &EEToCCBar::psi2S_total_width),

                make_observable("psi(3770)->D^0Dbar^0::decay_width", R"($\Gamma(\psi(3770) \to D^0\bar{D}^0)$)",
                        Unit::GeV(),
                        &EEToCCBar::psi3770_D0Dbar0_width),

                make_observable("psi(3770)->D^+D^-::decay_width", R"($\Gamma(\psi(3770) \to D^+\D^-)$)",
                        Unit::GeV(),
                        &EEToCCBar::psi3770_DpDm_width),

                make_observable("psi(3770)->eff::decay_width", R"($\Gamma(\psi(3770) \to \textrm{eff})$)",
                        Unit::GeV(),
                        &EEToCCBar::psi3770_eff_width),


                make_observable("psi(3770)::total_width", R"($\Gamma_{\psi(3770)}$)",
                        Unit::GeV(),
                        &EEToCCBar::psi3770_total_width),

                make_observable("psi(4040)::total_width", R"($\Gamma_{\psi(4040)}$)",
                        Unit::GeV(),
                        &EEToCCBar::psi4040_total_width),

                make_observable("psi(4160)::total_width", R"($\Gamma_{\psi(4160)}$)",
                        Unit::GeV(),
                        &EEToCCBar::psi4160_total_width),

                make_observable("psi(4415)::total_width", R"($\Gamma_{\psi(4415)}$)",
                        Unit::GeV(),
                        &EEToCCBar::psi4415_total_width),

                make_cacheable_observable("e^+e^-->e^+e^-::sigma(E)", R"($\sigma(e^+e^- \to e^+e^-)|_{s\textrm{-channel}}$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoee,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->eff::sigma(E)", R"($\sigma(e^+e^- \to \textrm{eff})$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoeff,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->D^0Dbar^0::sigma(E)", R"($\sigma(e^+e^- \to D^0 \bar{D}^0)$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoD0Dbar0,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->D^+D^-::sigma(E)", R"($\sigma(e^+e^- \to D^+ D^-)$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDpDm,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->D^0Dbar^*0::sigma(E)", R"($\sigma(e^+e^- \to D^0 \bar{D}^{0*})$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoD0Dbarst0,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->D^+D^*-::sigma(E)", R"($\sigma(e^+e^- \to D^+ D^{-*})$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDpDstm,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->D_s^+D_s^-::sigma(E)", R"($\sigma(e^+e^- \to D_s^+ D_s^-)$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDspDsm,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->D^*0Dbar^*0::sigma(E)", R"($\sigma(e^+e^- \to D^{*0} \bar{D}^{*-})$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDst0Dbarst0,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->D^*+D^*-::sigma(E)", R"($\sigma(e^+e^- \to D^[*+} D^{*-})$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDstpDstm,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->D_T^*+D_T^*-::sigma(E)", R"($\sigma(e^+e^- \to D^{*+} D^{*-})$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDstpTDstmT,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->D_T^*+D_L^*-::sigma(E)", R"($\sigma(e^+e^- \to D^{*+} D^{*-})$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDstpTDstmL,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->D_L^*+D_L^*-::sigma(E)", R"($\sigma(e^+e^- \to D^{*+} D^{*-})$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDstpLDstmL,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->D_s^+D_s^*-::sigma(E)", R"($\sigma(e^+e^- \to D_s^+ D_s^{*-})$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDspDsstm,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->D_s^*+D_s^*-::sigma(E)", R"($\sigma(e^+e^- \to D_s^{*+} D_s^{*-})$)",
                        Unit::InverseGeV2(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDsstpDsstm,
                        std::make_tuple("E")),

                make_cacheable_observable("e^+e^-->ccbar::R(E)", R"($R$)",
                        Unit::None(),
                        &EEToCCBar::prepare,
                        &EEToCCBar::R,
                        std::make_tuple("E")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    ObservableSection
    make_scattering_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Observables in scattering processes",
            "",
            {
                // ee -> ccbar
                make_ee_to_ccbar_group(),

            }
        );

        return ObservableSection(imp);
    }
}
