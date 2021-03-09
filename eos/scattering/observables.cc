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
                make_observable("ee->ccbar::psi2S_ee_width", R"($\Gamma(\psi(2S) \to ee)$)",
                        &EEToCCBar::psi2S_ee_width),

                make_observable("ee->ccbar::psi2S_eff_width", R"($\Gamma(\psi(2S) \to \textrm{eff})$)",
                        &EEToCCBar::psi2S_eff_width),

                make_observable("ee->ccbar::psi2S_total_width", R"($\Gamma_{\psi(2S)}$)",
                        &EEToCCBar::psi2S_total_width),

                make_cacheable_observable("ee->ccbar::sigma_eetoee(E)", R"($\sigma(ee \to ee)$)",
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoee,
                        std::make_tuple("E")),

                make_cacheable_observable("ee->ccbar::sigma_eetoeff(E)", R"($\sigma(ee \to \textrm{eff})$)",
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoeff,
                        std::make_tuple("E")),

                make_cacheable_observable("ee->ccbar::sigma_eetoD0Dbar0(E)", R"($\sigma(ee \to D^0 \bar{D}^0)$)",
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoD0Dbar0,
                        std::make_tuple("E")),

                make_cacheable_observable("ee->ccbar::sigma_eetoDpDm(E)", R"($\sigma(ee \to D^+ D^-)$)",
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDpDm,
                        std::make_tuple("E")),

                make_cacheable_observable("ee->ccbar::sigma_eetoD0Dbarst0(E)", R"($\sigma(ee \to D^0 \bar{D}^{0*})$)",
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoD0Dbarst0,
                        std::make_tuple("E")),

                make_cacheable_observable("ee->ccbar::sigma_eetoDpDstm(E)", R"($\sigma(ee \to D^+ D^{-*})$)",
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDpDstm,
                        std::make_tuple("E")),

                make_cacheable_observable("ee->ccbar::sigma_eetoDspDsm(E)", R"($\sigma(ee \to D_s^+ D_s^-)$)",
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDspDsm,
                        std::make_tuple("E")),

                make_cacheable_observable("ee->ccbar::sigma_eetoDst0Dbarst0(E)", R"($\sigma(ee \to D^{*0} \bar{D}^{*-})$)",
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDst0Dbarst0,
                        std::make_tuple("E")),

                make_cacheable_observable("ee->ccbar::sigma_eetoDstpDstm(E)", R"($\sigma(ee \to D^[*+} D^{*-})$)",
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDstpDstm,
                        std::make_tuple("E")),

                make_cacheable_observable("ee->ccbar::sigma_eetoDspDsstm(E)", R"($\sigma(ee \to D_s^+ D_s^{*-})$)",
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDspDsstm,
                        std::make_tuple("E")),

                make_cacheable_observable("ee->ccbar::sigma_eetoDsstpDsstm(E)", R"($\sigma(ee \to D_s^{*+} D_s^{*-})$)",
                        &EEToCCBar::prepare,
                        &EEToCCBar::sigma_eetoDsstpDsstm,
                        std::make_tuple("E")),

                make_cacheable_observable("ee->ccbar::R_c(E)", R"($R_c$)",
                        &EEToCCBar::prepare,
                        &EEToCCBar::Rc,
                        std::make_tuple("E")),

                make_cacheable_observable("ee->ccbar::R_udsc_prior(E)", R"($R_{udsc}$ prior)",
                        &EEToCCBar::prepare,
                        &EEToCCBar::Rudsc_prior,
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
