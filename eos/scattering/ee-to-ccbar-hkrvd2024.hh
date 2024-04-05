/*
 * Copyright (c) 2023 Méril Reboud
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

#ifndef EOS_GUARD_EOS_SCATTERING_EETOCCBAR_HKRvD2024_HH
#define EOS_GUARD_EOS_SCATTERING_EETOCCBAR_HKRvD2024_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/kmatrix.hh>
#include <eos/utils/concrete-cacheable-observable.hh>

#include <complex>
#include <vector>

namespace eos

{
    /*
    Channels follow the following convention
    #   name          type         Nf      copy
    0   ee            ee (S)       3       -
    1   effJpsi       PP (P)       3       -
    2   eff(2S)       PP (P)       3       -
    3   eff(3770)     PP (P)       3       -
    4   D0   D0bar    PP (P)       3       -
    5   D+   D-       PP (P)       3       4 (isospin)
    */

    class EEToCCBar :
        public ParameterUser,
        public PrivateImplementationPattern<EEToCCBar>
    {
        public:

            const static long unsigned nchannels = 6;
            const static long unsigned nresonances = 3;

            struct IntermediateResult :
                public CacheableObservable::IntermediateResult
            {
                std::shared_ptr<KMatrix<nchannels, nresonances>> K;
                // Amplitude on the first RS
                std::array<complex<double>, nchannels> tmatrix_row_0;
                // Amplitude on the second RS
                std::array<complex<double>, nchannels> tmatrix2_row_0;

                complex<double> E;
                complex<double> s;
            };

            EEToCCBar(const Parameters & parameters, const Options & options);
            ~EEToCCBar();

            // Observables
            const IntermediateResult * prepare(const double & E) const;
            const IntermediateResult * prepare_complex(const double & reE, const double & imE) const;

            // double evaluate(const IntermediateResult *) const;

            // resonances widths
            double Jpsi_ee_width() const;
            double Jpsi_eff_width() const;
            double Jpsi_total_width() const;
            double psi2S_ee_width() const;
            double psi2S_eff_width() const;
            double psi2S_total_width() const;
            double psi3770_D0Dbar0_width() const;
            double psi3770_DpDm_width() const;
            double psi3770_eff_width() const;
            double psi3770_total_width() const;

            // Phase space factor
            double rho_ee(const IntermediateResult *) const;
            double rho_eff(const IntermediateResult *) const;
            double rho_D0Dbar0(const IntermediateResult *) const;
            double rho_DpDm(const IntermediateResult *) const;

            // Chew-Mandelstam function on the first Riemann sheet
            double re_chew_mandelstam_ee(const IntermediateResult *) const;
            double im_chew_mandelstam_ee(const IntermediateResult *) const;
            double re_chew_mandelstam_eff(const IntermediateResult *) const;
            double im_chew_mandelstam_eff(const IntermediateResult *) const;
            double re_chew_mandelstam_D0Dbar0(const IntermediateResult *) const;
            double im_chew_mandelstam_D0Dbar0(const IntermediateResult *) const;
            double re_chew_mandelstam_DpDm(const IntermediateResult *) const;
            double im_chew_mandelstam_DpDm(const IntermediateResult *) const;

            // Chew-Mandelstam function on the second Riemann sheet
            double re_chew_mandelstam_II_ee(const IntermediateResult *) const;
            double im_chew_mandelstam_II_ee(const IntermediateResult *) const;
            double re_chew_mandelstam_II_eff(const IntermediateResult *) const;
            double im_chew_mandelstam_II_eff(const IntermediateResult *) const;
            double re_chew_mandelstam_II_D0Dbar0(const IntermediateResult *) const;
            double im_chew_mandelstam_II_D0Dbar0(const IntermediateResult *) const;
            double re_chew_mandelstam_II_DpDm(const IntermediateResult *) const;
            double im_chew_mandelstam_II_DpDm(const IntermediateResult *) const;


            // amplitudes on the first RS
            double re_T_eetoee(const IntermediateResult *) const;
            double im_T_eetoee(const IntermediateResult *) const;
            double re_T_eetoeff(const IntermediateResult *) const;
            double im_T_eetoeff(const IntermediateResult *) const;
            double re_T_eetoDpDm(const IntermediateResult *) const;
            double im_T_eetoDpDm(const IntermediateResult *) const;
            double re_T_eetoD0Dbar0(const IntermediateResult *) const;
            double im_T_eetoD0Dbar0(const IntermediateResult *) const;

            // amplitudes on the second RS
            double re_T_II_eetoee(const IntermediateResult *) const;
            double im_T_II_eetoee(const IntermediateResult *) const;
            double re_T_II_eetoeff(const IntermediateResult *) const;
            double im_T_II_eetoeff(const IntermediateResult *) const;
            double re_T_II_eetoDpDm(const IntermediateResult *) const;
            double im_T_II_eetoDpDm(const IntermediateResult *) const;
            double re_T_II_eetoD0Dbar0(const IntermediateResult *) const;
            double im_T_II_eetoD0Dbar0(const IntermediateResult *) const;

            // Spectral function
            double psi3770_spectral_function(const double & E) const;

            // sigma(ee -> channel)
            double sigma_eetoee(const IntermediateResult *) const;
            double sigma_eetoeff(const IntermediateResult *) const;
            double sigma_eetoD0Dbar0(const IntermediateResult *) const;
            double sigma_eetoDpDm(const IntermediateResult *) const;

            // R ratios
            double R(const IntermediateResult *) const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };
}
#endif
