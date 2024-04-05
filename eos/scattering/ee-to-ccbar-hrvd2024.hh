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

#ifndef EOS_GUARD_EOS_SCATTERING_EETOCCBAR_HRvD2024_HH
#define EOS_GUARD_EOS_SCATTERING_EETOCCBAR_HRvD2024_HH 1

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
    1   eff(2S)       PP (P)       3       -
    2   D0   D0bar    PP (P)       3       -
    3   D+   D-       PP (P)       3       2  (isospin option)
    4   eff(3770)     PP (P)       3       -
    5   D0   D*0bar   VP (P)       3       -
    6   D*0  D0bar    VP (P)       3       5  (c.c.)
    7   D+   D*-      VP (P)       3       5  (isospin option)
    8   D*+  D-       VP (P)       3       7  (c.c.)
    9   Ds+  Ds-      PP (P)       3       2  (u-spin option)
    10  D*0  D*0bar   VV (P, S=0)  3       -
    11  D*0  D*0bar   VV (P, S=2)  3       -
    12  D*0  D*0bar   VV (F, S=2)  7       -
    13  D*+  D*-      VV (P, S=0)  3       10 (isospin option)
    14  D*+  D*-      VV (P, S=2)  3       11 (isospin option)
    15  D*+  D*-      VV (F, S=2)  7       12 (isospin option)
    16  eff(4040)     PP (P)       3       -
    17  Ds+  Ds*-     VP (P)       3       5  (u-spin option)
    18  Ds*+ Ds-      VP (P)       3       17 (c.c.)
    19  eff(4160)     PP (P)       3       -
    // 20  Ds*+ Ds*-     VV (P, S=0)  3       10 (u-spin option)
    // 21  Ds*+ Ds*-     VV (P, S=2)  3       11 (u-spin option)
    // 22  Ds*+ Ds*-     VV (F, S=2)  7       12 (u-spin option)
    */

    class EEToCCBarHRvD :
        public ParameterUser,
        public PrivateImplementationPattern<EEToCCBarHRvD>
    {
        public:

            const static long unsigned nchannels = 20;
            const static long unsigned nresonances = 4;

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

            EEToCCBarHRvD(const Parameters & parameters, const Options & options);
            ~EEToCCBarHRvD();

            // Observables
            const IntermediateResult * prepare(const double & E) const;
            const IntermediateResult * prepare_complex(const double & reE, const double & imE) const;

            // double evaluate(const IntermediateResult *) const;

            // resonances widths
            double psi2S_ee_width() const;
            double psi2S_eff_width() const;
            double psi2S_total_width() const;
            double psi3770_D0Dbar0_width() const;
            double psi3770_DpDm_width() const;
            double psi3770_eff_width() const;
            double psi3770_total_width() const;
            double psi4040_D0Dbar0_width() const;
            double psi4040_DpDm_width() const;
            double psi4040_D0Dbarst0_width() const;
            double psi4040_DpDstm_width() const;
            double psi4040_DspDsm_width() const;
            double psi4040_Dst0Dbarst0_width() const;
            double psi4040_DstpDstm_width() const;
            double psi4040_eff_width() const;
            double psi4040_total_width() const;
            double psi4160_D0Dbar0_width() const;
            double psi4160_DpDm_width() const;
            double psi4160_D0Dbarst0_width() const;
            double psi4160_DpDstm_width() const;
            double psi4160_DspDsm_width() const;
            double psi4160_Dst0Dbarst0_width() const;
            double psi4160_DstpDstm_width() const;
            double psi4160_DspDsstm_width() const;
            double psi4160_eff_width() const;
            double psi4160_total_width() const;

            // Phase space factor
            double rho_ee(const IntermediateResult *) const;
            double rho_eff(const IntermediateResult *) const;
            double rho_D0Dbar0(const IntermediateResult *) const;
            double rho_DpDm(const IntermediateResult *) const;
            double rho_D0Dbarst0(const IntermediateResult *) const;
            double rho_DpDstm(const IntermediateResult *) const;
            double rho_DspDsm(const IntermediateResult *) const;
            double rho_Dst0Dbarst0(const IntermediateResult *) const;
            double rho_DstpDstm(const IntermediateResult *) const;
            double rho_DspDsstm(const IntermediateResult *) const;

            // Chew-Mandelstam function on the first Riemann sheet
            double re_chew_mandelstam_ee(const IntermediateResult *) const;
            double im_chew_mandelstam_ee(const IntermediateResult *) const;
            double re_chew_mandelstam_eff(const IntermediateResult *) const;
            double im_chew_mandelstam_eff(const IntermediateResult *) const;
            double re_chew_mandelstam_D0Dbar0(const IntermediateResult *) const;
            double im_chew_mandelstam_D0Dbar0(const IntermediateResult *) const;
            double re_chew_mandelstam_DpDm(const IntermediateResult *) const;
            double im_chew_mandelstam_DpDm(const IntermediateResult *) const;
            double re_chew_mandelstam_D0Dbarst0(const IntermediateResult *) const;
            double im_chew_mandelstam_D0Dbarst0(const IntermediateResult *) const;
            double re_chew_mandelstam_DpDstm(const IntermediateResult *) const;
            double im_chew_mandelstam_DpDstm(const IntermediateResult *) const;
            double re_chew_mandelstam_DspDsm(const IntermediateResult *) const;
            double im_chew_mandelstam_DspDsm(const IntermediateResult *) const;
            double re_chew_mandelstam_Dst0Dbarst0(const IntermediateResult *) const;
            double im_chew_mandelstam_Dst0Dbarst0(const IntermediateResult *) const;
            double re_chew_mandelstam_DstpDstm(const IntermediateResult *) const;
            double im_chew_mandelstam_DstpDstm(const IntermediateResult *) const;
            double re_chew_mandelstam_DspDsstm(const IntermediateResult *) const;
            double im_chew_mandelstam_DspDsstm(const IntermediateResult *) const;

            // Chew-Mandelstam function on the second Riemann sheet
            double re_chew_mandelstam_II_ee(const IntermediateResult *) const;
            double im_chew_mandelstam_II_ee(const IntermediateResult *) const;
            double re_chew_mandelstam_II_eff(const IntermediateResult *) const;
            double im_chew_mandelstam_II_eff(const IntermediateResult *) const;
            double re_chew_mandelstam_II_D0Dbar0(const IntermediateResult *) const;
            double im_chew_mandelstam_II_D0Dbar0(const IntermediateResult *) const;
            double re_chew_mandelstam_II_DpDm(const IntermediateResult *) const;
            double im_chew_mandelstam_II_DpDm(const IntermediateResult *) const;
            double re_chew_mandelstam_II_D0Dbarst0(const IntermediateResult *) const;
            double im_chew_mandelstam_II_D0Dbarst0(const IntermediateResult *) const;
            double re_chew_mandelstam_II_DpDstm(const IntermediateResult *) const;
            double im_chew_mandelstam_II_DpDstm(const IntermediateResult *) const;
            double re_chew_mandelstam_II_DspDsm(const IntermediateResult *) const;
            double im_chew_mandelstam_II_DspDsm(const IntermediateResult *) const;
            double re_chew_mandelstam_II_Dst0Dbarst0(const IntermediateResult *) const;
            double im_chew_mandelstam_II_Dst0Dbarst0(const IntermediateResult *) const;
            double re_chew_mandelstam_II_DstpDstm(const IntermediateResult *) const;
            double im_chew_mandelstam_II_DstpDstm(const IntermediateResult *) const;
            double re_chew_mandelstam_II_DspDsstm(const IntermediateResult *) const;
            double im_chew_mandelstam_II_DspDsstm(const IntermediateResult *) const;


            // amplitudes on the first RS
            double re_T_eetoee(const IntermediateResult *) const;
            double im_T_eetoee(const IntermediateResult *) const;
            double re_T_eetoeff(const IntermediateResult *) const;
            double im_T_eetoeff(const IntermediateResult *) const;
            double re_T_eetoD0Dbar0(const IntermediateResult *) const;
            double im_T_eetoD0Dbar0(const IntermediateResult *) const;
            double re_T_eetoDpDm(const IntermediateResult *) const;
            double im_T_eetoDpDm(const IntermediateResult *) const;
            double re_T_eetoD0Dbarst0(const IntermediateResult *) const;
            double im_T_eetoD0Dbarst0(const IntermediateResult *) const;
            double re_T_eetoDpDstm(const IntermediateResult *) const;
            double im_T_eetoDpDstm(const IntermediateResult *) const;
            double re_T_eetoDspDsm(const IntermediateResult *) const;
            double im_T_eetoDspDsm(const IntermediateResult *) const;
            double re_T_eetoDst0Dbarst0(const IntermediateResult *) const;
            double im_T_eetoDst0Dbarst0(const IntermediateResult *) const;
            double re_T_eetoDstpDstm(const IntermediateResult *) const;
            double im_T_eetoDstpDstm(const IntermediateResult *) const;
            double re_T_eetoDspDsstm(const IntermediateResult *) const;
            double im_T_eetoDspDsstm(const IntermediateResult *) const;

            // amplitudes on the second RS
            double re_T_II_eetoee(const IntermediateResult *) const;
            double im_T_II_eetoee(const IntermediateResult *) const;
            double re_T_II_eetoeff(const IntermediateResult *) const;
            double im_T_II_eetoeff(const IntermediateResult *) const;
            double re_T_II_eetoD0Dbar0(const IntermediateResult *) const;
            double im_T_II_eetoD0Dbar0(const IntermediateResult *) const;
            double re_T_II_eetoDpDm(const IntermediateResult *) const;
            double im_T_II_eetoDpDm(const IntermediateResult *) const;
            double re_T_II_eetoD0Dbarst0(const IntermediateResult *) const;
            double im_T_II_eetoD0Dbarst0(const IntermediateResult *) const;
            double re_T_II_eetoDpDstm(const IntermediateResult *) const;
            double im_T_II_eetoDpDstm(const IntermediateResult *) const;
            double re_T_II_eetoDspDsm(const IntermediateResult *) const;
            double im_T_II_eetoDspDsm(const IntermediateResult *) const;
            double re_T_II_eetoDst0Dbarst0(const IntermediateResult *) const;
            double im_T_II_eetoDst0Dbarst0(const IntermediateResult *) const;
            double re_T_II_eetoDstpDstm(const IntermediateResult *) const;
            double im_T_II_eetoDstpDstm(const IntermediateResult *) const;
            double re_T_II_eetoDspDsstm(const IntermediateResult *) const;
            double im_T_II_eetoDspDsstm(const IntermediateResult *) const;

            // Spectral function
            double psi3770_spectral_function(const double & E) const;
            double psi4040_spectral_function(const double & E) const;
            double psi4160_spectral_function(const double & E) const;

            // sigma(ee -> channel)
            double sigma_eetoee(const IntermediateResult *) const;
            double sigma_eetoeff(const IntermediateResult *) const;
            double sigma_eetoD0Dbar0(const IntermediateResult *) const;
            double sigma_eetoDpDm(const IntermediateResult *) const;
            double sigma_eetoD0Dbarst0(const IntermediateResult *) const;
            double sigma_eetoDpDstm(const IntermediateResult *) const;
            double sigma_eetoDspDsm(const IntermediateResult *) const;
            double sigma_eetoDst0Dbarst0(const IntermediateResult *) const;
            double sigma_eetoDstpDstm(const IntermediateResult *) const;
            double sigma_eetoDspDsstm(const IntermediateResult *) const;

            // R ratio
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
