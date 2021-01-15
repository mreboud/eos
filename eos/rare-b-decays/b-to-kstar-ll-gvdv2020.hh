/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_GVDV2020_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_GVDV2020_HH 1

#include <eos/rare-b-decays/b-to-kstar-ll-base.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    template <>
    class BToKstarDileptonAmplitudes<tag::GvDV2020> :
        public BToKstarDilepton::AmplitudeGenerator
    {
        public:
            UsedParameter hbar;

            UsedParameter m_b_MSbar;
            UsedParameter m_s_MSbar;

            UsedParameter mu;
            UsedParameter alpha_e;
            UsedParameter g_fermi;
            UsedParameter tau;

            SwitchOption q;

            SwitchOption formfactor;
            NonlocalFormFactorPtr<nc::PToV> nonlocal_formfactor;

            BToKstarDileptonAmplitudes(const Parameters & p, const Options & o);
            ~BToKstarDileptonAmplitudes() = default;

            virtual BToKstarDilepton::Amplitudes amplitudes(const double & q2) const;

            WilsonCoefficients<BToS> wilson_coefficients() const;
    };
}

#endif
