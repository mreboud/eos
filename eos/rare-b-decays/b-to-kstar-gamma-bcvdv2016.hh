/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_GAMMA_BCVDV2016_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_GAMMA_BCVDV2016_HH 1

#include <eos/rare-b-decays/b-to-kstar-gamma-base.hh>
#include <eos/rare-b-decays/b-to-kstar-gamma-bcvdv2016.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    template <>
    class BToKstarGammaAmplitudes<tag::BCvDV2016> :
        public BToKstarGamma::AmplitudeGenerator
    {
        public:
            UsedParameter m_b_MSbar;
            UsedParameter m_s_MSbar;

            SwitchOption formfactor;
            NonlocalFormFactorPtr<nff::PToV> nonlocal_formfactor;

            BToKstarGammaAmplitudes(const Parameters & p, const Options & o);
            ~BToKstarGammaAmplitudes() = default;

            virtual BToKstarGamma::Amplitudes amplitudes() const;

            WilsonCoefficients<BToS> wilson_coefficients() const;
    };
}

#endif
