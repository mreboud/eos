/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <eos/rare-b-decays/b-to-kstar-ll-base.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>

namespace eos
{
    BToKstarDilepton::AmplitudeGenerator::AmplitudeGenerator(const Parameters & p, const Options & o) :
        model(Model::make(o.get("model", "SM"), p, o)),
        form_factors(FormFactorFactory<PToV>::create("B->K^*::" + o.get("form-factors", "KMPW2010"), p)),
        mu(p["mu"], *this),
        alpha_e(p["QED::alpha_e(m_b)"], *this),
        g_fermi(p["G_Fermi"], *this),
        tau(p["life_time::B_" + o.get("q", "d")], *this),
        m_B(p["mass::B_" + o.get("q", "d")], *this),
        m_Kstar(p["mass::K_d^*"], *this),
        m_l(p["mass::" + o.get("l", "mu")], *this),
        cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false"))),
        lepton_flavour(o.get("l", "mu"))
    {
        if (! form_factors.get())
            throw InternalError("Form factors not found!");

        if (0.0 == m_l())
        {
            throw InternalError("Zero lepton mass leads to NaNs in timelike amplitudes. Use tiny lepton mass > 0!");
        }

        this->uses(*form_factors);
        this->uses(*model);
    }

    BToKstarDilepton::AmplitudeGenerator::~AmplitudeGenerator()
    {
    }

    double
    BToKstarDilepton::AmplitudeGenerator::beta_l(const double & s) const
    {
        return std::sqrt(1.0 - 4.0 * m_l * m_l / s);
    }

    double
    BToKstarDilepton::AmplitudeGenerator::lambda(const double & s) const
    {
        return eos::lambda(m_B() * m_B(), m_Kstar() * m_Kstar(), s);
    }

    double
    BToKstarDilepton::AmplitudeGenerator::energy(const double & s) const
    {
        return (m_B() * m_B() + m_Kstar() * m_Kstar() - s) / (2.0 * m_B());
    }

    double
    BToKstarDilepton::AmplitudeGenerator::s_hat(const double & s) const
    {
        return s / m_B() / m_B();
    }

}
