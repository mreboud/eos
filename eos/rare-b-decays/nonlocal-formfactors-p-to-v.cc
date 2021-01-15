/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2019 Danny van Dyk
 * Copyright (c) 2019 Nico Gubernari
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

#include <eos/rare-b-decays/nonlocal-formfactors.hh>
#include <eos/rare-b-decays/charm-loops.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>
#include <eos/form-factors/b-lcdas.hh>
#include <eos/form-factors/kstar-lcdas.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/integrate-impl.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/kinematic.hh>

#include <cmath>
#include <map>

#include <gsl/gsl_sf.h>

namespace eos
{
    namespace nc_p_to_v
    {
        class Naive :
            public NonlocalFormFactor<nc::PToV>
        {
            public:
                Naive(const Parameters &, const Options &)
                {
                }

                ~Naive() = default;

                virtual complex<double> H_perp(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> H_para(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> H_long(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> Hhat_perp(const double & q2) const
                {
                    return 0.0;
                }
                virtual complex<double> Hhat_para(const double & q2) const
                {
                    return 0.0;
                }
                virtual complex<double> Hhat_long(const double & q2) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V1(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V2(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V23(const double &) const
                {
                    return 0.0;
                }

                static NonlocalFormFactorPtr<nc::PToV> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nc::PToV>(new Naive(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    return {};
                }
        };

        /*
         * Calculates the QCD-improved factorisation results for the formfactor
         * as given in [BFS2001] and [BFS2004].
         */
        class QCDF :
            public NonlocalFormFactor<nc::PToV>
        {
            public:
                std::shared_ptr<Model> model;
                std::shared_ptr<FormFactors<PToV>> form_factors;
                KstarLCDAs kstar;

                UsedParameter mu_f;

                SwitchOption q;

                UsedParameter m_B;
                UsedParameter m_Kstar;

                UsedParameter f_B;
                UsedParameter lambda_B_p;

                UsedParameter mu;

                double e_q;

                std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                        const double &, const double &, const double &, const double &,
                        const double &, const double &)> qcdf_dilepton_massless_case;
                std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                        const double &, const double &, const double &, const double &,
                        const double &, const double &, const double &)> qcdf_dilepton_charm_case;
                std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                        const double &, const double &, const double &, const double &,
                        const double &, const double &, const double &)> qcdf_dilepton_bottom_case;

                QCDF(const Parameters & p, const Options & o) :
                    model(Model::make(o.get("model", "SM"), p, o)),
                    form_factors(FormFactorFactory<PToV>::create("B->K^*@" + o.get("form-factors", "KMPW2010"), p)),
                    kstar(p, o),
                    mu_f(p["B->K^*ccbar::mu_f"], *this),
                    q(o, "q", { "d", "u" }),
                    m_B(p["mass::B_" + q.value()], *this),
                    m_Kstar(p["mass::K^*_" + q.value()], *this),
                    f_B(p["decay-constant::B_" + o.get("q", "d")], *this),
                    lambda_B_p(p["lambda_B_p"], *this),
                    e_q((q.value() == "d") ? -1.0 / 3.0 : +2.0 / 3.0),
                    mu(p["mu"], *this)
                {
                    // Select the appropriate calculator for the QCDF integrals
                    std::string qcdf_integrals(o.get("qcdf-integrals", "mixed"));
                    using std::placeholders::_1;
                    using std::placeholders::_2;
                    using std::placeholders::_3;
                    using std::placeholders::_4;
                    using std::placeholders::_5;
                    using std::placeholders::_6;
                    using std::placeholders::_7;
                    using std::placeholders::_8;
                    using std::placeholders::_9;
                    if ("mixed" == qcdf_integrals)
                    {
                        qcdf_dilepton_massless_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::dilepton_massless_case,
                                    _1, _2, _3, _4, _5, _6, _7, _8);
                        qcdf_dilepton_charm_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::dilepton_charm_case,
                                    _1, _2, _3, _4, _5, _6, _7, _8, _9);
                        qcdf_dilepton_bottom_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::dilepton_bottom_case,
                                    _1, _2, _3, _4, _5, _6, _7, _8, _9);
                    }
                    else if ("numerical" == qcdf_integrals)
                    {
                        qcdf_dilepton_massless_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::dilepton_massless_case,
                                    _1, _2, _3, _4, _5, _6, _7, _8);
                        qcdf_dilepton_charm_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::dilepton_charm_case,
                                    _1, _2, _3, _4, _5, _6, _7, _8, _9);
                        qcdf_dilepton_bottom_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::dilepton_bottom_case,
                                    _1, _2, _3, _4, _5, _6, _7, _8, _9);
                    }
                    else if ("analytical" == qcdf_integrals)
                    {
                        qcdf_dilepton_massless_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_massless_case,
                                    _1, _2, _3, _4, _5, _6, _7, _8);
                        qcdf_dilepton_charm_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_charm_case,
                                    _1, _2, _3, _4, _5, _6, _7, _8, _9);
                        qcdf_dilepton_bottom_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_bottom_case,
                                    _1, _2, _3, _4, _5, _6, _7, _8, _9);
                    }
                    else
                    {
                        throw InvalidOptionValueError("qcdf-integrals", qcdf_integrals, "mixed, numerical, analytical");
                    }

                    this->uses(*form_factors);
                    this->uses(*model);
                }

                ~QCDF() = default;

                // scale at which the WC C9/C7 would be evaluated
                inline double mu_b() const
                {
                    return 4.2;
                }

                // the functions F_{1,2}^{7,9} expect m_c to be in the pole mass scheme
                inline double m_c() const
                {
                    return model->m_c_pole();
                }

                // cf. eq. (15), [BFS2001], C_1 and C_2 contributions only
                std::tuple<complex<double>, complex<double>>
                calT(const double & s ) const
                {
                    // charges of down- and up-type quarks
                    //static const double e_d = -1.0/3.0;
                    static const double e_u = +2.0/3.0;

                    // Wilson coefficients
                    const WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), "mu");

                    // spectator contributions
                    //double delta_qu = (q.value() == "u" ? 1.0 : 0.0);

                    // quark masses and kinematics
                    double m_c    = this->m_c();
                    double m_b_PS = model->m_b_ps(1.5);//, m_b_PS2 = m_b_PS * m_b_PS;
                    double energy = (m_B * m_B + m_Kstar * m_Kstar - s) / (2.0 * m_B);
                    //double L = -1.0 * (m_b_PS2 - s) / s * std::log(1.0 - s / m_b_PS2);

                    // couplings
                    double alpha_s_mu = model->alpha_s(mu_f()); // alpha_s at the hard scale
                    double a_mu = alpha_s_mu * QCD::casimir_f / 4.0 / M_PI;
                    double alpha_s_mu_f = model->alpha_s(std::sqrt(mu_f() * 0.5)); // alpha_s at the factorization scale
                    double a_mu_f = alpha_s_mu_f * QCD::casimir_f / 4.0 / M_PI;

                    // Compute the QCDF Integrals
                    //double invm1_par = 3.0 * (1.0 + a_1_para + a_2_para); // <ubar^-1>_par
                    //double invm1_perp = 3.0 * (1.0 + a_1_perp + a_2_perp); // <ubar^-1>_perp
                    const double a_1_perp = kstar.a_1_perp(mu_f());
                    const double a_2_perp = kstar.a_2_perp(mu_f());
                    const double a_1_para = kstar.a_1_para(mu_f());
                    const double a_2_para = kstar.a_2_para(mu_f());
                    QCDFIntegrals<BToKstarDilepton> qcdf_c = this->qcdf_dilepton_charm_case(s, m_c, m_B, m_Kstar, mu_b(), a_1_perp, a_2_perp, a_1_para, a_2_para);

                    // inverse of the "negative" moment of the B meson LCDA
                    // cf. [BFS2001], Eq. (54), p. 15
                    double omega_0 = lambda_B_p, lambda_B_p_inv = 1.0 / lambda_B_p;
                    complex<double> lambda_B_m_inv = std::numeric_limits<double>::quiet_NaN();
                    if (abs(s) > 1e-4)
                    {
                        lambda_B_m_inv = complex<double>(-gsl_sf_expint_Ei(s / m_B / omega_0), M_PI) * (std::exp(-s / m_B / omega_0) / omega_0);
                    }

                    /* Y(s) for the up and the top sector */
                    // cf. [BFS2001], Eq. (10), p. 4
                    complex<double> Y_top_c = 4.0 / 3.0 * wc.c1() + wc.c2();

                    // Use b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses,
                    // then replace b pole mass by the PS mass.
                    complex<double> Y_top = Y_top_c * CharmLoops::h(mu_b(), s, m_c);

                    // Charm-loop functions at order alpha_s
                    complex<double> F27 = memoise(CharmLoops::F27_massive, mu_b(), s, m_b_PS, m_c);
                    complex<double> F19 = 0.0;
                    complex<double> F29 = 0.0;
                    if (abs(s) >= 1e-4) // only evaluate F19 and F29 if s != 0
                    {
                        F19 = memoise(CharmLoops::F19_massive, mu_b(), s, m_b_PS, m_c);
                        F29 = memoise(CharmLoops::F29_massive, mu_b(), s, m_b_PS, m_c);
                    }

                    /* perpendicular, top sector */
                    // cf. [BFS2001], Eqs. (12), (15), p. 5, in comparison with \delta_1 = 1
                    complex<double> C0_top_perp = s / (2.0 * m_b_PS * m_B) * Y_top;
                    // cf. [BFS2001], Eqs. (34), (37), p. 9
                    complex<double> C1nf_top_perp = (-1.0 / QCD::casimir_f) * (
                            (wc.c2() - wc.c1() / 6.0) * F27
                            + (s / (2.0 * m_b_PS * m_B)) * (wc.c1() * F19 + wc.c2() * F29));

                    /* parallel, top sector */
                    // cf. [BFS2001], Eqs. (14), (15), p. 5, in comparison with \delta_{2,3} = 1
                    complex<double> C0_top_para   = -1.0 * (m_B / (2.0 * m_b_PS) * Y_top);
                    // cf. [BFS2001], Eqs. (38), p. 9
                    complex<double> C1nf_top_para = (+1.0 / QCD::casimir_f) * (
                            (wc.c2() - wc.c1() / 6.0) * F27
                            + (m_B / (2.0 * m_b_PS)) * (wc.c1() * F19 + wc.c2() * F29));

                    // compute the factorizing contributions
                    complex<double> C_perp = C0_top_perp + a_mu * C1nf_top_perp;
                    complex<double> C_para = C0_top_para + a_mu * C1nf_top_para;

                    /* perpendicular, top sector */
                    // T0_top_perp_{p,m} = 0, cf. [BFS2001], Eq. (17), p. 6
                    // cf. [BFS2001], Eq. (23), p. 7
                    complex<double> T1nf_top_perp_p = m_B / (2.0 * m_b_PS) * e_u * (-wc.c1() / 6.0 + wc.c2()) * qcdf_c.jtilde1_perp * lambda_B_p_inv;
                    // T1nf_top_perp_m = 0, cf. [BFS2001], Eq. (17), p. 6

                    /* parallel, top sector */
                    // T0_top_par_p = 0, cf. [BFS2001], Eq. (17), p. 6
                    // cf. [BFS2004], Eqs. (46)-(47), p. 25 without the \omega term.
                    complex<double> T0_top_para_m = -e_q * 4.0 * m_B / m_b_PS * (wc.c3() + 4.0/3.0 * wc.c4() + 16.0 * wc.c5() + 64.0/3.0 * wc.c6()) * lambda_B_m_inv;
                    // cf. [BFS2001], Eq. (25), p. 7
                    complex<double> T1nf_top_para_p = m_B / m_b_PS * e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde2_parallel * lambda_B_p_inv;
                    // cf. [BFS2001], Eq. (26), pp. 7-8
                    complex<double> T1nf_top_para_m = e_q * 6.0 * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2() + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j4_parallel * lambda_B_m_inv;

                    // Compute the nonfactorizing contributions
                    complex<double> T_perp = a_mu_f * T1nf_top_perp_p;
                    complex<double> T_para = a_mu_f * (T1nf_top_para_p + T1nf_top_para_m) + (T0_top_para_m);

                    // Compute the numerically leading power-suppressed hard spectator interaction contributions to order alpha_s^1
                    // cf. [BFS2004], Eqs. (52), (53)
                    complex<double> Delta_T_hsa_top_perp = e_q * a_mu_f * (M_PI * M_PI * f_B / (3.0 * m_b_PS * m_B)) * (
                            + 8.0 * kstar.f_perp(mu_f()) * (3.0 / 4.0) * (
                                (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j5_perp
                            )
                            - (4.0 * m_Kstar * kstar.f_para() / (1.0 - s / (m_B * m_B)) / lambda_B_p) * (3.0 / 4.0) * (
                                (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j6_perp
                            ));
                    // Compute the sum of the numerically leading power-suppressed contributions
                    complex<double> Delta_T_perp = Delta_T_hsa_top_perp;

                    const double xi_perp = m_B / (m_B + m_Kstar) * form_factors->v(s);
                    const double xi_para = (m_B + m_Kstar) / (2.0 * energy) * form_factors->a_1(s) - (1.0 - m_Kstar / m_B) * form_factors->a_2(s);

                    // cf. [BFS2001], Eq. (15), and [BHP2008], Eq. (C.4)
                    auto result = std::make_tuple(
                        xi_perp * C_perp
                            + pow(M_PI, 2) / 3.0 * (f_B * kstar.f_perp(mu_f)) / m_B * T_perp
                            + Delta_T_perp,
                        xi_para * C_para
                            + pow(M_PI, 2) / 3.0 * (f_B * kstar.f_para() * m_Kstar) / (m_B * energy) * T_para
                            );

                    return result;
                }

                virtual complex<double> H_perp(const double & s) const
                {
                    const double m_b = model->m_b_msbar(mu_b());
                    const double m_B = this->m_B(), m_B2 = m_B * m_B, m_B3 = m_B * m_B2;
                    const double m_Kstar = this->m_Kstar(), m_Kstar2 = m_Kstar * m_Kstar;
                    const double sqrt_lambda = sqrt(eos::lambda(m_B2, m_Kstar2, s));

                    auto calT = this->calT(s);
                    const complex<double> T_1 = std::get<0>(calT); // T_perp

                    return -m_b * sqrt_lambda / (8.0 * sqrt(2.0) * M_PI * M_PI * m_B3)
                            * T_1;
                }

                virtual complex<double> H_para(const double & s) const
                {
                    const double m_b = model->m_b_msbar(mu_b());
                    const double m_B = this->m_B(), m_B2 = m_B * m_B, m_B3 = m_B * m_B2;
                    const double m_Kstar = this->m_Kstar(), m_Kstar2 = m_Kstar * m_Kstar;

                    auto calT = this->calT(s);
                    const complex<double> T_2 = (m_B2 + m_Kstar2 - s) / m_B2 * std::get<0>(calT);

                    return -m_b * (m_B2 - m_Kstar2) / (8.0 * sqrt(2.0) * M_PI * M_PI * m_B3)
                            * T_2;
                }

                // Htilde_long = H_long * sqrt(q2) / M_K^*
                virtual complex<double> Htilde_long(const double & s) const
                {
                    const double m_b = model->m_b_msbar(mu_b());
                    const double m_B = this->m_B(), m_B2 = m_B * m_B;
                    const double m_Kstar = this->m_Kstar(), m_Kstar2 = m_Kstar * m_Kstar;
                    const double lambda = eos::lambda(m_B2, m_Kstar2, s);

                    auto calT = this->calT(s);
                    const complex<double> T_2 = (m_B2 + m_Kstar2 - s) / m_B2 * std::get<0>(calT); // T_perp
                    const complex<double> T_3 = std::get<1>(calT) + std::get<0>(calT); // T_para + T_perp;

                    return -m_b * s / (32 * M_PI * M_PI * m_B2 * m_Kstar2 * (m_B + m_Kstar))
                            * ((m_B2 + 3.0 * m_Kstar2 - s) * T_2 - lambda / (m_B2 - m_Kstar2) * T_3);
                }

                virtual complex<double> H_long(const double & s) const
                {
                    return Htilde_long(s) / sqrt(s) * m_Kstar();
                }


                virtual complex<double> Hhat_perp(const double & q2) const
                {
                    return 0.0;
                }
                virtual complex<double> Hhat_para(const double & q2) const
                {
                    return 0.0;
                }
                virtual complex<double> Hhat_long(const double & q2) const
                {
                    return 0.0;
                }

                static NonlocalFormFactorPtr<nc::PToV> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nc::PToV>(new QCDF(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    return {};
                }
        };

        class LCSR :
            public NonlocalFormFactor<nc::PToV>
        {
            private:
                std::shared_ptr<Model> model;

                // spectator quark option
                SwitchOption opt_q;

                // B-meson parameters
                UsedParameter m_B;
                UsedParameter f_B;

                // final state meson parameters
                UsedParameter m_V;
                UsedParameter f_V;

                // sum rule parameters
                UsedParameter M2;
                UsedParameter s0_0_V1;
                UsedParameter s0_1_V1;
                UsedParameter s0_0_V2;
                UsedParameter s0_1_V2;
                UsedParameter s0_0_V23;
                UsedParameter s0_1_V23;

                // properties of the virtual quark in the sum rule's formfactor
                std::function<double ()> m_v;

                // renormalization scale for the sum rule
                UsedParameter mu_sr;

                // renormalization scale for the Wilson coefficients
                UsedParameter mu_ren;

                // renormalization scale for the virtual charm quark's mass
                UsedParameter mu_c;

                BMesonLCDAs b_lcdas;

                std::string _final_state() const
                {
                    switch (opt_q.value()[0])
                    {
                        case 's':
                            return "phi";
                            break;

                        default:
                            return "K_d^*";
                    }
                }

                std::string _transition() const
                {
                    switch (opt_q.value()[0])
                    {
                        case 's':
                            return "B_s->phi";
                            break;

                        default:
                            return "B->K^*";
                    }
                }

            public:
                LCSR(const Parameters & p, const Options & o) :
                    model(Model::make(o.get("model", "SM"), p, o)),
                    opt_q(o, "q", { "d", "s" }, "d"),
                    m_B(p["mass::B_" + opt_q.value()], *this),
                    f_B(p["decay-constant::B_" + opt_q.value()], *this),
                    m_V(p["mass::" + _final_state()], *this),
                    f_V(p["decay-constant::"  + _final_state()], *this),
                    M2(p[_transition() + "::M^2@B-LCSR"], *this),
                    s0_0_V1(p[_transition() + "::s_0^V,0@B-LCSR"], *this),
                    s0_1_V1(p[_transition() + "::s_0^V,1@B-LCSR"], *this),
                    s0_0_V2(p[_transition() + "::s_0^A1,0@B-LCSR"], *this),
                    s0_1_V2(p[_transition() + "::s_0^A1,1@B-LCSR"], *this),
                    s0_0_V23(p[_transition() + "::s_0^A2,0@B-LCSR"], *this),
                    s0_1_V23(p[_transition() + "::s_0^A2,1@B-LCSR"], *this),
                    mu_sr(p[_transition() + "::mu@B-LCSR"], *this),
                    mu_ren(p["b->sccbar::mu"], *this),
                    mu_c(p["b->sccbar::mu_c"], *this),
                    b_lcdas(p, o + Options{ { "q", opt_q.value() } })
                {
                  this->uses(b_lcdas);

                  m_v = std::bind(&LCSR::m_virtual_s, this);
                }

                ~LCSR() = default;

                // mass of the virtual quark in the LCSR setup; v = s
                double m_virtual_s() const
                {
                    return model->m_s_msbar(mu_ren());
                }

                // mass of the virtual charm quark
                double m_c() const
                {
                    return model->m_c_msbar(mu_c());
                }

                /* forwarding the LCDAs */
                inline
                double phi_3(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.phi_3(omega_1, omega_2);
                }

                inline
                double phi_bar_3(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.phi_bar_3(omega_1, omega_2);
                }

                inline
                double phi_double_bar_3(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.phi_double_bar_3(omega_1, omega_2);
                }

                inline
                double phi_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.phi_4(omega_1, omega_2);
                }

                inline
                double phi_bar_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.phi_bar_4(omega_1, omega_2);
                }

                inline
                double phi_double_bar_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.phi_double_bar_4(omega_1, omega_2);
                }

                inline
                double psi_bar_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.psi_bar_4(omega_1, omega_2);
                }

                inline
                double psi_double_bar_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.psi_double_bar_4(omega_1, omega_2);
                }

                inline
                double chi_bar_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.chi_bar_4(omega_1, omega_2);
                }

                inline
                double chi_double_bar_4(const double omega_1, const double omega_2) const
                {
                    return b_lcdas.chi_double_bar_4(omega_1, omega_2);
                }

                /* auxilliary functions */

                double s(const double & sigma, const double & q2) const
                {
                    const double sigmabar = 1.0 - sigma;

                    return sigma * pow(m_B(), 2) + (pow(m_v(), 2) - sigma * q2) / sigmabar;
                }

                double sigma(const double & s, const double & q2) const
                {
                    const double m_B2 = pow(m_B(), 2);
                    const double m_v2 = pow(m_v(), 2);

                    return (m_B2 - q2 + s - std::sqrt(4.0 * (m_v2 - s) * m_B2 + pow(m_B2 - q2 + s, 2))) / (2.0 * m_B2);
                }

                double sigma_0(const double & q2, const double & s0_0, const double & s0_1) const
                {
                    const double s0 = s0_0 + s0_1 * q2;

                    return sigma(s0, q2);
                }

                /* Integrands */

                inline
                double I1_V1_phi_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_3 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_3 = this->phi_3(omega_1, omega_2);

                    const double C_1 = ((-(4.0*m_B*m_v*t*tbar*u) - q2*(3.0 - 2.0*(3.0 + 8.0*t2 - 8.0*t)*u)*pow(sigmabar,-1) +
                                       4.0*u2*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u)*pow(omega,2.0)*pow(sigmabar,-1) +
                                       4.0*m_v*(m_v2 - q2)*t*tbar*u*pow(m_B,-1)*pow(sigmabar,-2.0) +
                                       2.0*omega*u*(m_B*(3.0 - 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) +
                                       8.0*m_v*t*tbar*u*pow(sigmabar,-1) +
                                       (m_v2 - q2)*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u)*pow(m_B,-1)*pow(sigmabar,-2.0))))/8.0;

                    return C_1 * phi_3;
                }

                inline
                double I1_V1_phi_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_3 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*pow(sigmabar,-3.0)*(4.0*omega*u*
                                       (-(4.0*m_B2*sigmabar2*t*tbar*u) + q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_v2*(1 + (-2.0 - 4.0*t2 + 4.0*t)*u)) -
                                       4.0*m_B*u2*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(omega,2.0) -
                                       m_B*sigmabar*(-(4.0*(-m_v2 + m_B2*sigmabar2 - 2.0*m_B*m_v*sigmabar)*t*tbar*u) +
                                       q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0)))))/8.0;

                    return C_1 * phi_bar_3;
                }

                inline
                double I2_V1_phi_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_3 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

                    const double C_2 = (-(8.0*m_B2*m_v*t*tbar*u) + m_B*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(sigmabar,-1) +
                                       8.0*m_v3*t*tbar*u*pow(sigmabar,-2.0) +
                                       (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_B,-1)*pow(sigmabar,-3.0) +
                                       4.0*u2*(-(8.0*m_B*m_v*sigmabar*t*tbar*u) +
                                       m_B2*sigmabar2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       (m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))*pow(m_B,-1)*pow(omega,2.0)*
                                       pow(sigmabar,-3.0) + 2.0*omega*u*pow(m_B,-2.0)*
                                       (-(m_B4*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) + 16.0*m_B3*m_v*t*tbar*u*pow(sigmabar,-1) +
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)*pow(sigmabar,-4.0)))/8.0;

                    return C_2 * phi_bar_3;
                }

                inline
                double I1_V1_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to phi_double_bar_3 with a single pole in k2
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_1 = t*tbar*u*pow(sigmabar,-1);

                    return C_1 * phi_double_bar_3;
                }

                inline
                double I2_V1_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to phi_double_bar_3 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_2 = t*tbar*u*pow(m_B,-1)*(-(m_B3*sigmabar3) + 2.0*m_v*omega*u*(m_v - 2.0*omega*u) +
                                       m_B2*sigmabar2*(m_v + 4.0*omega*u) - 4.0*m_B*u2*sigmabar*pow(omega,2.0))*
                                       pow(sigmabar,-3.0);

                    return C_2 * phi_double_bar_3;
                }

                inline
                double I3_V1_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_3 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);

                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_3 = -(m_v*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*t*tbar*u*
                                       (-(m_B*sigmabar) + 2.0*omega*u)*(-m_v - m_B*sigmabar + 2.0*omega*u)*pow(m_B,-1)*
                                       pow(sigmabar,-4.0));

                    return C_3 * phi_double_bar_3;
                }

                inline
                double I3d1_V1_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_3 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3);
                    const double u2        = pow(u, 2);

                    const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);
                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_3 = -(m_v*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*t*tbar*u*
                                       (-(m_B*sigmabar) + 2.0*omega*u)*(-m_v - m_B*sigmabar + 2.0*omega*u)*pow(sigmabar,-4.0));

                    const double C_3d1 = m_v*t*tbar*u*pow(m_B,-1)*(8.0*omega*(m_v2 - q2)*u*(m_v - 2.0*omega*u) +
                                         m_B3*sigmabar3*(m_v + 4.0*omega*u) +
                                         2.0*m_B2*sigmabar2*(m_v2 + q2 - 6.0*m_v*omega*u - 4.0*u2*pow(omega,2.0)) -
                                         3.0*m_B*sigmabar*(m_v3 + 4.0*omega*q2*u - m_v*(q2 + 8.0*u2*pow(omega,2.0))))*
                                         pow(sigmabar,-5.0);

                    return C_3 * phi_bar_3 + C_3d1 * phi_double_bar_3;
                }

                inline
                double I1_V1_phi_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_4 = this->phi_4(omega_1, omega_2);

                    const double C_1 = (-((q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u - 4.0*(-m_v2 + m_B2*sigmabar2)*t*tbar*u)*
                                       pow(sigmabar,-1)) - 4.0*u2*(1 - 2.0*u)*pow(omega,2.0)*pow(sigmabar,-1) +
                                       2.0*omega*u*pow(m_B,-1)*(m_B2*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) +
                                       (m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(sigmabar,-2.0)))/8.0;

                    return C_1 * phi_4;
                }

                inline
                double I1_V1_phi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*(-(m_B*sigmabar*(4.0*m_B*sigmabar*(-m_v - m_B*sigmabar)*t*tbar*u +
                                       q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))) +
                                       4.0*omega*u*(-((m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) -
                                       2.0*m_B*m_v*sigmabar*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) -
                                       4.0*m_B*u2*sigmabar*(-1 + 2.0*u)*pow(omega,2.0))*pow(sigmabar,-3.0))/8.0;

                    return C_1 * phi_bar_4;
                }

                inline
                double I2_V1_phi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

                    const double C_2 = (pow(m_B,-2.0)*pow(sigmabar,-4.0)*(4.0*m_B*u2*sigmabar*
                                       (m_B2*sigmabar2*(-1 + 2.0*u) + (m_v2 - q2)*(-1 + 2.0*u) +
                                       4.0*m_B*m_v*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))*pow(omega,2.0) +
                                       2.0*omega*u*(8.0*m_B2*sigmabar2*q2*t*tbar*u -
                                       4.0*m_B3*m_v*sigmabar3*(-1 + (2.0 + 3.0*t2 - 3.0*t)*u) +
                                       m_B4*sigmabar4*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) +
                                       4.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)) -
                                       m_B*sigmabar*(-(4.0*m_B3*m_v*sigmabar3*t*tbar*u) - 4.0*m_B4*sigmabar4*t*tbar*u -
                                       (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       4.0*m_B*m_v*sigmabar*(-(m_v2*t*tbar*u) + q2*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) -
                                       m_B2*sigmabar2*(-(4.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/8.0;

                    return C_2 * phi_bar_4;
                }

                inline
                double I1_V1_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to phi_double_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double t2        = pow(t, 2);

                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_1 = -(u*(2.0*m_v*omega*(1 + (-2.0 - 4.0*t2 + 4.0*t)*u) -
                                       m_B*sigmabar*(-(m_v*t*tbar) + 2.0*omega*(1 + (-2.0 - 6.0*t2 + 6.0*t)*u)))*
                                       pow(m_B,-2.0)*pow(sigmabar,-3.0))/2.0;

                    return C_1 * phi_double_bar_4;
                }

                inline
                double I2_V1_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_2 = (pow(m_B,-2.0)*(m_B*sigmabar*(-(m_v3*t*tbar*u) - m_B2*m_v*sigmabar2*t*tbar*u +
                                       m_v*q2*(-1 + (2.0 + 3.0*t2 - 3.0*t)*u) -
                                       m_B*sigmabar*(-(2.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))) -
                                       2.0*omega*u*(2.0*m_B2*m_v*sigmabar2*t*tbar*u -
                                       2.0*m_v*(m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B3*sigmabar3*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) -
                                       m_B*sigmabar*(m_v2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))) +
                                       4.0*m_B*u2*sigmabar*(m_v*(-1 + 2.0*(1 + t2 - t)*u) -
                                       m_B*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))*pow(omega,2.0))*pow(sigmabar,-4.0))/2.0;

                    return C_2 * phi_double_bar_4;
                }

                inline
                double I3_V1_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(m_v*pow(m_B,-2.0)*(-(m_B*sigmabar*(-(2.0*m_B3*m_v*sigmabar3*t*tbar*u) +
                                       2.0*m_B4*sigmabar4*t*tbar*u + (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       2.0*m_B*m_v*sigmabar*(m_v2*t*tbar*u + q2*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) +
                                       m_B2*sigmabar2*(-(2.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)))) -
                                       4.0*m_B*u2*sigmabar*(m_B2*sigmabar2*(-1 + 2.0*(1 + t2 - t)*u) +
                                       (m_v2 - q2)*(-1 + 2.0*(1 + t2 - t)*u) +
                                       2.0*m_B*m_v*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))*pow(omega,2.0) +
                                       2.0*omega*u*(-(4.0*m_B2*sigmabar2*q2*t*tbar*u) + m_B4*sigmabar4*(-1 + 2.0*u) -
                                       2.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                       2.0*m_B3*m_v*sigmabar3*(-1 + (2.0 + 7.0*t2 - 7.0*t)*u) -
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)))*pow(sigmabar,-5.0))/2.0;

                    return C_3 * phi_double_bar_4;
                }

                inline
                double I3d1_V1_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);
                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(m_v*(-(sigmabar*(-(2.0*m_B3*m_v*sigmabar3*t*tbar*u) + 2.0*m_B4*sigmabar4*t*tbar*u +
                                       (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       2.0*m_B*m_v*sigmabar*(m_v2*t*tbar*u + q2*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) +
                                       m_B2*sigmabar2*(-(2.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)))) -
                                       4.0*u2*sigmabar*(m_B2*sigmabar2*(-1 + 2.0*(1 + t2 - t)*u) +
                                       (m_v2 - q2)*(-1 + 2.0*(1 + t2 - t)*u) +
                                       2.0*m_B*m_v*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))*pow(omega,2.0) +
                                       2.0*omega*u*pow(m_B,-1)*(-(4.0*m_B2*sigmabar2*q2*t*tbar*u) +
                                       m_B4*sigmabar4*(-1 + 2.0*u) -
                                       2.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                       2.0*m_B3*m_v*sigmabar3*(-1 + (2.0 + 7.0*t2 - 7.0*t)*u) -
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)))*pow(sigmabar,-5.0))/2.0;

                    const double C_3d1 = m_v*pow(m_B,-2.0)*(-(2.0*m_B*omega*sigmabar*u*
                                         (4.0*m_B*q2*sigmabar*t*tbar*u - 2.0*m_B3*sigmabar3*(-1 + 2.0*u) +
                                         m_v*(m_v2 - q2)*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) -
                                         3.0*m_B2*m_v*sigmabar2*(-1 + (2.0 + 7.0*t2 - 7.0*t)*u))) +
                                         2.0*m_B*sigmabar*(-(2.0*m_B3*m_v*sigmabar3*t*tbar*u) + 2.0*m_B4*sigmabar4*t*tbar*u +
                                         (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                         2.0*m_B*m_v*sigmabar*(m_v2*t*tbar*u + q2*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) +
                                         m_B2*sigmabar2*(-(2.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))) +
                                         m_B2*sigmabar2*(3.0*m_B2*m_v*sigmabar2*t*tbar*u - 4.0*m_B3*sigmabar3*t*tbar*u +
                                         m_v*(q2 + q2*(-2.0 - 5.0*t2 + 5.0*t)*u - m_v2*t*tbar*u) -
                                         m_B*sigmabar*(-(2.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))) +
                                         4.0*m_B2*sigmabar2*u2*(m_v + m_v*(-2.0 - 6.0*t2 + 6.0*t)*u -
                                         m_B*sigmabar*(-1 + 2.0*(1 + t2 - t)*u))*pow(omega,2.0) +
                                         8.0*m_B*u2*sigmabar*(m_B2*sigmabar2*(-1 + 2.0*(1 + t2 - t)*u) +
                                         (m_v2 - q2)*(-1 + 2.0*(1 + t2 - t)*u) +
                                         2.0*m_B*m_v*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))*pow(omega,2.0) -
                                         5.0*omega*u*(-(4.0*m_B2*sigmabar2*q2*t*tbar*u) + m_B4*sigmabar4*(-1 + 2.0*u) -
                                         2.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                         2.0*m_B3*m_v*sigmabar3*(-1 + (2.0 + 7.0*t2 - 7.0*t)*u) -
                                         (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)))*pow(sigmabar,-6.0);

                    return C_3 * phi_bar_4 + C_3d1 * phi_double_bar_4;
                }

                inline
                double I1_V1_psi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*(-(4.0*omega*u*(-(2.0*m_B2*sigmabar2*t*tbar*u) -
                                       (m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B*m_v*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))) +
                                       m_B*sigmabar*(-(2.0*(-m_v2 + m_B2*sigmabar2)*t*tbar*u) +
                                       q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) +
                                       4.0*m_B*u2*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(omega,2.0))*pow(sigmabar,-3.0))
                                       /4.0;

                    return C_1 * psi_bar_4;
                }

                inline
                double I2_V1_psi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

                    const double C_2 = ((m_v2 + m_B2*sigmabar2 - q2 + 2.0*m_B*m_v*sigmabar)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*
                                       pow(m_B,-2.0)*(2.0*m_B2*omega*sigmabar2*u + 2.0*omega*(-m_v2 + q2)*u -
                                       m_B*sigmabar*(q2 + 4.0*u2*pow(omega,2.0)))*pow(sigmabar,-4.0))/4.0;

                    return C_2 * psi_bar_4;
                }

                inline
                double I1_V1_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to psi_double_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double t2        = pow(t, 2);

                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_1 = -(u*(-(m_B*sigmabar*(m_v - 2.0*m_B*sigmabar)*t*tbar) +
                                       2.0*omega*(m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)))*pow(m_B,-2.0)*pow(sigmabar,-3.0))/2.0;

                    return C_1 * psi_double_bar_4;
                }

                inline
                double I2_V1_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_double_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2),m_B3 = pow(m_B, 3);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_2 = (pow(m_B,-2.0)*(m_B*sigmabar*(-(m_B2*m_v*sigmabar2*t*tbar*u) + 2.0*m_B3*sigmabar3*t*tbar*u +
                                       m_v*(q2 + q2*(-2.0 - 3.0*t2 + 3.0*t)*u + m_v2*t*tbar*u) +
                                       m_B*sigmabar*(-(2.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))) +
                                       2.0*omega*u*(2.0*m_B2*m_v*sigmabar2*t*tbar*u -
                                       m_B3*sigmabar3*(-1 + 2.0*(1 + t2 - t)*u) -
                                       2.0*m_v*(m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B*sigmabar*(m_v2*(-1 + 2.0*(1 + t2 - t)*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))) -
                                       4.0*m_B*u2*sigmabar*(m_v - m_B*sigmabar)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*
                                       pow(omega,2.0))*pow(sigmabar,-4.0))/2.0;

                    return C_2 * psi_double_bar_4;
                }

                inline
                double I3_V1_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(m_v*(m_v2 + m_B2*sigmabar2 - q2 + 2.0*m_B*m_v*sigmabar)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*
                                       pow(m_B,-2.0)*(-(2.0*m_B2*omega*sigmabar2*u) + 2.0*omega*(m_v2 - q2)*u +
                                       m_B*sigmabar*(q2 + 4.0*u2*pow(omega,2.0)))*pow(sigmabar,-5.0))/2.0;

                    return C_3 * psi_double_bar_4;
                }

                inline
                double I3d1_V1_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);
                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(m_v*(m_v2 + m_B2*sigmabar2 - q2 + 2.0*m_B*m_v*sigmabar)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*
                                       pow(m_B,-1)*(-(2.0*m_B2*omega*sigmabar2*u) + 2.0*omega*(m_v2 - q2)*u +
                                       m_B*sigmabar*(q2 + 4.0*u2*pow(omega,2.0)))*pow(sigmabar,-5.0))/2.0;

                    const double C_3d1 = m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_B,-2.0)*
                                         (m_B4*omega*sigmabar4*u - m_B3*sigmabar3*(q2 + 4.0*omega*u*(-m_v + omega*u)) -
                                         2.0*m_B*(m_v2 - q2)*sigmabar*(q2 + 4.0*omega*u*(m_v + omega*u)) -
                                         3.0*m_B2*m_v*sigmabar2*(q2 + 4.0*u2*pow(omega,2.0)) - 5.0*omega*u*pow(m_v2 - q2,2.0))*
                                         pow(sigmabar,-6.0);

                    return C_3 * psi_bar_4 + C_3d1 * psi_double_bar_4;
                }

                inline
                double I1_V1_chi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);

                    const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

                    const double C_1 = (u*(-((-m_v2 + m_B2*sigmabar2 + q2 - 2.0*m_B*m_v*sigmabar)*t*tbar) +
                                       omega*(m_v + m_v*(-2.0 - 4.0*t2 + 4.0*t)*u + 4.0*m_B*sigmabar*t*tbar*u))*pow(m_B,-1)*
                                       pow(sigmabar,-2.0))/2.0;

                    return C_1 * chi_bar_4;
                }

                inline
                double I2_V1_chi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

                    const double C_2 = -(m_v*(-(sigmabar*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u - 4.0*(-m_v2 + m_B2*sigmabar2)*t*tbar*u)) +
                                       2.0*omega*u*(m_B2*sigmabar2*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) +
                                       (m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))*pow(m_B,-1) -
                                       4.0*u2*sigmabar*(1 - 2.0*u)*pow(omega,2.0))*pow(sigmabar,-3.0))/4.0;

                    return C_2 * chi_bar_4;
                }

                inline
                double I1_V1_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to chi_double_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);

                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_1 = -(u*(-(2.0*m_B2*sigmabar2*t*tbar) + 2.0*m_v*omega*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_B*sigmabar*(-(m_v*t*tbar) + 2.0*omega*(1 + (-2.0 - 6.0*t2 + 6.0*t)*u)))*
                                       pow(m_B,-2.0)*pow(sigmabar,-3.0))/2.0;

                    return C_1 * chi_double_bar_4;
                }

                inline
                double I2_V1_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_double_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_2 = (pow(m_B,-2.0)*pow(sigmabar,-4.0)*(2.0*omega*u*
                                       (2.0*m_B2*m_v*sigmabar2*t*tbar*u -
                                       2.0*m_v*(m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B3*sigmabar3*(-1 + 2.0*(1 + 5.0*t2 - 5.0*t)*u) -
                                       m_B*(m_v2 + q2)*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) -
                                       m_B*sigmabar*(-(m_v3*t*tbar*u) - 3.0*m_B2*m_v*sigmabar2*t*tbar*u +
                                       2.0*m_B3*sigmabar3*t*tbar*u + m_v*q2*(-1 + (2.0 + 3.0*t2 - 3.0*t)*u) -
                                       m_B*sigmabar*(-(2.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))) -
                                       4.0*m_B*u2*sigmabar*pow(omega,2.0)*
                                       (m_v*(-1 + 2.0*u) - m_B*sigmabar*(-1 + 2.0*u*pow(1 - 2.0*t,2.0)))))/2.0;

                    return C_2 * chi_double_bar_4;
                }

                inline
                double I3_V1_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(m_v*pow(m_B,-2.0)*pow(sigmabar,-5.0)*
                                       (2.0*omega*u*(8.0*m_B2*sigmabar2*q2*t*tbar*u +
                                       m_B4*sigmabar4*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) -
                                       2.0*m_B3*m_v*sigmabar3*(-1 + 2.0*(1 + 5.0*t2 - 5.0*t)*u) +
                                       2.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)) +
                                       4.0*m_B*u2*sigmabar*pow(omega,2.0)*
                                       (m_B2*sigmabar2*(-1 + 2.0*u) + (m_v2 - q2)*(-1 + 2.0*u) +
                                       2.0*m_B*m_v*sigmabar*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))) -
                                       m_B*sigmabar*(4.0*m_B3*m_v*sigmabar3*t*tbar*u - 4.0*m_B4*sigmabar4*t*tbar*u -
                                       (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       2.0*m_B*m_v*sigmabar*(2.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) -
                                       m_B2*sigmabar2*(-(4.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/2.0;

                    return C_3 * chi_double_bar_4;
                }

                inline
                double I3d1_V1_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);
                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(m_v*pow(sigmabar,-5.0)*(2.0*omega*u*pow(m_B,-1)*
                                       (8.0*m_B2*sigmabar2*q2*t*tbar*u + m_B4*sigmabar4*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) -
                                       2.0*m_B3*m_v*sigmabar3*(-1 + 2.0*(1 + 5.0*t2 - 5.0*t)*u) +
                                       2.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)) +
                                       4.0*u2*sigmabar*pow(omega,2.0)*
                                       (m_B2*sigmabar2*(-1 + 2.0*u) + (m_v2 - q2)*(-1 + 2.0*u) +
                                       2.0*m_B*m_v*sigmabar*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))) -
                                       sigmabar*(4.0*m_B3*m_v*sigmabar3*t*tbar*u - 4.0*m_B4*sigmabar4*t*tbar*u -
                                       (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       2.0*m_B*m_v*sigmabar*(2.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) -
                                       m_B2*sigmabar2*(-(4.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/2.0;

                    const double C_3d1 = m_v*pow(m_B,-2.0)*pow(sigmabar,-6.0)*(-(2.0*m_B*omega*sigmabar*u*
                                         (-(8.0*m_B*q2*sigmabar*t*tbar*u) - 2.0*m_B3*sigmabar3*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) +
                                         3.0*m_B2*m_v*sigmabar2*(-1 + 2.0*(1 + 5.0*t2 - 5.0*t)*u) -
                                         m_v*(m_v2 - q2)*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))) -
                                         5.0*omega*u*(8.0*m_B2*sigmabar2*q2*t*tbar*u +
                                         m_B4*sigmabar4*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) -
                                         2.0*m_B3*m_v*sigmabar3*(-1 + 2.0*(1 + 5.0*t2 - 5.0*t)*u) +
                                         2.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                         (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)) -
                                         4.0*m_B2*sigmabar2*u2*pow(omega,2.0)*
                                         (-(m_B*sigmabar*(-1 + 2.0*u)) + m_v*(1 - 2.0*u*pow(1 - 2.0*t,2.0))) -
                                         8.0*m_B*u2*sigmabar*pow(omega,2.0)*
                                         (m_B2*sigmabar2*(-1 + 2.0*u) + (m_v2 - q2)*(-1 + 2.0*u) +
                                         2.0*m_B*m_v*sigmabar*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))) +
                                         2.0*m_B*sigmabar*(4.0*m_B3*m_v*sigmabar3*t*tbar*u - 4.0*m_B4*sigmabar4*t*tbar*u -
                                         (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                         2.0*m_B*m_v*sigmabar*(2.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) -
                                         m_B2*sigmabar2*(-(4.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0)))) +
                                         m_B2*sigmabar2*(2.0*m_v3*t*tbar*u - 6.0*m_B2*m_v*sigmabar2*t*tbar*u +
                                         8.0*m_B3*sigmabar3*t*tbar*u + m_v*q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                         m_B*sigmabar*(-(4.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0)))));

                    return C_3 * chi_bar_4 + C_3d1 * chi_double_bar_4;
                }

                double integrand_V1(const std::array<double, 4> & args, const double & q2) const
                {
                    const double sigma = args[0];
                    const double x     = args[1];
                    const double t     = args[2];
                    const double u     = args[3];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = - f_B() * m_B / (8.0  * M_PI * M_PI * f_V() * m_V());
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_V2 = pow(m_V(), 2);
                    const double M4   = pow(M2, 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_V2) / M2());

                    const double I1 = I1_V1_phi_3(sigma, omega, t, u, q2)             + I1_V1_phi_4(sigma, omega, t, u, q2)
                                    + I1_V1_phi_bar_3(sigma, omega, t, u, q2)         + I1_V1_phi_bar_4(sigma, omega, t, u, q2)
                                    + I1_V1_phi_double_bar_3(sigma, omega, t, u, q2)  + I1_V1_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I1_V1_psi_bar_4(sigma, omega, t, u, q2)         + I1_V1_chi_bar_4(sigma, omega, t, u, q2)
                                    + I1_V1_psi_double_bar_4(sigma, omega, t, u, q2)  + I1_V1_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I2 = I2_V1_phi_bar_3(sigma, omega, t, u, q2)         + I2_V1_phi_bar_4(sigma, omega, t, u, q2)
                                    + I2_V1_phi_double_bar_3(sigma, omega, t, u, q2)  + I2_V1_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I2_V1_psi_bar_4(sigma, omega, t, u, q2)         + I2_V1_chi_bar_4(sigma, omega, t, u, q2)
                                    + I2_V1_psi_double_bar_4(sigma, omega, t, u, q2)  + I2_V1_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3 = I3_V1_phi_double_bar_3(sigma, omega, t, u, q2)  + I3_V1_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I3_V1_psi_double_bar_4(sigma, omega, t, u, q2)  + I3_V1_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result = 0.0;
                    result += - I1;
                    result +=   I2 / M2;
                    result += - I3 / (2.0 * M4);
                    result *=   prefactor * exp;

                    return result;
                }

                double surface_V1(const std::array<double, 3> & args, const double & sigma, const double & q2) const
                {
                    const double x     = args[0];
                    const double t     = args[1];
                    const double u     = args[2];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = - f_B() * m_B / (8.0  * M_PI * M_PI * f_V() * m_V());
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_B2 = pow(m_B(), 2);
                    const double m_v2 = pow(m_v(), 2);
                    const double m_V2 = pow(m_V(), 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_V2) / M2());

                    const double sigmabar = 1.0 - sigma, sigmabar2 = pow(sigmabar, 2);
                    const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
                    const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;

                    const double I2   = I2_V1_phi_bar_3(sigma, omega, t, u, q2)           + I2_V1_phi_bar_4(sigma, omega, t, u, q2)
                                      + I2_V1_phi_double_bar_3(sigma, omega, t, u, q2)    + I2_V1_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I2_V1_psi_bar_4(sigma, omega, t, u, q2)           + I2_V1_chi_bar_4(sigma, omega, t, u, q2)
                                      + I2_V1_psi_double_bar_4(sigma, omega, t, u, q2)    + I2_V1_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3   = I3_V1_phi_double_bar_3(sigma, omega, t, u, q2)    + I3_V1_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3_V1_psi_double_bar_4(sigma, omega, t, u, q2)    + I3_V1_chi_double_bar_4(sigma, omega, t, u, q2);

                    const double I3d1 = I3d1_V1_phi_double_bar_3(sigma, omega, t, u, q2)  + I3d1_V1_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3d1_V1_psi_double_bar_4(sigma, omega, t, u, q2)  + I3d1_V1_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result = 0.0;
                    result +=  eta * I2 / m_B2;
                    result += -0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
                    result *=  prefactor * exp;

                    return result;
                }

                double integrand_V1_m1(const std::array<double, 4> & args, const double & q2) const
                {
                    const double sigma = args[0];
                    const double x     = args[1];
                    const double t     = args[2];
                    const double u     = args[3];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = - f_B() * m_B / (8.0  * M_PI * M_PI * f_V() * m_V());
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_V2 = pow(m_V(), 2);
                    const double M4   = pow(M2, 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_V2) / M2());

                    const double I1 = I1_V1_phi_3(sigma, omega, t, u, q2)             + I1_V1_phi_4(sigma, omega, t, u, q2)
                                    + I1_V1_phi_bar_3(sigma, omega, t, u, q2)         + I1_V1_phi_bar_4(sigma, omega, t, u, q2)
                                    + I1_V1_phi_double_bar_3(sigma, omega, t, u, q2)  + I1_V1_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I1_V1_psi_bar_4(sigma, omega, t, u, q2)         + I1_V1_chi_bar_4(sigma, omega, t, u, q2)
                                    + I1_V1_psi_double_bar_4(sigma, omega, t, u, q2)  + I1_V1_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I2 = I2_V1_phi_bar_3(sigma, omega, t, u, q2)         + I2_V1_phi_bar_4(sigma, omega, t, u, q2)
                                    + I2_V1_phi_double_bar_3(sigma, omega, t, u, q2)  + I2_V1_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I2_V1_psi_bar_4(sigma, omega, t, u, q2)         + I2_V1_chi_bar_4(sigma, omega, t, u, q2)
                                    + I2_V1_psi_double_bar_4(sigma, omega, t, u, q2)  + I2_V1_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3 = I3_V1_phi_double_bar_3(sigma, omega, t, u, q2)  + I3_V1_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I3_V1_psi_double_bar_4(sigma, omega, t, u, q2)  + I3_V1_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result1 = 0.0;
                    result1 += - I1;
                    result1 +=   I2 / M2;
                    result1 += - I3 / (2.0 * M4);
                    result1 *=   exp * s(sigma, q2);

                    double result2 = 0.0;
                    result2 += - I2;
                    result2 +=   I3 / M2;
                    result2 *=   exp;

                    return prefactor * (result1 + result2);
                }

                double surface_V1_m1(const std::array<double, 3> & args, const double & sigma, const double & q2) const
                {
                    const double x     = args[0];
                    const double t     = args[1];
                    const double u     = args[2];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = - f_B() * m_B / (8.0  * M_PI * M_PI * f_V() * m_V());
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_B2 = pow(m_B(), 2);
                    const double m_v2 = pow(m_v(), 2);
                    const double m_V2 = pow(m_V(), 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_V2) / M2());

                    const double sigmabar = 1.0 - sigma, sigmabar2 = pow(sigmabar, 2);
                    const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
                    const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;

                    const double I2   = I2_V1_phi_bar_3(sigma, omega, t, u, q2)           + I2_V1_phi_bar_4(sigma, omega, t, u, q2)
                                      + I2_V1_phi_double_bar_3(sigma, omega, t, u, q2)    + I2_V1_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I2_V1_psi_bar_4(sigma, omega, t, u, q2)           + I2_V1_chi_bar_4(sigma, omega, t, u, q2)
                                      + I2_V1_psi_double_bar_4(sigma, omega, t, u, q2)    + I2_V1_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3   = I3_V1_phi_double_bar_3(sigma, omega, t, u, q2)    + I3_V1_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3_V1_psi_double_bar_4(sigma, omega, t, u, q2)    + I3_V1_chi_double_bar_4(sigma, omega, t, u, q2);

                    const double I3d1 = I3d1_V1_phi_double_bar_3(sigma, omega, t, u, q2)  + I3d1_V1_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3d1_V1_psi_double_bar_4(sigma, omega, t, u, q2)  + I3d1_V1_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result1 =   0.0;
                           result1 +=  eta * I2 / m_B2;
                           result1 += -0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
                           result1 *=  exp * s(sigma, q2);

                    double result2 =  0.0;
                           result2 += 0.5 * eta * I3 / m_B2;
                           result2 *= exp;

                    return prefactor * (result1 + result2);
                }

                double V1(const double & q2) const
                {
                    const double sigma_0 = this->sigma_0(q2, s0_0_V1(), s0_1_V1());

                    double value_integral = 0.0;
                    double value_surface  = 0.0;

                    const std::function<double (const std::array<double, 4> &)> integrand = std::bind(&LCSR::integrand_V1, this, std::placeholders::_1, q2);
                    const std::function<double (const std::array<double, 3> &)> surface = std::bind(&LCSR::surface_V1, this, std::placeholders::_1, sigma_0, q2);

                    value_integral = integrate(integrand, { 0.0, 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0, 1.0 }, cubature::Config());
                    value_surface  = integrate(surface, { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, cubature::Config());

                    return value_integral + value_surface;
                }

                double normalized_first_moment_V1(const double & q2) const
                {
                    const double sigma_0 = this->sigma_0(q2, s0_0_V1(), s0_1_V1());

                    double value_integral = 0.0;
                    double value_surface  = 0.0;

                    const std::function<double (const std::array<double, 4> &)> integrand = std::bind(&LCSR::integrand_V1, this, std::placeholders::_1, q2);
                    const std::function<double (const std::array<double, 3> &)> surface = std::bind(&LCSR::surface_V1, this, std::placeholders::_1, sigma_0, q2);

                    value_integral = integrate(integrand, { 0.0, 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0, 1.0 }, cubature::Config());
                    value_surface  = integrate(surface, { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, cubature::Config());

                    const double denominator = value_integral + value_surface;

                    double value_integral_m1 = 0.0;
                    double value_surface_m1  = 0.0;

                    const std::function<double (const std::array<double, 4> &)> integrand_m1 = std::bind(&LCSR::integrand_V1_m1, this, std::placeholders::_1, q2);
                    const std::function<double (const std::array<double, 3> &)> surface_m1 = std::bind(&LCSR::surface_V1_m1, this, std::placeholders::_1, sigma_0, q2);

                    value_integral_m1 = integrate(integrand_m1, { 0.0, 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0, 1.0 }, cubature::Config());
                    value_surface_m1  = integrate(surface_m1, { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, cubature::Config());

                    const double numerator = value_integral_m1 + value_surface_m1;

                    return numerator / denominator;
                }

                inline
                double I1_V2_phi_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_3 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_3 = this->phi_3(omega_1, omega_2);

                    const double C_1 = -((4.0*m_B3*m_v*sigmabar4*t*tbar*u - m_B2*sigmabar3*q2*
                                       (-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) -
                                       (m_v2 - q2)*q2*sigmabar*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) +
                                       2.0*m_B*m_v*sigmabar2*(-(4.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       2.0*omega*(-m_v2 + m_B2*sigmabar2 + q2)*u*
                                       (-(2.0*m_B*m_v*sigmabar*(-1 + 2.0*u)) +
                                       m_B2*sigmabar2*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) +
                                       (m_v2 - q2)*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u))*pow(m_B,-1) -
                                       4.0*u2*(1 - sigma)*(-(2.0*m_B*m_v*sigmabar*(-1 + 2.0*u)) +
                                       m_B2*sigmabar2*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) +
                                       (m_v2 - q2)*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u))*pow(omega,2.0) +
                                       4.0*m_v*t*tbar*u*pow(m_B,-1)*pow(m_v2 - q2,2.0))*pow(sigmabar,-3.0))/16.0;

                    return C_1 * phi_3;
                }

                inline
                double I1_V2_phi_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_3 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

                    const double C_1 = (pow(m_B,-2.0)*pow(sigmabar,-4.0)*(-(4.0*m_B*u2*sigmabar*
                                       (-(4.0*m_B2*sigmabar2*t*tbar*u) + m_B*m_v*(-sigmabar + 2.0*u - 2.0*sigma*u) -
                                       (m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))*pow(omega,2.0)) +
                                       omega*u*(-(16.0*m_B3*m_v*sigmabar3*t*tbar*u) -
                                       4.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_B4*sigmabar4*(-1 + 2.0*(1 + 6.0*t2 - 6.0*t)*u) +
                                       8.0*m_B2*sigmabar2*(-(q2*t*tbar*u) + m_v2*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) +
                                       3.0*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)) -
                                       m_B*sigmabar*(4.0*m_B2*m_v2*sigmabar2*t*tbar*u - 8.0*m_B3*m_v*sigmabar3*t*tbar*u -
                                       2.0*m_B4*sigmabar4*t*tbar*u + (m_v2 - q2)*
                                       (q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u - 2.0*m_v2*t*tbar*u) +
                                       m_B*m_v*sigmabar*(8.0*m_v2*t*tbar*u + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/8.0;

                    return C_1 * phi_bar_3;
                }

                inline
                double I2_V2_phi_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_3 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4), m_B5 = pow(m_B, 5), m_B6 = pow(m_B, 6);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4), sigmabar5 = pow(sigmabar, 5), sigmabar6 = pow(sigmabar, 6);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

                    const double C_2 = -(pow(m_B,-2.0)*pow(sigmabar,-5.0)*(-(4.0*m_B*u2*sigmabar*pow(omega,2.0)*
                                       (-(2.0*m_B3*m_v*sigmabar3*(-1 + 2.0*u)) -
                                       2.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + 2.0*u) +
                                       m_B4*sigmabar4*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       2.0*m_B2*sigmabar2*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u +
                                       m_v2*(-5.0 + 2.0*(5.0 + 14.0*t2 - 14.0*t)*u)) +
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0))) +
                                       2.0*omega*u*(16.0*m_B3*m_v*sigmabar3*(m_v2 + q2)*t*tbar*u +
                                       2.0*m_B5*m_v*sigmabar5*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_B6*sigmabar6*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_B4*sigmabar4*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u +
                                       m_v2*(-9.0 + 2.0*(9.0 + 26*t2 - 26*t)*u)) -
                                       m_B2*sigmabar2*(m_v2 - q2)*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u +
                                       m_v2*(-9.0 + 2.0*(9.0 + 26*t2 - 26*t)*u)) +
                                       2.0*m_B*m_v*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0) -
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,3.0)) -
                                       m_B*sigmabar*(-(8.0*m_B5*m_v*sigmabar5*t*tbar*u) +
                                       m_B4*sigmabar4*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       2.0*m_B*m_v*(m_v2 - q2)*sigmabar*
                                       (q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u - 4.0*m_v2*t*tbar*u) +
                                       2.0*m_B2*sigmabar2*q2*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u +
                                       m_v2*(-5.0 + 2.0*(5.0 + 14.0*t2 - 14.0*t)*u)) +
                                       q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0) -
                                       2.0*m_B3*m_v*sigmabar3*(-(8.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/
                                       16.0;

                    return C_2 * phi_bar_3;
                }

                inline
                double I1_V2_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to phi_double_bar_3 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_1 = (m_v*u*(-(m_B2*sigmabar2*t*tbar) - m_B*sigmabar*
                                       (-(m_v*t*tbar) + omega*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       2.0*omega*(-(2.0*omega*u2*t*tbar) + m_v*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)))*pow(m_B,-1)*
                                       pow(sigmabar,-3.0))/2.0;

                    return C_1 * phi_double_bar_3;
                }

                inline
                double I2_V2_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_3 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_2 = -(m_v*pow(sigmabar,-4.0)*(-(sigmabar*(-(4.0*m_B2*m_v*sigmabar2*t*tbar*u) +
                                       4.0*m_B3*sigmabar3*t*tbar*u -
                                       2.0*m_v*(q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u - 2.0*m_v2*t*tbar*u) -
                                       m_B*sigmabar*(4.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)))) -
                                       4.0*u2*(-(4.0*(m_v2 - q2)*t*tbar*u) - m_B2*sigmabar2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       2.0*m_B*m_v*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))*pow(m_B,-1)*pow(omega,2.0) +
                                       2.0*omega*u*pow(m_B,-1)*(-(4.0*m_B2*m_v*sigmabar2*t*tbar*u) -
                                       4.0*m_v*(m_v2 - q2)*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) -
                                       m_B3*sigmabar3*(-1 + 2.0*u*pow(1 - 2.0*t,2.0)) +
                                       m_B*sigmabar*(m_v2*(-5.0 + 2.0*(5.0 + 12.0*t2 - 12.0*t)*u) +
                                       q2*(1 - 2.0*u*pow(1 - 2.0*t,2.0))))))/4.0;

                    return C_2 * phi_double_bar_3;
                }

                inline
                double I3_V2_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_3 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_3 = (m_v*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                       (-(sigmabar*(m_B*m_v2*t*(-sigmabar + t - sigma*t)*u - m_v3*t*tbar*u +
                                       m_B*(m_B2*sigmabar2 - q2)*sigmabar*t*tbar*u +
                                       m_v*(q2 + q2*(-2.0 - 5.0*t2 + 5.0*t)*u + m_B2*sigmabar2*t*tbar*u))) +
                                       2.0*omega*u*(2.0*m_B3*sigmabar3*t*tbar*u - 2.0*m_B*q2*sigmabar*t*tbar*u +
                                       m_v*(m_v2 - q2)*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) -
                                       m_B2*m_v*sigmabar2*(-1 + (2.0 + 7.0*t2 - 7.0*t)*u))*pow(m_B,-1) +
                                       4.0*u2*(-(m_B2*sigmabar2*t*tbar*u) - (m_v2 - q2)*t*tbar*u +
                                       m_B*m_v*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))*pow(m_B,-1)*pow(omega,2.0))*
                                       pow(sigmabar,-5.0))/2.0;

                    return C_3 * phi_double_bar_3;
                }

                inline
                double I3d1_V2_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_3 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);
                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_3 = -(m_v*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                       (-(2.0*omega*u*(2.0*m_B3*sigmabar3*t*tbar*u - 2.0*m_B*q2*sigmabar*t*tbar*u +
                                       m_v*(m_v2 - q2)*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) -
                                       m_B2*m_v*sigmabar2*(-1 + (2.0 + 7.0*t2 - 7.0*t)*u))) +
                                       m_B*sigmabar*(m_B*m_v2*t*(-sigmabar + t - sigma*t)*u - m_v3*t*tbar*u +
                                       m_B*(m_B2*sigmabar2 - q2)*sigmabar*t*tbar*u +
                                       m_v*(q2 + q2*(-2.0 - 5.0*t2 + 5.0*t)*u + m_B2*sigmabar2*t*tbar*u)) -
                                       4.0*u2*(-(m_B2*sigmabar2*t*tbar*u) - (m_v2 - q2)*t*tbar*u +
                                       m_B*m_v*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))*pow(omega,2.0))*pow(sigmabar,-5.0))/
                                       2.0;

                    const double C_3d1 = (m_v*(-(m_B*sigmabar2*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                         (m_v2 - 3.0*m_B2*sigmabar2 + q2 - 2.0*m_B*m_v*sigmabar)*t*tbar*u) +
                                         4.0*omega*(1 - sigma)*(m_v - m_B*sigmabar)*u*
                                         (2.0*m_B3*sigmabar3*t*tbar*u - 2.0*m_B*q2*sigmabar*t*tbar*u +
                                         m_v*(m_v2 - q2)*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) -
                                         m_B2*m_v*sigmabar2*(-1 + (2.0 + 7.0*t2 - 7.0*t)*u)) -
                                         2.0*m_B*sigmabar2*(m_v - m_B*sigmabar)*
                                         (m_B*m_v2*t*(-sigmabar + t - sigma*t)*u - m_v3*t*tbar*u +
                                         m_B*(m_B2*sigmabar2 - q2)*sigmabar*t*tbar*u +
                                         m_v*(q2 + q2*(-2.0 - 5.0*t2 + 5.0*t)*u + m_B2*sigmabar2*t*tbar*u)) -
                                         4.0*(1 - sigma)*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                         (m_B*m_v2*t*(-sigmabar + t - sigma*t)*u - m_v3*t*tbar*u +
                                         m_B*(m_B2*sigmabar2 - q2)*sigmabar*t*tbar*u +
                                         m_v*(q2 + q2*(-2.0 - 5.0*t2 + 5.0*t)*u + m_B2*sigmabar2*t*tbar*u)) +
                                         10.0*omega*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*u*
                                         (2.0*m_B3*sigmabar3*t*tbar*u - 2.0*m_B*q2*sigmabar*t*tbar*u +
                                         m_v*(m_v2 - q2)*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) -
                                         m_B2*m_v*sigmabar2*(-1 + (2.0 + 7.0*t2 - 7.0*t)*u))*pow(m_B,-1) +
                                         2.0*omega*(1 - sigma)*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*u*
                                         (-(6.0*m_B3*sigmabar2*t*tbar*u) + 2.0*m_B*q2*t*tbar*u +
                                         2.0*m_B2*m_v*sigmabar*(-1 + (2.0 + 7.0*t2 - 7.0*t)*u))*pow(m_B,-1) +
                                         4.0*u2*(1 - sigma)*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                         (m_v + m_v*(-2.0 - 6.0*t2 + 6.0*t)*u + 2.0*m_B*sigmabar*t*tbar*u)*pow(omega,2.0) +
                                         8.0*u2*(1 - sigma)*(m_v - m_B*sigmabar)*
                                         (-(m_B2*sigmabar2*t*tbar*u) - (m_v2 - q2)*t*tbar*u +
                                         m_B*m_v*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))*pow(omega,2.0) +
                                         20.0*u2*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                         (-(m_B2*sigmabar2*t*tbar*u) - (m_v2 - q2)*t*tbar*u +
                                         m_B*m_v*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))*pow(m_B,-1)*pow(omega,2.0))*
                                         pow(sigmabar,-6.0))/2.0;

                    return C_3 * phi_bar_3 + C_3d1 * phi_double_bar_3;
                }

                inline
                double I1_V2_phi_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_4 = this->phi_4(omega_1, omega_2);

                    const double C_1 = (pow(sigmabar,-3.0)*(4.0*u2*sigmabar*(m_B2*sigmabar2*(-1 + 2.0*u) + (m_v2 - q2)*(-1 + 2.0*u) -
                                       2.0*m_B*m_v*sigmabar*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u))*pow(omega,2.0) +
                                       2.0*omega*u*pow(m_B,-1)*(8.0*m_B2*sigmabar2*(m_v2 + q2)*t*tbar*u +
                                       m_B4*sigmabar4*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) +
                                       2.0*m_B3*m_v*sigmabar3*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) -
                                       2.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) +
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)) -
                                       sigmabar*(-(4.0*m_B4*sigmabar4*t*tbar*u) +
                                       2.0*m_B*m_v*q2*sigmabar*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) +
                                       (m_v2 - q2)*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u - 4.0*m_v2*t*tbar*u) -
                                       m_B2*sigmabar2*(-(8.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/16.0;

                    return C_1 * phi_4;
                }

                inline
                double I1_V2_phi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*pow(sigmabar,-4.0)*(-(m_B*sigmabar*
                                       (4.0*m_B2*m_v2*sigmabar2*t*tbar*u - 4.0*m_B3*m_v*sigmabar3*t*tbar*u -
                                       4.0*m_B4*sigmabar4*t*tbar*u + (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B*m_v*sigmabar*(-(4.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)))) +
                                       omega*u*(-(8.0*m_B3*m_v*sigmabar3*t*tbar*u) +
                                       m_B4*sigmabar4*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) +
                                       4.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       8.0*m_B2*sigmabar2*(-(q2*t*tbar*u) + m_v2*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) -
                                       3.0*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)) -
                                       4.0*m_B*u2*sigmabar*pow(omega,2.0)*
                                       ((m_v2 - q2)*(-1 + 2.0*u) - m_B*m_v*sigmabar*(-1 + 2.0*u*pow(1 - 2.0*t,2.0)))))/8.0;

                    return C_1 * phi_bar_4;
                }

                inline
                double I2_V2_phi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

                    const double C_2 = ((m_v2 + m_B2*sigmabar2 - q2 + 2.0*m_B*m_v*sigmabar)*pow(m_B,-2.0)*pow(sigmabar,-5.0)*
                                       (4.0*m_B*u2*sigmabar*(m_B2*sigmabar2*(-1 + 2.0*u) + (m_v2 - q2)*(-1 + 2.0*u) -
                                       4.0*m_B*m_v*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))*pow(omega,2.0) +
                                       2.0*omega*u*(8.0*m_B2*sigmabar2*q2*t*tbar*u +
                                       4.0*m_B3*m_v*sigmabar3*(-1 + (2.0 + 3.0*t2 - 3.0*t)*u) +
                                       m_B4*sigmabar4*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) -
                                       4.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)) -
                                       m_B*sigmabar*(4.0*m_B3*m_v*sigmabar3*t*tbar*u - 4.0*m_B4*sigmabar4*t*tbar*u -
                                       (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       4.0*m_B*m_v*sigmabar*(-(m_v2*t*tbar*u) + q2*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) -
                                       m_B2*sigmabar2*(-(4.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/16.0;

                    return C_2 * phi_bar_4;
                }

                inline
                double I1_V2_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double m_v2      = pow(m_v, 2);

                    const double C_1 = (m_v*pow(m_B,-2.0)*(-(m_B*sigmabar*(q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u -
                                       2.0*m_v*(m_v - m_B*sigmabar)*t*tbar*u)) +
                                       2.0*omega*u*(3.0*(m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       2.0*m_B*m_v*sigmabar*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) +
                                       4.0*m_B*u2*sigmabar*(-1 + 2.0*(1 + t2 - t)*u)*pow(omega,2.0))*pow(sigmabar,-4.0))/4.0;

                    return C_1 * phi_double_bar_4;
                }

                inline
                double I2_V2_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_2 = (m_v*pow(m_B,-2.0)*pow(sigmabar,-5.0)*(-(m_B*sigmabar*
                                       (-(4.0*m_B3*m_v*sigmabar3*t*tbar*u) + 3.0*m_B4*sigmabar4*t*tbar*u -
                                       m_B2*sigmabar2*(2.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) +
                                       2.0*m_B*m_v*sigmabar*(2.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) +
                                       (m_v2 - q2)*(-(m_v2*t*tbar*u) + q2*(2.0 + (-4.0 - 9.0*t2 + 9.0*t)*u)))) +
                                       4.0*m_B*u2*sigmabar*(2.0*(m_v2 - q2)*(-1 + 2.0*(1 + t2 - t)*u) +
                                       m_B2*sigmabar2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       2.0*m_B*m_v*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))*pow(omega,2.0) +
                                       2.0*omega*u*(-(4.0*m_B3*m_v*sigmabar3*t*tbar*u) + 4.0*m_B4*sigmabar4*t*tbar*u -
                                       4.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) +
                                       3.0*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0) +
                                       m_B2*sigmabar2*(m_v2 - q2)*(-1 + 2.0*u*pow(1 - 2.0*t,2.0)))))/4.0;

                    return C_2 * phi_double_bar_4;
                }

                inline
                double I3_V2_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_3 = (m_v*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*pow(m_B,-2.0)*
                                       (-(m_B*sigmabar*(-(2.0*m_B3*m_v*sigmabar3*t*tbar*u) - 2.0*m_B4*sigmabar4*t*tbar*u +
                                       2.0*m_B*m_v*(m_v2 - q2)*sigmabar*t*tbar*u -
                                       (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B2*sigmabar2*(-(2.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)))) +
                                       4.0*m_B*u2*sigmabar*(4.0*m_B*m_v*sigmabar*t*tbar*u +
                                       m_B2*sigmabar2*(-1 + 2.0*(1 + t2 - t)*u) + (m_v2 - q2)*(-1 + 2.0*(1 + t2 - t)*u))*
                                       pow(omega,2.0) - 2.0*omega*u*(6.0*m_B3*m_v*sigmabar3*t*tbar*u -
                                       4.0*m_B2*sigmabar2*q2*t*tbar*u - 2.0*m_B*m_v*(m_v2 - q2)*sigmabar*t*tbar*u +
                                       m_B4*sigmabar4*(-1 + 2.0*u) - (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)))*
                                       pow(sigmabar,-6.0))/4.0;

                    return C_3 * phi_double_bar_4;
                }

                inline
                double I3d1_V2_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);
                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_3 = (m_v*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                       (-(sigmabar*(-(2.0*m_B3*m_v*sigmabar3*t*tbar*u) - 2.0*m_B4*sigmabar4*t*tbar*u +
                                       2.0*m_B*m_v*(m_v2 - q2)*sigmabar*t*tbar*u -
                                       (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B2*sigmabar2*(-(2.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)))) +
                                       4.0*u2*sigmabar*(4.0*m_B*m_v*sigmabar*t*tbar*u +
                                       m_B2*sigmabar2*(-1 + 2.0*(1 + t2 - t)*u) + (m_v2 - q2)*(-1 + 2.0*(1 + t2 - t)*u))*
                                       pow(omega,2.0) - 2.0*omega*u*pow(m_B,-1)*
                                       (6.0*m_B3*m_v*sigmabar3*t*tbar*u - 4.0*m_B2*sigmabar2*q2*t*tbar*u -
                                       2.0*m_B*m_v*(m_v2 - q2)*sigmabar*t*tbar*u + m_B4*sigmabar4*(-1 + 2.0*u) -
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)))*pow(sigmabar,-6.0))/4.0;

                    const double C_3d1 = -(m_v*pow(m_B,-2.0)*(4.0*m_B*omega*sigmabar*
                                         (m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*u*
                                         (-(9.0*m_B2*m_v*sigmabar2*t*tbar*u) + m_v*(m_v2 - q2)*t*tbar*u +
                                         4.0*m_B*q2*sigmabar*t*tbar*u - 2.0*m_B3*sigmabar3*(-1 + 2.0*u)) +
                                         2.0*m_B2*sigmabar2*(m_v - m_B*sigmabar)*
                                         (-(2.0*m_B3*m_v*sigmabar3*t*tbar*u) - 2.0*m_B4*sigmabar4*t*tbar*u +
                                         2.0*m_B*m_v*(m_v2 - q2)*sigmabar*t*tbar*u -
                                         (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                         m_B2*sigmabar2*(-(2.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))) +
                                         5.0*m_B*sigmabar*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                         (-(2.0*m_B3*m_v*sigmabar3*t*tbar*u) - 2.0*m_B4*sigmabar4*t*tbar*u +
                                         2.0*m_B*m_v*(m_v2 - q2)*sigmabar*t*tbar*u -
                                         (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                         m_B2*sigmabar2*(-(2.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))) +
                                         2.0*m_B2*sigmabar2*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                         (3.0*m_B2*m_v*sigmabar2*t*tbar*u + 4.0*m_B3*sigmabar3*t*tbar*u -
                                         m_v*(m_v2 - q2)*t*tbar*u + m_B*sigmabar*
                                         (-(2.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u))) -
                                         8.0*m_B2*sigmabar2*u2*(m_v - m_B*sigmabar)*
                                         (4.0*m_B*m_v*sigmabar*t*tbar*u + m_B2*sigmabar2*(-1 + 2.0*(1 + t2 - t)*u) +
                                         (m_v2 - q2)*(-1 + 2.0*(1 + t2 - t)*u))*pow(omega,2.0) -
                                         20.0*m_B*u2*sigmabar*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                         (4.0*m_B*m_v*sigmabar*t*tbar*u + m_B2*sigmabar2*(-1 + 2.0*(1 + t2 - t)*u) +
                                         (m_v2 - q2)*(-1 + 2.0*(1 + t2 - t)*u))*pow(omega,2.0) -
                                         8.0*m_B2*sigmabar2*u2*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                         (-(2.0*m_v*t*tbar*u) - m_B*sigmabar*(-1 + 2.0*(1 + t2 - t)*u))*pow(omega,2.0) +
                                         4.0*m_B*omega*sigmabar*(m_v - m_B*sigmabar)*u*
                                         (6.0*m_B3*m_v*sigmabar3*t*tbar*u - 4.0*m_B2*sigmabar2*q2*t*tbar*u -
                                         2.0*m_B*m_v*(m_v2 - q2)*sigmabar*t*tbar*u + m_B4*sigmabar4*(-1 + 2.0*u) -
                                         (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)) +
                                         12.0*omega*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*u*
                                         (6.0*m_B3*m_v*sigmabar3*t*tbar*u - 4.0*m_B2*sigmabar2*q2*t*tbar*u -
                                         2.0*m_B*m_v*(m_v2 - q2)*sigmabar*t*tbar*u + m_B4*sigmabar4*(-1 + 2.0*u) -
                                         (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)))*pow(sigmabar,-7.0))/4.0;

                    return C_3 * phi_bar_4 + C_3d1 * phi_double_bar_4;
                }

                inline
                double I1_V2_psi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*(-(m_B*sigmabar*(-(m_B4*sigmabar4*t*tbar*u) +
                                       (m_v2 - q2)*(q2 + q2*(-2.0 - 5.0*t2 + 5.0*t)*u - m_v2*t*tbar*u) -
                                       m_B2*sigmabar2*(-(2.0*m_v2*t*tbar*u) + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)))) +
                                       4.0*m_B*u2*(m_v2 + m_B2*sigmabar2 - q2)*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*
                                       pow(omega,2.0) - omega*u*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*
                                       (m_B4*sigmabar4 + 2.0*m_B2*sigmabar2*(m_v2 + q2) - 3.0*pow(m_v2 - q2,2.0)))*
                                       pow(sigmabar,-4.0))/4.0;

                    return C_1 * psi_bar_4;
                }

                inline
                double I2_V2_psi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

                    const double C_2 = ((-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_B,-2.0)*
                                       (2.0*m_B2*omega*sigmabar2*u + 2.0*omega*(-m_v2 + q2)*u -
                                       m_B*sigmabar*(q2 + 4.0*u2*pow(omega,2.0)))*
                                       (m_B4*sigmabar4 - 2.0*m_B2*sigmabar2*(m_v2 + q2) + pow(m_v2 - q2,2.0))*pow(sigmabar,-5.0))/
                                       8.0;

                    return C_2 * psi_bar_4;
                }

                inline
                double I1_V2_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_double_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_1 = -(m_v*pow(m_B,-2.0)*(-(2.0*omega*(-(3.0*m_v2) + m_B2*sigmabar2 + 3.0*q2)*u*
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) -
                                       m_B*sigmabar*(q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u +
                                       2.0*(-m_v2 + m_B2*sigmabar2)*t*tbar*u) +
                                       4.0*m_B*u2*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(omega,2.0))*pow(sigmabar,-4.0))
                                       /4.0;

                    return C_1 * psi_double_bar_4;
                }

                inline
                double I2_V2_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_double_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_2 = -(m_v*pow(m_B,-2.0)*(8.0*m_B*u2*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*
                                       pow(omega,2.0) + 2.0*omega*u*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*
                                       (-(4.0*m_B2*m_v2*sigmabar2) + m_B4*sigmabar4 + 3.0*pow(m_v2 - q2,2.0)) -
                                       m_B*sigmabar*(-(m_v4*t*tbar*u) - m_B4*sigmabar4*t*tbar*u + 2.0*m_B2*sigmabar2*q2*t*tbar*u -
                                       2.0*m_v2*(-(m_B2*sigmabar2*t*tbar*u) + q2*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) +
                                       (-2.0 + (4.0 + 9.0*t2 - 9.0*t)*u)*pow(q2,2.0)))*pow(sigmabar,-5.0))/4.0;

                    return C_2 * psi_double_bar_4;
                }

                inline
                double I3_V2_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_B,-2.0)*
                                       (-(2.0*m_B2*omega*sigmabar2*u) + 2.0*omega*(m_v2 - q2)*u +
                                       m_B*sigmabar*(q2 + 4.0*u2*pow(omega,2.0)))*
                                       (m_B4*sigmabar4 - 2.0*m_B2*sigmabar2*(m_v2 + q2) + pow(m_v2 - q2,2.0))*pow(sigmabar,-6.0))/
                                       4.0;

                    return C_3 * psi_double_bar_4;
                }

                inline
                double I3d1_V2_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;

                    const double m_v2      = pow(m_v, 2), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4), m_B5 = pow(m_B, 5);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4), sigmabar5 = pow(sigmabar, 5);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);
                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_3 = (m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_B,-1)*
                                       (2.0*m_B2*omega*sigmabar2*u + 2.0*omega*(-m_v2 + q2)*u -
                                       m_B*sigmabar*(q2 + 4.0*u2*pow(omega,2.0)))*
                                       (m_B4*sigmabar4 - 2.0*m_B2*sigmabar2*(m_v2 + q2) + pow(m_v2 - q2,2.0))*pow(sigmabar,-6.0))/
                                       4.0;

                    const double C_3d1 = -(m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_B,-2.0)*
                                         (4.0*m_B4*omega*sigmabar4*(3.0*m_v2 + q2)*u +
                                         m_B5*sigmabar5*(q2 + 4.0*u2*pow(omega,2.0)) -
                                         6.0*m_B3*sigmabar3*(m_v2 + q2)*(q2 + 4.0*u2*pow(omega,2.0)) +
                                         5.0*m_B*sigmabar*(q2 + 4.0*u2*pow(omega,2.0))*pow(m_v2 - q2,2.0) +
                                         12.0*omega*u*pow(m_v2 - q2,3.0) -
                                         8.0*m_B2*omega*sigmabar2*u*(3.0*m_v4 - 2.0*m_v2*q2 - pow(q2,2.0)))*pow(sigmabar,-7.0))/
                                         4.0;

                    return C_3 * psi_bar_4 + C_3d1 * psi_double_bar_4;
                }

                inline
                double I1_V2_chi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

                    const double C_1 = -(pow(sigmabar,-3.0)*(-(8.0*m_B2*m_v*sigmabar3*t*tbar*u) - 2.0*m_B3*sigmabar4*t*tbar*u +
                                       2.0*m_B*sigmabar2*(2.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) -
                                       2.0*t*tbar*u*pow(m_B,-1)*pow(m_v2 - q2,2.0) -
                                       m_v*sigmabar*(-(8.0*m_v2*t*tbar*u) + q2*(1 - 2.0*u*pow(1 - 2.0*t,2.0))) -
                                       4.0*u2*sigmabar*pow(omega,2.0)*
                                       (m_v - 2.0*m_v*u - 2.0*m_B*sigmabar*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))) +
                                       4.0*omega*u*pow(m_B,-1)*(4.0*m_B2*m_v*sigmabar2*t*tbar*u +
                                       m_v*(m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       m_B3*sigmabar3*(-1 + 2.0*u*pow(1 - 2.0*t,2.0)) -
                                       m_B*sigmabar*(m_v2*(-2.0 + (4.0 + 8.0*t2 - 8.0*t)*u) +
                                       q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/8.0;

                    return C_1 * chi_bar_4;
                }

                inline
                double I2_V2_chi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

                    const double C_2 = -(m_v*pow(sigmabar,-4.0)*(4.0*u2*sigmabar*
                                       (m_B2*sigmabar2*(-1 + 2.0*u) + (m_v2 - q2)*(-1 + 2.0*u) -
                                       2.0*m_B*m_v*sigmabar*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u))*pow(omega,2.0) +
                                       2.0*omega*u*pow(m_B,-1)*(8.0*m_B2*sigmabar2*(m_v2 + q2)*t*tbar*u +
                                       m_B4*sigmabar4*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) +
                                       2.0*m_B3*m_v*sigmabar3*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) -
                                       2.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) +
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)) -
                                       sigmabar*(-(4.0*m_B4*sigmabar4*t*tbar*u) +
                                       2.0*m_B*m_v*q2*sigmabar*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) +
                                       (m_v2 - q2)*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u - 4.0*m_v2*t*tbar*u) -
                                       m_B2*sigmabar2*(-(8.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/8.0;

                    return C_2 * chi_bar_4;
                }

                inline
                double I1_V2_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_double_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_1 = -(m_v*pow(m_B,-2.0)*pow(sigmabar,-4.0)*
                                       (2.0*omega*u*(m_B2*sigmabar2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       3.0*(m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       4.0*m_B*m_v*sigmabar*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) +
                                       4.0*m_B*u2*sigmabar*(-1 + 2.0*u)*pow(omega,2.0) -
                                       m_B*sigmabar*(q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u -
                                       2.0*t*tbar*u*pow(m_v - m_B*sigmabar,2.0))))/4.0;

                    return C_1 * chi_double_bar_4;
                }

                inline
                double I2_V2_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_double_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_2 = (m_v*pow(m_B,-2.0)*pow(sigmabar,-5.0)*(-(m_B*sigmabar*
                                       (8.0*m_B3*m_v*sigmabar3*t*tbar*u - 7.0*m_B4*sigmabar4*t*tbar*u +
                                       2.0*m_B2*sigmabar2*(3.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u)) -
                                       4.0*m_B*m_v*sigmabar*(2.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) -
                                       (m_v2 - q2)*(-(m_v2*t*tbar*u) + q2*(2.0 + (-4.0 - 9.0*t2 + 9.0*t)*u)))) -
                                       8.0*m_B*u2*sigmabar*((m_v2 - q2)*(-1 + 2.0*u) +
                                       m_B2*sigmabar2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       2.0*m_B*m_v*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))*pow(omega,2.0) -
                                       2.0*omega*u*(-(8.0*m_B3*m_v*sigmabar3*t*tbar*u) -
                                       8.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 5.0*t2 - 5.0*t)*u) -
                                       m_B4*sigmabar4*(-1 + 2.0*(1 + 6.0*t2 - 6.0*t)*u) +
                                       3.0*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0) +
                                       2.0*m_B2*sigmabar2*(m_v2*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) +
                                       q2*(1 - 2.0*u*pow(1 - 2.0*t,2.0))))))/4.0;

                    return C_2 * chi_double_bar_4;
                }

                inline
                double I3_V2_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(m_v*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*pow(m_B,-2.0)*pow(sigmabar,-6.0)*
                                       (2.0*omega*u*(8.0*m_B2*sigmabar2*q2*t*tbar*u +
                                       m_B4*sigmabar4*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) +
                                       2.0*m_B3*m_v*sigmabar3*(-1 + 2.0*(1 + 5.0*t2 - 5.0*t)*u) -
                                       2.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)) +
                                       4.0*m_B*u2*sigmabar*pow(omega,2.0)*
                                       (m_B2*sigmabar2*(-1 + 2.0*u) + (m_v2 - q2)*(-1 + 2.0*u) -
                                       2.0*m_B*m_v*sigmabar*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))) -
                                       m_B*sigmabar*(-(4.0*m_B3*m_v*sigmabar3*t*tbar*u) - 4.0*m_B4*sigmabar4*t*tbar*u -
                                       (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       2.0*m_B*m_v*sigmabar*(q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u - 2.0*m_v2*t*tbar*u) -
                                       m_B2*sigmabar2*(-(4.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/4.0;

                    return C_3 * chi_double_bar_4;
                }

                inline
                double I3d1_V2_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3), sigmabar4 = pow(sigmabar, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);
                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_3 = (m_v*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*pow(sigmabar,-6.0)*
                                       (-(2.0*omega*u*pow(m_B,-1)*(8.0*m_B2*sigmabar2*q2*t*tbar*u +
                                       m_B4*sigmabar4*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) +
                                       2.0*m_B3*m_v*sigmabar3*(-1 + 2.0*(1 + 5.0*t2 - 5.0*t)*u) -
                                       2.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0))) -
                                       4.0*u2*sigmabar*pow(omega,2.0)*(m_B2*sigmabar2*(-1 + 2.0*u) + (m_v2 - q2)*(-1 + 2.0*u) -
                                       2.0*m_B*m_v*sigmabar*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))) +
                                       sigmabar*(-(4.0*m_B3*m_v*sigmabar3*t*tbar*u) - 4.0*m_B4*sigmabar4*t*tbar*u -
                                       (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                       2.0*m_B*m_v*sigmabar*(q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u - 2.0*m_v2*t*tbar*u) -
                                       m_B2*sigmabar2*(-(4.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/4.0;

                    const double C_3d1 = -(m_v*pow(m_B,-2.0)*pow(sigmabar,-7.0)*
                                         (4.0*m_B*omega*sigmabar*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*u*
                                         (-(8.0*m_B*q2*sigmabar*t*tbar*u) - 2.0*m_B3*sigmabar3*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) -
                                         3.0*m_B2*m_v*sigmabar2*(-1 + 2.0*(1 + 5.0*t2 - 5.0*t)*u) +
                                         m_v*(m_v2 - q2)*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) +
                                         4.0*m_B*omega*sigmabar*(m_v - m_B*sigmabar)*u*
                                         (8.0*m_B2*sigmabar2*q2*t*tbar*u + m_B4*sigmabar4*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) +
                                         2.0*m_B3*m_v*sigmabar3*(-1 + 2.0*(1 + 5.0*t2 - 5.0*t)*u) -
                                         2.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                         (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)) +
                                         12.0*omega*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*u*
                                         (8.0*m_B2*sigmabar2*q2*t*tbar*u + m_B4*sigmabar4*(1 + (-2.0 + 4.0*t2 - 4.0*t)*u) +
                                         2.0*m_B3*m_v*sigmabar3*(-1 + 2.0*(1 + 5.0*t2 - 5.0*t)*u) -
                                         2.0*m_B*m_v*(m_v2 - q2)*sigmabar*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                         (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_v2 - q2,2.0)) +
                                         8.0*m_B2*sigmabar2*u2*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                         pow(omega,2.0)*(-(m_B*sigmabar*(-1 + 2.0*u)) + m_v*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))) +
                                         8.0*m_B2*sigmabar2*u2*(m_v - m_B*sigmabar)*pow(omega,2.0)*
                                         (m_B2*sigmabar2*(-1 + 2.0*u) + (m_v2 - q2)*(-1 + 2.0*u) -
                                         2.0*m_B*m_v*sigmabar*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))) +
                                         20.0*m_B*u2*sigmabar*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                         pow(omega,2.0)*(m_B2*sigmabar2*(-1 + 2.0*u) + (m_v2 - q2)*(-1 + 2.0*u) -
                                         2.0*m_B*m_v*sigmabar*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))) -
                                         2.0*m_B2*sigmabar2*(m_v - m_B*sigmabar)*
                                         (-(4.0*m_B3*m_v*sigmabar3*t*tbar*u) - 4.0*m_B4*sigmabar4*t*tbar*u -
                                         (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                         2.0*m_B*m_v*sigmabar*(q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u - 2.0*m_v2*t*tbar*u) -
                                         m_B2*sigmabar2*(-(4.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0)))) -
                                         5.0*m_B*sigmabar*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                         (-(4.0*m_B3*m_v*sigmabar3*t*tbar*u) - 4.0*m_B4*sigmabar4*t*tbar*u -
                                         (m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) -
                                         2.0*m_B*m_v*sigmabar*(q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u - 2.0*m_v2*t*tbar*u) -
                                         m_B2*sigmabar2*(-(4.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0)))) -
                                         2.0*m_B2*sigmabar2*(m_v2 + m_B2*sigmabar2 - q2 - 2.0*m_B*m_v*sigmabar)*
                                         (6.0*m_B2*m_v*sigmabar2*t*tbar*u + 8.0*m_B3*sigmabar3*t*tbar*u +
                                         m_v*(q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u - 2.0*m_v2*t*tbar*u) +
                                         m_B*sigmabar*(-(4.0*m_v2*t*tbar*u) + q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/4.0;

                    return C_3 * chi_bar_4 + C_3d1 * chi_double_bar_4;
                }

                double integrand_V2(const std::array<double, 4> & args, const double & q2) const
                {
                    const double sigma = args[0];
                    const double x     = args[1];
                    const double t     = args[2];
                    const double u     = args[3];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = f_B() * m_B / ((8.0  * M_PI * M_PI * f_V() * m_V()) * (m_B * m_B - m_V() * m_V()));
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_V2 = pow(m_V(), 2);
                    const double M4   = pow(M2, 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_V2) / M2());

                    const double I1 = I1_V2_phi_3(sigma, omega, t, u, q2)             + I1_V2_phi_4(sigma, omega, t, u, q2)
                                    + I1_V2_phi_bar_3(sigma, omega, t, u, q2)         + I1_V2_phi_bar_4(sigma, omega, t, u, q2)
                                    + I1_V2_phi_double_bar_3(sigma, omega, t, u, q2)  + I1_V2_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I1_V2_psi_bar_4(sigma, omega, t, u, q2)         + I1_V2_chi_bar_4(sigma, omega, t, u, q2)
                                    + I1_V2_psi_double_bar_4(sigma, omega, t, u, q2)  + I1_V2_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I2 = I2_V2_phi_bar_3(sigma, omega, t, u, q2)         + I2_V2_phi_bar_4(sigma, omega, t, u, q2)
                                    + I2_V2_phi_double_bar_3(sigma, omega, t, u, q2)  + I2_V2_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I2_V2_psi_bar_4(sigma, omega, t, u, q2)         + I2_V2_chi_bar_4(sigma, omega, t, u, q2)
                                    + I2_V2_psi_double_bar_4(sigma, omega, t, u, q2)  + I2_V2_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3 = I3_V2_phi_double_bar_3(sigma, omega, t, u, q2)  + I3_V2_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I3_V2_psi_double_bar_4(sigma, omega, t, u, q2)  + I3_V2_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result = 0.0;
                    result += - I1;
                    result +=   I2 / M2;
                    result += - I3 / (2.0 * M4);
                    result *=   prefactor * exp;

                    return result;
                }

                double surface_V2(const std::array<double, 3> & args, const double & sigma, const double & q2) const
                {
                    const double x     = args[0];
                    const double t     = args[1];
                    const double u     = args[2];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = f_B() * m_B / ((8.0  * M_PI * M_PI * f_V() * m_V()) * (m_B * m_B - m_V() * m_V()));
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_B2 = pow(m_B(), 2);
                    const double m_v2 = pow(m_v(), 2);
                    const double m_V2 = pow(m_V(), 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_V2) / M2());

                    const double sigmabar = 1.0 - sigma, sigmabar2 = pow(sigmabar, 2);
                    const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
                    const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;

                    const double I2   = I2_V2_phi_bar_3(sigma, omega, t, u, q2)           + I2_V2_phi_bar_4(sigma, omega, t, u, q2)
                                      + I2_V2_phi_double_bar_3(sigma, omega, t, u, q2)    + I2_V2_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I2_V2_psi_bar_4(sigma, omega, t, u, q2)           + I2_V2_chi_bar_4(sigma, omega, t, u, q2)
                                      + I2_V2_psi_double_bar_4(sigma, omega, t, u, q2)    + I2_V2_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3   = I3_V2_phi_double_bar_3(sigma, omega, t, u, q2)    + I3_V2_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3_V2_psi_double_bar_4(sigma, omega, t, u, q2)    + I3_V2_chi_double_bar_4(sigma, omega, t, u, q2);

                    const double I3d1 = I3d1_V2_phi_double_bar_3(sigma, omega, t, u, q2)  + I3d1_V2_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3d1_V2_psi_double_bar_4(sigma, omega, t, u, q2)  + I3d1_V2_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result = 0.0;
                    result +=  eta * I2 / m_B2;
                    result += -0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
                    result *=  prefactor * exp;

                    return result;
                }

                double integrand_V2_m1(const std::array<double, 4> & args, const double & q2) const
                {
                    const double sigma = args[0];
                    const double x     = args[1];
                    const double t     = args[2];
                    const double u     = args[3];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = f_B() * m_B / ((8.0  * M_PI * M_PI * f_V() * m_V()) * (m_B * m_B - m_V() * m_V()));
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_V2 = pow(m_V(), 2);
                    const double M4   = pow(M2, 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_V2) / M2());

                    const double I1 = I1_V2_phi_3(sigma, omega, t, u, q2)             + I1_V2_phi_4(sigma, omega, t, u, q2)
                                    + I1_V2_phi_bar_3(sigma, omega, t, u, q2)         + I1_V2_phi_bar_4(sigma, omega, t, u, q2)
                                    + I1_V2_phi_double_bar_3(sigma, omega, t, u, q2)  + I1_V2_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I1_V2_psi_bar_4(sigma, omega, t, u, q2)         + I1_V2_chi_bar_4(sigma, omega, t, u, q2)
                                    + I1_V2_psi_double_bar_4(sigma, omega, t, u, q2)  + I1_V2_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I2 = I2_V2_phi_bar_3(sigma, omega, t, u, q2)         + I2_V2_phi_bar_4(sigma, omega, t, u, q2)
                                    + I2_V2_phi_double_bar_3(sigma, omega, t, u, q2)  + I2_V2_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I2_V2_psi_bar_4(sigma, omega, t, u, q2)         + I2_V2_chi_bar_4(sigma, omega, t, u, q2)
                                    + I2_V2_psi_double_bar_4(sigma, omega, t, u, q2)  + I2_V2_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3 = I3_V2_phi_double_bar_3(sigma, omega, t, u, q2)  + I3_V2_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I3_V2_psi_double_bar_4(sigma, omega, t, u, q2)  + I3_V2_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result1 =    0.0;
                           result1 += - I1;
                           result1 +=   I2 / M2;
                           result1 += - I3 / (2.0 * M4);
                           result1 *=   exp * s(sigma, q2);

                    double result2 =    0.0;
                           result2 += - I2;
                           result2 +=   I3 / M2;
                           result2 *=   exp;

                    return prefactor * (result1 + result2);
                }

                double surface_V2_m1(const std::array<double, 3> & args, const double & sigma, const double & q2) const
                {
                    const double x     = args[0];
                    const double t     = args[1];
                    const double u     = args[2];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = f_B() * m_B / ((8.0  * M_PI * M_PI * f_V() * m_V()) * (m_B * m_B - m_V() * m_V()));
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_B2 = pow(m_B(), 2);
                    const double m_v2 = pow(m_v(), 2);
                    const double m_V2 = pow(m_V(), 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_V2) / M2());

                    const double sigmabar = 1.0 - sigma, sigmabar2 = pow(sigmabar, 2);
                    const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
                    const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;

                    const double I2   = I2_V2_phi_bar_3(sigma, omega, t, u, q2)           + I2_V2_phi_bar_4(sigma, omega, t, u, q2)
                                      + I2_V2_phi_double_bar_3(sigma, omega, t, u, q2)    + I2_V2_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I2_V2_psi_bar_4(sigma, omega, t, u, q2)           + I2_V2_chi_bar_4(sigma, omega, t, u, q2)
                                      + I2_V2_psi_double_bar_4(sigma, omega, t, u, q2)    + I2_V2_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3   = I3_V2_phi_double_bar_3(sigma, omega, t, u, q2)    + I3_V2_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3_V2_psi_double_bar_4(sigma, omega, t, u, q2)    + I3_V2_chi_double_bar_4(sigma, omega, t, u, q2);

                    const double I3d1 = I3d1_V2_phi_double_bar_3(sigma, omega, t, u, q2)  + I3d1_V2_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3d1_V2_psi_double_bar_4(sigma, omega, t, u, q2)  + I3d1_V2_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result1 =   0.0;
                           result1 +=  eta * I2 / m_B2;
                           result1 += -0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
                           result1 *=  exp * s(sigma, q2);

                    double result2 =  0.0;
                           result2 += 0.5 * eta * I3 / m_B2;
                           result2 *= exp;

                    return prefactor * (result1 + result2);
                }

                double V2(const double & q2) const
                {
                    const double sigma_0 = this->sigma_0(q2, s0_0_V2(), s0_1_V2());

                    double value_integral = 0.0;
                    double value_surface  = 0.0;

                    const std::function<double (const std::array<double, 4> &)> integrand = std::bind(&LCSR::integrand_V2, this, std::placeholders::_1, q2);
                    const std::function<double (const std::array<double, 3> &)> surface = std::bind(&LCSR::surface_V2, this, std::placeholders::_1, sigma_0, q2);

                    value_integral = integrate(integrand, { 0.0, 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0, 1.0 }, cubature::Config());
                    value_surface  = integrate(surface, { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, cubature::Config());

                    return value_integral + value_surface;
                }

                double normalized_first_moment_V2(const double & q2) const
                {
                    const double sigma_0 = this->sigma_0(q2, s0_0_V2(), s0_1_V2());

                    double value_integral = 0.0;
                    double value_surface  = 0.0;

                    const std::function<double (const std::array<double, 4> &)> integrand = std::bind(&LCSR::integrand_V2, this, std::placeholders::_1, q2);
                    const std::function<double (const std::array<double, 3> &)> surface = std::bind(&LCSR::surface_V2, this, std::placeholders::_1, sigma_0, q2);

                    value_integral = integrate(integrand, { 0.0, 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0, 1.0 }, cubature::Config());
                    value_surface  = integrate(surface, { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, cubature::Config());

                    const double denominator = value_integral + value_surface;

                    double value_integral_m1 = 0.0;
                    double value_surface_m1  = 0.0;

                    const std::function<double (const std::array<double, 4> &)> integrand_m1 = std::bind(&LCSR::integrand_V2_m1, this, std::placeholders::_1, q2);
                    const std::function<double (const std::array<double, 3> &)> surface_m1 = std::bind(&LCSR::surface_V2_m1, this, std::placeholders::_1, sigma_0, q2);

                    value_integral_m1 = integrate(integrand_m1, { 0.0, 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0, 1.0 }, cubature::Config());
                    value_surface_m1  = integrate(surface_m1, { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, cubature::Config());

                    const double numerator = value_integral_m1 + value_surface_m1;

                    return numerator / denominator;
                }

                inline
                double I1_V23_phi_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_3 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_3 = this->phi_3(omega_1, omega_2);

                    const double C_1 = (pow(m_B,-2.0)*pow(sigmabar,-2.0)*(m_B*
                                       (-(2.0*m_B3*sigma3*t*tbar*(-3.0 + 2.0*u)) +
                                       m_B*sigma2*(q2*(3.0 - 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) +
                                       4.0*m_B*t*tbar*(m_v + m_B*(-3.0 + 2.0*u))) -
                                       2.0*m_v*(-(t*tbar*(m_B2 - 2.0*m_v2 + m_B*m_v*(3.0 - 2.0*u))) +
                                       q2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       sigma*(-(6.0*m_B2*m_v*t*tbar) - 2.0*m_B3*t*tbar*(-3.0 + 2.0*u) +
                                       2.0*m_v*q2*(1 + (-2.0 - 6.0*t2 + 6.0*t)*u) +
                                       m_B*(2.0*m_v2*t*tbar*(-3.0 + 2.0*u) + q2*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u)))
                                       ) + 4.0*m_B*u2*sigmabar*(2.0*m_v*(1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u) +
                                       m_B*sigma*(-3.0 + 6.0*u + t*(6.0 - 20.0*u) + t2*(-6.0 + 20.0*u)))*pow(omega,2.0) -
                                       2.0*omega*u*(m_B3*sigma3*(-3.0 + 6.0*u + t*(6.0 - 20.0*u) + t2*(-6.0 + 20.0*u)) -
                                       2.0*m_B2*sigma2*(m_B*(-3.0 + t*(9.0 - 22*u) + 6.0*u + t2*(-9.0 + 22*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) -
                                       2.0*m_v*(m_B*m_v*t*tbar*(-3.0 + 2.0*u) +
                                       (m_v2 - q2)*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_B2*(1 - 2.0*u + t2*(4.0 - 6.0*u) + t*(-4.0 + 6.0*u))) +
                                       m_B*sigma*(-(4.0*m_B*m_v*(1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u)) +
                                       q2*(-3.0 + 2.0*(3.0 + 8.0*t2 - 8.0*t)*u) +
                                       m_v2*(3.0 - 6.0*u + t2*(6.0 - 20.0*u) + t*(-6.0 + 20.0*u)) +
                                       3.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))))/4.0;

                    return C_1 * phi_3;
                }

                inline
                double I1_V23_phi_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_3 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigma2    = pow(sigma, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*(-(m_B*sigmabar*(q2*(-1 + (2.0 + 4.0*t2*(1 + sigma) - 4.0*(1 + sigma)*t)*u) -
                                       2.0*t*tbar*(m_B*m_v*(3.0 - 4.0*sigma) - m_v2*(1 + 2.0*u) -
                                       m_B2*sigmabar*(1 - 2.0*u + sigma*(-5.0 + 6.0*u))))) -
                                       2.0*omega*u*(-(q2*(1 + sigma)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) -
                                       4.0*m_B*m_v*sigmabar*(-(t*tbar) +
                                       sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       m_v2*(-1 + 2.0*t2 - 2.0*t + 2.0*u +
                                       sigma*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))) -
                                       m_B2*sigmabar*(-1 + 2.0*u + t*(4.0 - 12.0*u) + 4.0*t2*(-1 + 3.0*u) +
                                       sigma2*(-5.0 + t*(10.0 - 32*u) + 10.0*u + 2.0*t2*(-5.0 + 16.0*u)) -
                                       2.0*sigma*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)))) -
                                       4.0*u2*sigmabar*(2.0*m_v*t*tbar +
                                       m_B*(-1 + 2.0*u + t*(2.0 - 8.0*u - 4.0*sigma*(-2.0 + u)) +
                                       t2*(-2.0 + 8.0*u + 4.0*sigma*(-2.0 + u))))*pow(omega,2.0))*pow(sigmabar,-3.0))/4.0;

                    return C_1 * phi_bar_3;
                }

                inline
                double I2_V23_phi_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_3 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);

                    const double C_2 = (pow(m_B,-2.0)*pow(sigmabar,-4.0)*(-(2.0*omega*u*
                                       (m_B4*sigma5*(-5.0 + t*(10.0 - 32*u) + 10.0*u + 2.0*t2*(-5.0 + 16.0*u)) +
                                       2.0*m_v*t*tbar*(2.0*m_B3 - 2.0*m_B*(m_v2 + q2) - m_v*(m_v2 - q2)*(-1 + 2.0*u) +
                                       m_B2*m_v*(9.0 - 6.0*u)) -
                                       2.0*m_B3*sigma4*(m_B*(-8.0 + t*(25 - 58*u) + 16.0*u + t2*(-25 + 58*u)) +
                                       2.0*m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       2.0*m_B2*sigma3*(-(6.0*m_B*m_v*(1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u)) +
                                       3.0*m_B2*(-3.0 + t*(15.0 - 26*u) + 6.0*u + t2*(-15.0 + 26*u)) -
                                       2.0*m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       q2*(-2.0 + 4.0*u + t2*(5.0 + 6.0*u) - t*(5.0 + 6.0*u))) +
                                       sigma*(m_B4*(-1 + 20.0*t2*(-1 + u) - 20.0*t*(-1 + u) + 2.0*u) +
                                       2.0*m_B2*t*tbar*(m_v2*(-17.0 + 10.0*u) + 2.0*q2*(-5.0 + 4.0*u)) -
                                       4.0*m_B3*m_v*(-1 + 2.0*u + t*(7.0 - 6.0*u) + t2*(-7.0 + 6.0*u)) +
                                       4.0*m_B*m_v*(m_v2*(1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u) +
                                       q2*(1 + t - 2.0*u + 6.0*t*u - t2*(1 + 6.0*u))) -
                                       (m_v2 - q2)*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u +
                                       m_v2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)))) -
                                       2.0*m_B*sigma2*(m_B3*(-4.0 + t*(35 - 46*u) + 8.0*u + t2*(-35 + 46*u)) +
                                       2.0*m_B2*m_v*(3.0 - 6.0*u + t2*(13.0 - 18.0*u) + t*(-13.0 + 18.0*u)) +
                                       2.0*m_v*(q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u +
                                       m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       m_B*(q2*(-2.0 + 4.0*u + t2*(15.0 - 2.0*u) + t*(-15.0 + 2.0*u)) -
                                       2.0*m_v2*(-1 + 2.0*u + t*(6.0 - 8.0*u) + t2*(-6.0 + 8.0*u)))))) +
                                       4.0*u2*sigmabar*(2.0*m_v*t*tbar*(m_B2 + m_v2 - q2 - 2.0*m_B*m_v*(-2.0 + u)) +
                                       m_B3*sigma3*(-5.0 + t*(10.0 - 32*u) + 10.0*u + 2.0*t2*(-5.0 + 16.0*u)) -
                                       2.0*m_B2*sigma2*(m_B*(-3.0 + t*(10.0 - 22*u) + 6.0*u + 2.0*t2*(-5.0 + 11.0*u)) +
                                       m_v*(2.0 - 4.0*u + t2*(5.0 - 12.0*u) + t*(-5.0 + 12.0*u))) +
                                       m_B*sigma*(-(4.0*m_B*m_v*(1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u)) +
                                       m_B2*(-1 + 2.0*u + 2.0*t2*(-5.0 + 6.0*u) - 2.0*t*(-5.0 + 6.0*u)) +
                                       q2*(1 - 2.0*u - 2.0*t2*(-5.0 + 6.0*u) + 2.0*t*(-5.0 + 6.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))))*pow(omega,2.0) -
                                       m_B*sigmabar*(-(2.0*m_B4*sigma4*t*tbar*(-5.0 + 4.0*u)) -
                                       2.0*m_v*t*tbar*(m_B3 + m_v3 - m_B*(3.0*m_v2 + q2) + m_v*q2*(-1 + 2.0*u) +
                                       m_B2*m_v*(5.0 - 4.0*u)) +
                                       m_B2*sigma3*(q2*(5.0 - 2.0*(5.0 + 12.0*t2 - 12.0*t)*u) +
                                       2.0*m_B*t*tbar*(4.0*m_v + 3.0*m_B*(-5.0 + 4.0*u))) +
                                       2.0*m_B*sigma2*(-(9.0*m_B2*m_v*t*tbar) - 3.0*m_B3*t*tbar*(-5.0 + 4.0*u) +
                                       2.0*m_v*q2*(1 + (-2.0 - 6.0*t2 + 6.0*t)*u) +
                                       m_B*(4.0*m_v2*t*(-tbar + u - t*u) +
                                       q2*(-3.0 + 6.0*u + 5.0*t2*(1 + 2.0*u) - 5.0*t*(1 + 2.0*u)))) +
                                       sigma*(12.0*m_B3*m_v*t*tbar + 2.0*m_B4*t*tbar*(-5.0 + 4.0*u) +
                                       4.0*m_B*(-(2.0*m_v3*t*tbar) + m_v*q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)) +
                                       m_B2*(-(2.0*m_v2*t*tbar*(-9.0 + 8.0*u)) +
                                       q2*(1 - 2.0*u + 2.0*t2*(-5.0 + 2.0*u) + t*(10.0 - 4.0*u))) +
                                       q2*(q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) + m_v2*(1 - 2.0*u*pow(1 - 2.0*t,2.0)))))))/4.0;

                    return C_2 * phi_bar_3;
                }

                inline
                double I1_V23_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to phi_double_bar_3 with a single pole in k2
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_1 = -((sigmabar*t*tbar*(1 + sigma - u) + omega*sigma*u*
                                       ((1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       sigma*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)))*pow(m_B,-1) +
                                       4.0*u2*(1 + sigma)*t*tbar*(-1 + u)*pow(m_B,-2.0)*pow(omega,2.0))*pow(sigmabar,-3.0));

                    return C_1 * phi_double_bar_3;
                }

                inline
                double I2_V23_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_3 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_2 = (pow(sigmabar,-4.0)*(-(2.0*omega*u*(m_B2*sigma4*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) -
                                       4.0*t*tbar*(m_B*m_v + q2 + m_B2*(-1 + u) - q2*u + m_v2*(-2.0 + u)) +
                                       sigma*(8.0*m_B*m_v*t*tbar - q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) -
                                       m_B2*(1 + 6.0*t2 - 6.0*t)*(-1 + 2.0*u) +
                                       m_v2*(-3.0 + 6.0*u + 2.0*t2*(1 + 6.0*u) - 2.0*t*(1 + 6.0*u))) +
                                       m_B*sigma3*(m_B*(3.0 - 6.0*u + t2*(2.0 - 12.0*u) + 2.0*t*(-1 + 6.0*u)) +
                                       2.0*m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       sigma2*(-(2.0*m_B*m_v*(-1 + 2.0*u + t*(4.0 - 6.0*u) + t2*(-4.0 + 6.0*u))) +
                                       m_B2*(-3.0 + 6.0*u + t*(2.0 - 16.0*u) + 2.0*t2*(-1 + 8.0*u)) +
                                       q2*(3.0 - 6.0*u + t2*(2.0 - 16.0*u) + 2.0*t*(-1 + 8.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))))*pow(m_B,-1)) +
                                       4.0*u2*(2.0*m_B*sigma2*t*tbar*(m_v + 5.0*m_B*(-1 + u)) -
                                       2.0*t*tbar*(m_B*m_v + m_B2*(-1 + u) + (3.0*m_v2 - q2)*(-1 + u)) +
                                       m_B2*sigma3*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)) +
                                       sigma*(-(2.0*(m_v2 - 3.0*q2)*t*tbar*(-1 + u)) +
                                       m_B2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u))))*pow(m_B,-2.0)*pow(omega,2.0) -
                                       sigmabar*(-(2.0*m_B2*sigma3*t*tbar*(-1 + 2.0*u)) +
                                       2.0*t*tbar*(m_B*m_v + q2 + m_B2*(-1 + u) - q2*u - m_v2*(1 + u)) +
                                       sigma*(q2*(-1 + 2.0*t2 - 2.0*t + 2.0*u) -
                                       2.0*t*tbar*(3.0*m_B*m_v - m_v2*(1 + 2.0*u) + m_B2*(-3.0 + 4.0*u))) +
                                       sigma2*(2.0*m_B*t*tbar*(2.0*m_v + m_B*(-3.0 + 5.0*u)) +
                                       q2*(-1 + 2.0*u*pow(1 - 2.0*t,2.0))))))/2.0;

                    return C_2 * phi_double_bar_3;
                }

                inline
                double I3_V23_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_3 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5), sigma6 = pow(sigma, 6);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_3 = -(pow(sigmabar,-5.0)*(-(sigmabar*(-(2.0*m_B4*sigma5*t*tbar*(-1 + u)) +
                                       m_B2*sigma4*(q2 + 2.0*m_B*t*tbar*(m_v + 4.0*m_B*(-1 + u)) +
                                       q2*(-2.0 - 4.0*t2 + 4.0*t)*u) -
                                       2.0*m_v2*t*tbar*(-m_v2 + q2 + m_B2*(-1 + u) - q2*u - m_B*m_v*(-2.0 + u)) +
                                       m_B*sigma3*(-(6.0*m_B2*m_v*t*tbar) - 12.0*m_B3*t*tbar*(-1 + u) +
                                       m_v*q2*(1 + (-2.0 - 6.0*t2 + 6.0*t)*u) +
                                       2.0*m_B*(m_v2*t*(-tbar + 2.0*u - 2.0*t*u) +
                                       q2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)))) -
                                       sigma2*(-(6.0*m_B3*m_v*t*tbar) - 8.0*m_B4*t*tbar*(-1 + u) +
                                       q2*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u - 2.0*m_v2*t*tbar*u) -
                                       m_B2*(-(2.0*m_v2*t*tbar*(-3.0 + 5.0*u)) +
                                       q2*(1 - 2.0*u + 4.0*t2*(-2.0 + u) - 4.0*t*(-2.0 + u))) +
                                       m_B*(2.0*m_v3*t*tbar + m_v*q2*
                                       (1 - 2.0*u - 2.0*t2*(1 + 3.0*u) + t*(2.0 + 6.0*u)))) +
                                       sigma*(-(2.0*m_B*m_v*(m_B2 - q2)*t*tbar) - 2.0*m_v4*t*tbar*u -
                                       2.0*m_B*m_v3*t*tbar*(-3.0 + u) +
                                       m_v2*(-(q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)) + 2.0*m_B2*t*tbar*(-3.0 + 4.0*u)) -
                                       2.0*t*tbar*(-1 + u)*pow(m_B2 - q2,2.0)))) -
                                       2.0*omega*u*pow(m_B,-1)*(2.0*m_v2*t*tbar*
                                       (2.0*m_B2*(-1 + u) - 2.0*q2*(-1 + u) - 2.0*m_B*m_v*(-2.0 + u) + m_v2*(-2.0 + u)) +
                                       m_B4*sigma6*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       m_B3*sigma5*(m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       4.0*m_B*(-1 + 2.0*u + t*(3.0 - 7.0*u) + t2*(-3.0 + 7.0*u))) +
                                       m_B2*sigma4*(4.0*q2*t*(-tbar + u - t*u) +
                                       2.0*m_B2*(-3.0 + t*(14.0 - 26*u) + 6.0*u + 2.0*t2*(-7.0 + 13.0*u)) +
                                       m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_B*m_v*(3.0 - 6.0*u - 2.0*t2*(-5.0 + 9.0*u) + 2.0*t*(-5.0 + 9.0*u))) +
                                       m_B*sigma3*(3.0*m_B2*m_v*(-1 + 6.0*t2*(-1 + u) - 6.0*t*(-1 + u) + 2.0*u) -
                                       4.0*m_B3*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_v3*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_v*q2*(-1 + 2.0*u - 2.0*t*(1 + 3.0*u) + t2*(2.0 + 6.0*u)) +
                                       m_B*(-(16.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-3.0 + t*(10.0 - 24*u) + 6.0*u + 2.0*t2*(-5.0 + 12.0*u)))) +
                                       sigma2*(m_B4*(-1 + t*(18.0 - 22*u) + 2.0*u + 2.0*t2*(-9.0 + 11.0*u)) +
                                       m_B3*m_v*(1 - 2.0*u + 2.0*t*(-7.0 + 3.0*u) + t2*(14.0 - 6.0*u)) +
                                       m_B*(m_v3*(-1 + 10.0*t2*(-1 + u) - 10.0*t*(-1 + u) + 2.0*u) +
                                       m_v*q2*(1 - 2.0*u - 6.0*t2*(1 + u) + 6.0*t*(1 + u))) +
                                       m_B2*(20.0*q2*t*tbar*(-1 + u) +
                                       m_v2*(3.0 + t2*(18.0 - 34*u) - 6.0*u + 2.0*t*(-9.0 + 17.0*u))) +
                                       q2*(q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u)))) +
                                       sigma*(4.0*m_B3*m_v*t*tbar + 4.0*m_B4*t*tbar*(-1 + u) +
                                       m_v2*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       4.0*m_B*m_v*t*tbar*(-q2 + 2.0*m_v2*(-2.0 + u)) +
                                       m_v4*(1 - 2.0*u - 2.0*t2*(1 + 2.0*u) + t*(2.0 + 4.0*u)) +
                                       m_B2*(-(8.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-1 + 2.0*u + 2.0*t2*(-7.0 + 10.0*u) + t*(14.0 - 20.0*u))) +
                                       4.0*t*(-tbar + u - t*u)*pow(q2,2.0))) -
                                       4.0*u2*pow(m_B,-2.0)*pow(omega,2.0)*
                                       (-(2.0*m_v2*t*tbar*(m_B2*(-1 + u) + (m_v2 - q2)*(-1 + u) - m_B*m_v*(-2.0 + u))) +
                                       m_B4*sigma5*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       m_B3*sigma4*(m_B*(-3.0 + 6.0*u + t*(8.0 - 20.0*u) + 4.0*t2*(-2.0 + 5.0*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       m_B*sigma2*(2.0*m_v*t*tbar*(-q2 + m_v2*(-1 + u)) +
                                       m_B2*m_v*(-1 + 6.0*t2*(-1 + u) - 6.0*t*(-1 + u) + 2.0*u) +
                                       m_B3*(1 - 2.0*u + t2*(8.0 - 12.0*u) + 4.0*t*(-2.0 + 3.0*u)) +
                                       m_B*(q2*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)))) +
                                       sigma*(-(2.0*m_B3*m_v*t*tbar) - 2.0*m_B4*t*tbar*(-1 + u) +
                                       2.0*(m_v2 - q2)*q2*t*tbar*(-1 + u) - 2.0*m_B*m_v*t*tbar*(-q2 + m_v2*(-3.0 + 2.0*u)) +
                                       m_B2*(4.0*q2*t*(-tbar + u - t*u) +
                                       m_v2*(1 - 2.0*u + t2*(6.0 - 10.0*u) + 2.0*t*(-3.0 + 5.0*u)))) +
                                       m_B2*sigma3*(-(2.0*m_B*m_v*(1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u)) +
                                       2.0*m_v2*t*(-tbar + u - t*u) + 3.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0) -
                                       q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))));

                    return C_3 * phi_double_bar_3;
                }

                inline
                double I3d1_V23_phi_double_bar_3(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_3 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4), m_B5 = pow(m_B, 5), m_B6 = pow(m_B, 6);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5), sigma6 = pow(sigma, 6);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_3 = this->phi_bar_3(omega_1, omega_2);
                    const double phi_double_bar_3 = this->phi_double_bar_3(omega_1, omega_2);

                    const double C_3 = -(pow(sigmabar,-5.0)*(-(m_B*sigmabar*(-(2.0*m_B4*sigma5*t*tbar*(-1 + u)) +
                                       m_B2*sigma4*(q2 + 2.0*m_B*t*tbar*(m_v + 4.0*m_B*(-1 + u)) +
                                       q2*(-2.0 - 4.0*t2 + 4.0*t)*u) -
                                       2.0*m_v2*t*tbar*(-m_v2 + q2 + m_B2*(-1 + u) - q2*u - m_B*m_v*(-2.0 + u)) +
                                       m_B*sigma3*(-(6.0*m_B2*m_v*t*tbar) - 12.0*m_B3*t*tbar*(-1 + u) +
                                       m_v*q2*(1 + (-2.0 - 6.0*t2 + 6.0*t)*u) +
                                       2.0*m_B*(m_v2*t*(-tbar + 2.0*u - 2.0*t*u) +
                                       q2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)))) -
                                       sigma2*(-(6.0*m_B3*m_v*t*tbar) - 8.0*m_B4*t*tbar*(-1 + u) +
                                       q2*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u - 2.0*m_v2*t*tbar*u) -
                                       m_B2*(-(2.0*m_v2*t*tbar*(-3.0 + 5.0*u)) +
                                       q2*(1 - 2.0*u + 4.0*t2*(-2.0 + u) - 4.0*t*(-2.0 + u))) +
                                       m_B*(2.0*m_v3*t*tbar + m_v*q2*
                                       (1 - 2.0*u - 2.0*t2*(1 + 3.0*u) + t*(2.0 + 6.0*u)))) +
                                       sigma*(-(2.0*m_B*m_v*(m_B2 - q2)*t*tbar) - 2.0*m_v4*t*tbar*u -
                                       2.0*m_B*m_v3*t*tbar*(-3.0 + u) +
                                       m_v2*(-(q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)) + 2.0*m_B2*t*tbar*(-3.0 + 4.0*u)) -
                                       2.0*t*tbar*(-1 + u)*pow(m_B2 - q2,2.0)))) -
                                       2.0*omega*u*(2.0*m_v2*t*tbar*(2.0*m_B2*(-1 + u) - 2.0*q2*(-1 + u) -
                                       2.0*m_B*m_v*(-2.0 + u) + m_v2*(-2.0 + u)) +
                                       m_B4*sigma6*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       m_B3*sigma5*(m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       4.0*m_B*(-1 + 2.0*u + t*(3.0 - 7.0*u) + t2*(-3.0 + 7.0*u))) +
                                       m_B2*sigma4*(4.0*q2*t*(-tbar + u - t*u) +
                                       2.0*m_B2*(-3.0 + t*(14.0 - 26*u) + 6.0*u + 2.0*t2*(-7.0 + 13.0*u)) +
                                       m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_B*m_v*(3.0 - 6.0*u - 2.0*t2*(-5.0 + 9.0*u) + 2.0*t*(-5.0 + 9.0*u))) +
                                       m_B*sigma3*(3.0*m_B2*m_v*(-1 + 6.0*t2*(-1 + u) - 6.0*t*(-1 + u) + 2.0*u) -
                                       4.0*m_B3*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_v3*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_v*q2*(-1 + 2.0*u - 2.0*t*(1 + 3.0*u) + t2*(2.0 + 6.0*u)) +
                                       m_B*(-(16.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-3.0 + t*(10.0 - 24*u) + 6.0*u + 2.0*t2*(-5.0 + 12.0*u)))) +
                                       sigma2*(m_B4*(-1 + t*(18.0 - 22*u) + 2.0*u + 2.0*t2*(-9.0 + 11.0*u)) +
                                       m_B3*m_v*(1 - 2.0*u + 2.0*t*(-7.0 + 3.0*u) + t2*(14.0 - 6.0*u)) +
                                       m_B*(m_v3*(-1 + 10.0*t2*(-1 + u) - 10.0*t*(-1 + u) + 2.0*u) +
                                       m_v*q2*(1 - 2.0*u - 6.0*t2*(1 + u) + 6.0*t*(1 + u))) +
                                       m_B2*(20.0*q2*t*tbar*(-1 + u) +
                                       m_v2*(3.0 + t2*(18.0 - 34*u) - 6.0*u + 2.0*t*(-9.0 + 17.0*u))) +
                                       q2*(q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u)))) +
                                       sigma*(4.0*m_B3*m_v*t*tbar + 4.0*m_B4*t*tbar*(-1 + u) +
                                       m_v2*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       4.0*m_B*m_v*t*tbar*(-q2 + 2.0*m_v2*(-2.0 + u)) +
                                       m_v4*(1 - 2.0*u - 2.0*t2*(1 + 2.0*u) + t*(2.0 + 4.0*u)) +
                                       m_B2*(-(8.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-1 + 2.0*u + 2.0*t2*(-7.0 + 10.0*u) + t*(14.0 - 20.0*u))) +
                                       4.0*t*(-tbar + u - t*u)*pow(q2,2.0))) -
                                       4.0*u2*pow(m_B,-1)*pow(omega,2.0)*
                                       (-(2.0*m_v2*t*tbar*(m_B2*(-1 + u) + (m_v2 - q2)*(-1 + u) - m_B*m_v*(-2.0 + u))) +
                                       m_B4*sigma5*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       m_B3*sigma4*(m_B*(-3.0 + 6.0*u + t*(8.0 - 20.0*u) + 4.0*t2*(-2.0 + 5.0*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       m_B*sigma2*(2.0*m_v*t*tbar*(-q2 + m_v2*(-1 + u)) +
                                       m_B2*m_v*(-1 + 6.0*t2*(-1 + u) - 6.0*t*(-1 + u) + 2.0*u) +
                                       m_B3*(1 - 2.0*u + t2*(8.0 - 12.0*u) + 4.0*t*(-2.0 + 3.0*u)) +
                                       m_B*(q2*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)))) +
                                       sigma*(-(2.0*m_B3*m_v*t*tbar) - 2.0*m_B4*t*tbar*(-1 + u) +
                                       2.0*(m_v2 - q2)*q2*t*tbar*(-1 + u) - 2.0*m_B*m_v*t*tbar*(-q2 + m_v2*(-3.0 + 2.0*u)) +
                                       m_B2*(4.0*q2*t*(-tbar + u - t*u) +
                                       m_v2*(1 - 2.0*u + t2*(6.0 - 10.0*u) + 2.0*t*(-3.0 + 5.0*u)))) +
                                       m_B2*sigma3*(-(2.0*m_B*m_v*(1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u)) +
                                       2.0*m_v2*t*(-tbar + u - t*u) + 3.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0) -
                                       q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))));

                    const double C_3d1 = pow(m_B,-2.0)*pow(sigmabar,-6.0)*(-(2.0*m_B6*t*tbar*(-1 + u)) -
                                         2.0*m_B5*t*tbar*(m_v - 4.0*omega*(-1 + u)*u) -
                                         12.0*m_B5*sigma5*(-(m_B*t*tbar*(-1 + u)) +
                                         omega*u*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                         2.0*m_B5*sigma6*(-(m_B*t*tbar*(-1 + u)) +
                                         omega*u*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                         2.0*m_B4*t*tbar*(m_v2 + 4.0*m_v*omega*u + 2.0*(-1 + u)*(q2 - 2.0*u2*pow(omega,2.0))) -
                                         2.0*m_B3*(m_v2*omega*u*(1 - 6.0*t2 + 6.0*t - 2.0*u) - m_v3*t*tbar*(-5.0 + 3.0*u) +
                                         8.0*omega*q2*t*u*(-tbar + u - t*u) - m_v*t*tbar*(q2 - 4.0*u2*pow(omega,2.0))) +
                                         m_B2*(-(8.0*m_v*omega*q2*t*tbar*u) - 24*m_v3*omega*t*tbar*u*(-2.0 + u) -
                                         2.0*m_v4*t*tbar*(-4.0 + u) - 2.0*q2*t*tbar*(-1 + u)*(q2 - 8.0*u2*pow(omega,2.0)) +
                                         m_v2*(q2*(1 - 2.0*u - 2.0*t2*(-5.0 + 6.0*u) + 2.0*t*(-5.0 + 6.0*u)) +
                                         4.0*u2*(1 - 4.0*t2 + 4.0*t - 2.0*u)*pow(omega,2.0))) +
                                         m_B3*sigma4*(-30*m_B3*t*tbar*(-1 + u) + 8.0*omega*q2*t*u*(-tbar + u - t*u) -
                                         2.0*m_v2*omega*u*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                         2.0*m_B2*(-(m_v*t*tbar) + 2.0*omega*u*
                                         (7.0 + t2*(16.0 - 44*u) - 14.0*u + 4.0*t*(-4.0 + 11.0*u))) +
                                         2.0*m_B*(-(m_v2*t*tbar*(-1 + 2.0*u)) +
                                         2.0*m_v*omega*u*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                         q2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                         4.0*u2*(-tbar + 2.0*u - 5.0*t*u + t2*(-1 + 5.0*u))*pow(omega,2.0)) +
                                         m_v*(q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                         4.0*u2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))*pow(omega,2.0))) -
                                         8.0*u2*t*tbar*(-1 + u)*pow(omega,2.0)*(5.0*m_v4 - 6.0*m_v2*q2 + pow(q2,2.0)) +
                                         2.0*m_B*omega*u*(4.0*m_v*omega*q2*t*tbar*u + 4.0*m_v3*omega*t*tbar*u*(-7.0 + 3.0*u) +
                                         m_v2*q2*(-1 + t*(22 - 24*u) + 2.0*u + t2*(-22 + 24*u)) +
                                         m_v4*(1 - 2.0*u - 2.0*t2*(-9.0 + 7.0*u) + 2.0*t*(-9.0 + 7.0*u)) +
                                         4.0*t*(-tbar + u - t*u)*pow(q2,2.0)) -
                                         2.0*sigma*(-(6.0*m_B6*t*tbar*(-1 + u)) +
                                         2.0*m_B5*(-(2.0*m_v*t*tbar) + omega*u*
                                         (1 - 2.0*u - 2.0*t2*(-5.0 + 7.0*u) + 2.0*t*(-5.0 + 7.0*u))) -
                                         16.0*u2*(m_v2 - q2)*q2*t*tbar*(-1 + u)*pow(omega,2.0) +
                                         m_B2*(2.0*m_v4*t*tbar*(2.0 + u) +
                                         2.0*m_v3*omega*u*(1 - 2.0*u + t*(22 - 6.0*u) + t2*(-22 + 6.0*u)) +
                                         2.0*m_v*omega*q2*u*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                         q2*(q2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                         4.0*u2*(1 + (-2.0 - 4.0*t2 + 4.0*t)*u)*pow(omega,2.0)) +
                                         m_v2*(q2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)) +
                                         4.0*u2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u))*pow(omega,2.0)))
                                         + m_B3*(-(8.0*omega*q2*t*tbar*(-1 + u)*u) + 6.0*m_v3*t*tbar*(-2.0 + u) -
                                         2.0*m_v2*omega*u*(1 - 2.0*u + 2.0*t2*(-5.0 + 3.0*u) + t*(10.0 - 6.0*u)) +
                                         m_v*(q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u +
                                         4.0*u2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))*pow(omega,2.0))) -
                                         2.0*m_B*omega*u*(4.0*m_v*omega*q2*t*tbar*u -
                                         4.0*m_v3*omega*t*tbar*u*(-5.0 + 3.0*u) +
                                         m_v2*q2*(-3.0 + 6.0*u + t*(2.0 - 12.0*u) + 2.0*t2*(-1 + 6.0*u)) +
                                         m_v4*(2.0 - 4.0*u - 4.0*t2*(1 + 2.0*u) + t*(4.0 + 8.0*u)) +
                                         (1 - 2.0*u + t2*(6.0 - 10.0*u) + 2.0*t*(-3.0 + 5.0*u))*pow(q2,2.0)) +
                                         m_B4*(-(2.0*m_v2*t*tbar*(-2.0 + u)) + q2*(-1 + 4.0*t2 - 4.0*t + 2.0*u) +
                                         2.0*m_v*omega*u*(-1 + 6.0*t2*(-1 + u) - 6.0*t*(-1 + u) + 2.0*u) +
                                         4.0*u2*(-1 + 2.0*u)*pow(omega,2.0)*pow(1 - 2.0*t,2.0))) -
                                         2.0*m_B2*sigma3*(-(20.0*m_B4*t*tbar*(-1 + u)) +
                                         2.0*m_v2*omega*u*(-(4.0*omega*t*tbar*(-1 + u)*u) +
                                         m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                         4.0*m_B3*(-(m_v*t*tbar) - 4.0*omega*u*
                                         (-1 + 2.0*u + t*(3.0 - 7.0*u) + t2*(-3.0 + 7.0*u))) +
                                         2.0*m_B*(-(m_v3*t*tbar) - m_v2*omega*(1 + 6.0*t2 - 6.0*t)*u*(-1 + 2.0*u) +
                                         8.0*omega*q2*t*u*(-tbar + u - t*u) - m_v*t*tbar*(q2 - 4.0*u2*pow(omega,2.0))) +
                                         m_B2*(q2*(-3.0 + 6.0*u + t*(4.0 - 16.0*u) + 4.0*t2*(-1 + 4.0*u)) +
                                         2.0*(-(m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                         m_v*omega*u*(-3.0 + 6.0*u + t*(2.0 - 18.0*u) + 2.0*t2*(-1 + 9.0*u)) +
                                         2.0*u2*(-3.0 + 6.0*u + t*(4.0 - 16.0*u) + 4.0*t2*(-1 + 4.0*u))*pow(omega,2.0)))
                                         + (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(q2,2.0) +
                                         2.0*q2*u*(m_v2*t*tbar + m_v*omega*(1 - 2.0*u - 2.0*t2*(1 + 3.0*u) + t*(2.0 + 6.0*u)) +
                                         2.0*u*(-1 + 2.0*u)*pow(omega,2.0)*pow(1 - 2.0*t,2.0))) +
                                         3.0*m_B*sigma2*(-(10.0*m_B5*t*tbar*(-1 + u)) +
                                         2.0*m_B4*(-(2.0*m_v*t*tbar) + omega*u*
                                         (3.0 + t2*(14.0 - 26*u) - 6.0*u + 2.0*t*(-7.0 + 13.0*u))) +
                                         m_B2*(2.0*m_v3*t*tbar*(-3.0 + u) + 4.0*m_v2*omega*t*tbar*u*(-4.0 + 5.0*u) +
                                         8.0*omega*q2*t*u*(-tbar + u - t*u) +
                                         m_v*q2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                         4.0*m_v*u2*(1 + (-2.0 - 6.0*t2 + 6.0*t)*u)*pow(omega,2.0)) +
                                         2.0*m_B3*(-q2 - 2.0*m_v2*t*tbar*(-1 + u) + q2*(2.0 + 4.0*t2 - 4.0*t)*u +
                                         2.0*m_v*omega*u*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                         4.0*u2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))*pow(omega,2.0)) +
                                         2.0*omega*u*(-(4.0*m_v*omega*q2*t*tbar*u) + 4.0*m_v3*omega*t*tbar*(-1 + u)*u +
                                         m_v2*q2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u)) +
                                         (1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u))*pow(q2,2.0)) +
                                         m_B*(2.0*m_v4*t*tbar*u + 8.0*m_v*omega*q2*t*tbar*u -
                                         8.0*m_v3*omega*t*tbar*u*(-2.0 + u) +
                                         2.0*q2*t*tbar*(-1 + u)*(q2 - 8.0*u2*pow(omega,2.0)) +
                                         m_v2*(-1 + 2.0*u)*(q2*(1 + 2.0*t2 - 2.0*t) +
                                         4.0*u2*pow(omega,2.0)*pow(1 - 2.0*t,2.0)))));

                    return C_3 * phi_bar_3 + C_3d1 * phi_double_bar_3;
                }

                inline
                double I1_V23_phi_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_4 = this->phi_4(omega_1, omega_2);

                    const double C_1 = -(pow(sigmabar,-2.0)*(-(2.0*m_B2*sigma3*t*tbar*(-1 + 2.0*u)) -
                                       2.0*m_v*t*tbar*(m_v + m_B*(-3.0 + 2.0*u)) +
                                       sigma2*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u + 4.0*m_B2*t*(-tbar + 2.0*u - 2.0*t*u)) +
                                       4.0*u2*sigmabar*(m_B*sigma*(1 + 6.0*t2 - 6.0*t)*(-1 + 2.0*u) -
                                       2.0*m_v*t*tbar*(-3.0 + 2.0*u))*pow(m_B,-1)*pow(omega,2.0) +
                                       sigma*(-(2.0*t*tbar*(m_v2*(1 - 2.0*u) + m_B*m_v*(3.0 - 2.0*u) + m_B2*(-1 + 2.0*u))) +
                                       q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)) -
                                       2.0*omega*u*pow(m_B,-1)*(-(2.0*m_B2*sigma2*
                                       (-1 + 2.0*u + t*(5.0 - 12.0*u) + t2*(-5.0 + 12.0*u))) -
                                       2.0*m_v*t*tbar*(m_v + m_B*(-6.0 + 4.0*u)) +
                                       m_B2*sigma3*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)) +
                                       sigma*(m_B2*(1 + 8.0*t2 - 8.0*t)*(-1 + 2.0*u) + 4.0*m_B*m_v*t*tbar*(-3.0 + 2.0*u) +
                                       m_v2*(1 - 2.0*u + t2*(2.0 - 8.0*u) + t*(-2.0 + 8.0*u)) +
                                       q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))))/4.0;

                    return C_1 * phi_4;
                }

                inline
                double I1_V23_phi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigma2    = pow(sigma, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*pow(sigmabar,-3.0)*(-(4.0*u2*sigmabar*
                                       (-(m_B*(1 + 6.0*t2 - 6.0*t)*(-1 + 2.0*u)) + 2.0*m_v*t*tbar*(-3.0 + 2.0*u) +
                                       m_B*sigma*(-2.0 + (4.0 + 8.0*t2 - 8.0*t)*u))*pow(omega,2.0)) -
                                       m_B*sigmabar*(2.0*t*tbar*(m_B*(m_v - 2.0*m_v*sigma) + m_v2*(1 - 2.0*u) -
                                       m_B2*sigmabar*(1 - 2.0*u + sigma*(-5.0 + 6.0*u))) +
                                       q2*(sigma*(-2.0 + (4.0 + 8.0*t2 - 8.0*t)*u) - (-1 + 2.0*u)*pow(1 - 2.0*t,2.0))) +
                                       2.0*omega*u*(m_v2*(1 - 2.0*u - sigma*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) -
                                       2.0*t2*(1 + 2.0*u) + t*(2.0 + 4.0*u)) -
                                       m_B2*sigmabar*(-(2.0*sigma*(1 + 6.0*t2 - 6.0*t)*(-1 + 2.0*u)) +
                                       (1 + 8.0*t2 - 8.0*t)*(-1 + 2.0*u) +
                                       sigma2*(-5.0 + t*(10.0 - 32*u) + 10.0*u + 2.0*t2*(-5.0 + 16.0*u))) -
                                       2.0*m_B*m_v*sigmabar*(-(2.0*t*tbar*(-1 + u)) +
                                       sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       q2*(sigma*(-3.0 + 6.0*u + t*(4.0 - 16.0*u) + 4.0*t2*(-1 + 4.0*u)) -
                                       (-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))))/4.0;

                    return C_1 * phi_bar_4;
                }

                inline
                double I2_V23_phi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);

                    const double C_2 = -(pow(m_B,-2.0)*pow(sigmabar,-4.0)*(-(2.0*omega*u*
                                       (m_B4*sigma5*(-3.0 + 6.0*u + t*(6.0 - 20.0*u) + t2*(-6.0 + 20.0*u)) +
                                       2.0*m_v*(-(2.0*m_B3*t*tbar*(-3.0 + 2.0*u)) -
                                       2.0*m_B*t*tbar*(m_v2*(-1 + u) + q2*(3.0 - 2.0*u)) +
                                       m_B2*m_v*(1 + t - 2.0*u + 2.0*t*u - t2*(1 + 2.0*u)) +
                                       m_v*(m_v2 - q2)*(-tbar + 2.0*u - 6.0*t*u + t2*(-1 + 6.0*u))) -
                                       2.0*m_B3*sigma4*(m_B*(-4.0 + t*(15.0 - 36*u) + 8.0*u +
                                       3.0*t2*(-5.0 + 12.0*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       2.0*m_B2*sigma3*(3.0*m_B2*
                                       (-1 + 2.0*u + t*(9.0 - 16.0*u) + t2*(-9.0 + 16.0*u)) +
                                       q2*(-2.0 + 4.0*u + t2*(3.0 + 4.0*u) - t*(3.0 + 4.0*u)) +
                                       m_B*m_v*(3.0 - 6.0*u - 4.0*t2*(-3.0 + 5.0*u) + 4.0*t*(-3.0 + 5.0*u)) +
                                       m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) -
                                       2.0*m_B*sigma2*(-(7.0*m_B3*t*tbar*(-3.0 + 4.0*u)) +
                                       m_B2*m_v*(3.0 + t2*(24 - 26*u) - 6.0*u + t*(-24 + 26*u)) +
                                       m_v*q2*(1 - 2.0*u - 2.0*t2*(1 + 2.0*u) + t*(2.0 + 4.0*u)) +
                                       m_v3*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_B*(m_v2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)) +
                                       q2*(-2.0 + 4.0*u + t2*(9.0 - 4.0*u) + t*(-9.0 + 4.0*u)))) +
                                       sigma*(m_B4*(1 + 12.0*t2*(-1 + u) - 12.0*t*(-1 + u) - 2.0*u) +
                                       2.0*m_B3*m_v*(1 - 2.0*u - 4.0*t2*(-5.0 + 4.0*u) + 4.0*t*(-5.0 + 4.0*u)) +
                                       2.0*m_B2*(2.0*q2*t*tbar*(-3.0 + 4.0*u) +
                                       m_v2*(-3.0 + t2 - t + 6.0*u + 10.0*t2*u - 10.0*t*u)) +
                                       2.0*m_B*(m_v*q2*(1 - 8.0*t2 + 8.0*t - 2.0*u) +
                                       m_v3*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) -
                                       (m_v2 - q2)*(q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_v2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)))))) +
                                       4.0*u2*sigmabar*(-(2.0*m_B2*sigma2*
                                       (m_v*(-1 - 3.0*t2 + 3.0*t) + m_B*(1 + 6.0*t2 - 6.0*t))*(-1 + 2.0*u)) +
                                       m_B3*sigma3*(-3.0 + 6.0*u + t*(6.0 - 20.0*u) + t2*(-6.0 + 20.0*u)) +
                                       2.0*m_v*(-(m_B2*t*tbar*(-3.0 + 2.0*u)) - (m_v2 - q2)*t*tbar*(-3.0 + 2.0*u) +
                                       m_B*m_v*(1 + (-2.0 - 4.0*t2 + 4.0*t)*u)) +
                                       m_B*sigma*(m_v2*(1 + 6.0*t2 - 6.0*t)*(-1 + 2.0*u) +
                                       m_B2*(1 - 2.0*u + t*(6.0 - 4.0*u) + t2*(-6.0 + 4.0*u)) +
                                       q2*(-1 + 2.0*u + t2*(6.0 - 4.0*u) + t*(-6.0 + 4.0*u)) -
                                       2.0*m_B*m_v*(-1 + 2.0*u + t*(6.0 - 8.0*u) + t2*(-6.0 + 8.0*u))))*pow(omega,2.0)
                                       - m_B*sigmabar*(-(2.0*m_B4*sigma4*t*tbar*(-3.0 + 4.0*u)) +
                                       2.0*m_B*sigma2*(m_v*q2 + m_v*q2*(-2.0 - 6.0*t2 + 6.0*t)*u +
                                       m_B2*m_v*t*tbar*(-7.0 + 2.0*u) + m_B*q2*(-1 + 3.0*t2 - 3.0*t + 2.0*u) -
                                       3.0*m_B3*t*tbar*(-3.0 + 4.0*u) + 2.0*m_B*m_v2*t*(-tbar + 2.0*u - 2.0*t*u)) +
                                       m_B2*sigma3*(-(3.0*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       2.0*m_B*t*tbar*(2.0*m_v + 3.0*m_B*(-3.0 + 4.0*u))) -
                                       2.0*m_v*(-(m_B*m_v2*t*tbar) + m_v3*t*tbar - m_B2*m_v*t*tbar*(-1 + 2.0*u) -
                                       m_B*(m_B2 - q2)*t*tbar*(-3.0 + 2.0*u) +
                                       m_v*q2*(1 - 2.0*u + t2*(1 - 6.0*u) + t*(-1 + 6.0*u))) +
                                       sigma*(-(8.0*m_B3*m_v*t*tbar*(-2.0 + u)) +
                                       m_B2*(q2*(1 + 6.0*t2 - 6.0*t) - 6.0*m_v2*t*tbar)*(-1 + 2.0*u) +
                                       2.0*m_B4*t*tbar*(-3.0 + 4.0*u) +
                                       m_B*(-(4.0*m_v3*t*tbar) + 2.0*m_v*q2*
                                       (-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) +
                                       q2*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u + m_v2*(1 - 2.0*u*pow(1 - 2.0*t,2.0)))))))/4.0;

                    return C_2 * phi_bar_4;
                }

                inline
                double I1_V23_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to phi_double_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double sigma2    = pow(sigma, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_1 = pow(m_B,-2.0)*(m_B*sigmabar*t*tbar*(m_v + m_B*(1 + sigma - u)) +
                                       omega*u*(m_v*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)) +
                                       m_B*sigma2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)) +
                                       sigma*(m_B*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)))) +
                                       4.0*u2*(1 + sigma)*t*tbar*(-1 + u)*pow(omega,2.0))*pow(sigmabar,-3.0);

                    return C_1 * phi_double_bar_4;
                }

                inline
                double I2_V23_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_2 = -(pow(m_B,-2.0)*pow(sigmabar,-4.0)*(-(2.0*omega*u*
                                       (-(4.0*m_B3*t*tbar*(-1 + u)) + m_B2*m_v*
                                       (1 - 6.0*t2*(-1 + u) + 6.0*t*(-1 + u) - 2.0*u) +
                                       m_B3*sigma4*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_v*q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) -
                                       4.0*m_B*t*tbar*(q2 - q2*u + m_v2*(-2.0 + u)) +
                                       m_v3*(-3.0 + 6.0*u + 2.0*t2*(1 + 5.0*u) - 2.0*t*(1 + 5.0*u)) -
                                       m_B2*sigma3*(m_B*(-3.0 + 6.0*u + t*(2.0 - 12.0*u) + 2.0*t2*(-1 + 6.0*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       m_B*sigma2*(m_B*m_v*(3.0 - 6.0*u + t2*(6.0 - 14.0*u) +
                                       2.0*t*(-3.0 + 7.0*u)) +
                                       m_B2*(-3.0 + 6.0*u + t*(2.0 - 16.0*u) + 2.0*t2*(-1 + 8.0*u)) +
                                       q2*(3.0 - 6.0*u + t2*(2.0 - 16.0*u) + 2.0*t*(-1 + 8.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))) +
                                       sigma*(-(m_B3*(1 + 6.0*t2 - 6.0*t)*(-1 + 2.0*u)) +
                                       m_v3*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_v*q2*(3.0 - 6.0*u + t2*(2.0 - 14.0*u) + 2.0*t*(-1 + 7.0*u)) +
                                       m_B2*m_v*(-3.0 + 6.0*u + 2.0*t2*(-5.0 + 7.0*u) - 2.0*t*(-5.0 + 7.0*u)) +
                                       m_B*(-(q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)) +
                                       m_v2*(-3.0 + 6.0*u + 2.0*t2*(1 + 6.0*u) - 2.0*t*(1 + 6.0*u)))))) +
                                       4.0*u2*(-(2.0*m_B2*t*tbar*(-1 + u)) - 2.0*(3.0*m_v2 - q2)*t*tbar*(-1 + u) -
                                       m_B*m_v*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_B*sigma2*(10.0*m_B*t*tbar*(-1 + u) + m_v*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)) +
                                       m_B2*sigma3*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)) +
                                       sigma*(-(2.0*(m_v2 - 3.0*q2)*t*tbar*(-1 + u)) +
                                       m_B2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u))))*pow(omega,2.0) -
                                       m_B*sigmabar*(2.0*m_v3*t*tbar + 2.0*m_B*(m_B2 - q2)*t*tbar*(-1 + u) -
                                       2.0*m_B*m_v2*t*tbar*(1 + u) - 2.0*m_B2*m_v*t*tbar*(-2.0 + u) -
                                       2.0*m_B3*sigma3*t*tbar*(-1 + 2.0*u) +
                                       m_v*q2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)) +
                                       sigma*(2.0*m_B2*m_v*t*tbar*(-3.0 + u) - 2.0*m_B3*t*tbar*(-3.0 + 4.0*u) +
                                       m_v*q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       m_B*(2.0*m_v2*t*tbar*(1 + 2.0*u) + q2*(-1 + 2.0*t2 - 2.0*t + 2.0*u))) -
                                       m_B*sigma2*(-(2.0*m_B*t*tbar*(m_v + m_B*(-3.0 + 5.0*u))) +
                                       q2*(1 - 2.0*u*pow(1 - 2.0*t,2.0))))))/2.0;

                    return C_2 * phi_double_bar_4;
                }

                inline
                double I3_V23_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4), m_B5 = pow(m_B, 5);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5), sigma6 = pow(sigma, 6);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_3 = pow(m_B,-2.0)*pow(sigmabar,-5.0)*(-(m_B*sigmabar*
                                       (-(2.0*m_B5*sigma5*t*tbar*(-1 + u)) +
                                       m_B3*sigma4*(q2 + 2.0*m_B*t*tbar*(m_v + 4.0*m_B*(-1 + u)) +
                                       q2*(-2.0 - 4.0*t2 + 4.0*t)*u) +
                                       m_v2*(2.0*m_B*m_v2*t*tbar - 2.0*m_B*(m_B2 - q2)*t*tbar*(-1 + u) +
                                       m_v*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u + 2.0*m_B2*t*tbar*(-2.0 + u))) +
                                       m_B2*sigma3*(-(6.0*m_B2*m_v*t*tbar) - 12.0*m_B3*t*tbar*(-1 + u) +
                                       m_v*q2*(1 + (-2.0 - 6.0*t2 + 6.0*t)*u) +
                                       2.0*m_B*(m_v2*t*(-tbar + 2.0*u - 2.0*t*u) +
                                       q2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)))) +
                                       m_B*sigma2*(6.0*m_B3*m_v*t*tbar + 8.0*m_B4*t*tbar*(-1 + u) +
                                       2.0*m_B*(-(m_v3*t*tbar) + m_v*q2*(-1 + t2 - t + 2.0*u + 5.0*t2*u - 5.0*t*u)) +
                                       q2*(2.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       m_B2*(-(2.0*m_v2*t*tbar*(-3.0 + 5.0*u)) +
                                       q2*(1 - 2.0*u + 4.0*t2*(-2.0 + u) - 4.0*t*(-2.0 + u)))) +
                                       sigma*(-(2.0*m_B4*m_v*t*tbar) - 2.0*m_B5*t*tbar*(-1 + u) +
                                       2.0*m_B3*t*tbar*(2.0*q2*(-1 + u) + m_v2*(-3.0 + 4.0*u)) +
                                       m_B2*m_v*(-(2.0*m_v2*t*tbar*(-3.0 + u)) +
                                       q2*(1 - 2.0*u - 2.0*t2*(1 + 2.0*u) + t*(2.0 + 4.0*u))) +
                                       m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(q2,2.0) +
                                       m_B*(-(2.0*m_v4*t*tbar*u) - m_v2*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) -
                                       2.0*t*tbar*(-1 + u)*pow(q2,2.0))))) -
                                       2.0*omega*u*(m_B5*sigma6*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_v2*(4.0*m_B3*t*tbar*(-1 + u) +
                                       m_B2*m_v*(-1 + 8.0*t2*(-1 + u) - 8.0*t*(-1 + u) + 2.0*u) -
                                       m_v*(m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       2.0*m_B*t*tbar*(-(2.0*q2*(-1 + u)) + m_v2*(-2.0 + u))) +
                                       m_B3*sigma4*(4.0*q2*t*(-tbar + u - t*u) -
                                       2.0*m_B*m_v*(-2.0 + 4.0*u + t*(5.0 - 11.0*u) + t2*(-5.0 + 11.0*u)) +
                                       2.0*m_B2*(-3.0 + t*(14.0 - 26*u) + 6.0*u + 2.0*t2*(-7.0 + 13.0*u)) +
                                       m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) -
                                       m_B4*sigma5*(m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       4.0*m_B*(-1 + 2.0*u + t*(3.0 - 7.0*u) + t2*(-3.0 + 7.0*u))) +
                                       m_B2*sigma3*(-(2.0*m_v*q2*t*tbar*(1 + u)) -
                                       4.0*m_B3*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       6.0*m_B2*m_v*(-1 + 2.0*u + t*(3.0 - 5.0*u) + t2*(-3.0 + 5.0*u)) +
                                       m_v3*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_B*(-(16.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-3.0 + t*(10.0 - 24*u) + 6.0*u + 2.0*t2*(-5.0 + 12.0*u)))) +
                                       m_B*sigma2*(m_B4*(-1 + t*(18.0 - 22*u) + 2.0*u + 2.0*t2*(-9.0 + 11.0*u)) -
                                       2.0*m_B3*m_v*(-2.0 + 4.0*u + t*(7.0 - 9.0*u) + t2*(-7.0 + 9.0*u)) +
                                       m_B2*(20.0*q2*t*tbar*(-1 + u) +
                                       m_v2*(3.0 + t2*(18.0 - 34*u) - 6.0*u + 2.0*t*(-9.0 + 17.0*u))) +
                                       q2*(q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) +
                                       m_B*(2.0*m_v*q2*t*tbar*(3.0 + u) +
                                       m_v3*(-3.0 + 6.0*u + 2.0*t2*(-5.0 + 9.0*u) - 2.0*t*(-5.0 + 9.0*u)))) -
                                       sigma*(-(4.0*m_B5*t*tbar*(-1 + u)) +
                                       m_B4*m_v*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u) -
                                       m_v*(m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_B3*(8.0*q2*t*(-tbar + u - t*u) +
                                       m_v2*(1 - 2.0*u + 2.0*t*(-7.0 + 10.0*u) + t2*(14.0 - 20.0*u))) +
                                       m_B2*(4.0*m_v*q2*t*tbar + m_v3*
                                       (-3.0 + 6.0*u + 4.0*t2*(-4.0 + 5.0*u) - 4.0*t*(-4.0 + 5.0*u))) +
                                       m_B*(-(m_v2*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)) +
                                       m_v4*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u)) -
                                       4.0*t*tbar*(-1 + u)*pow(q2,2.0)))) -
                                       4.0*u2*pow(omega,2.0)*(m_B4*sigma5*
                                       (-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       m_B3*sigma4*(m_B*(-3.0 + 6.0*u + t*(8.0 - 20.0*u) + 4.0*t2*(-2.0 + 5.0*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       m_v2*(-(2.0*m_B2*t*tbar*(-1 + u)) - 2.0*(m_v2 - q2)*t*tbar*(-1 + u) +
                                       m_B*m_v*(1 - 2.0*u + t2*(4.0 - 6.0*u) + t*(-4.0 + 6.0*u))) +
                                       m_B*sigma2*(2.0*m_v3*t*tbar*(-1 + u) - m_v*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_B3*(1 - 2.0*u + t2*(8.0 - 12.0*u) + 4.0*t*(-2.0 + 3.0*u)) +
                                       m_B2*m_v*(-3.0 + 6.0*u + t*(6.0 - 14.0*u) + 2.0*t2*(-3.0 + 7.0*u)) +
                                       m_B*(q2*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)))) +
                                       sigma*(-(2.0*m_B4*t*tbar*(-1 + u)) + 2.0*(m_v2 - q2)*q2*t*tbar*(-1 + u) -
                                       m_B3*m_v*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_B2*(4.0*q2*t*(-tbar + u - t*u) +
                                       m_v2*(1 - 2.0*u + t2*(6.0 - 10.0*u) + 2.0*t*(-3.0 + 5.0*u))) +
                                       m_B*m_v*(q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_v2*(-1 + 2.0*u + t*(6.0 - 8.0*u) + t2*(-6.0 + 8.0*u)))) +
                                       m_B2*sigma3*(2.0*m_v2*t*(-tbar + u - t*u) +
                                       m_B*m_v*(3.0 - 6.0*u + t2*(6.0 - 16.0*u) + 2.0*t*(-3.0 + 8.0*u)) +
                                       3.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0) - q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))));

                    return C_3 * phi_double_bar_4;
                }

                inline
                double I3d1_V23_phi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to phi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4), m_B5 = pow(m_B, 5), m_B6 = pow(m_B, 6);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5), sigma6 = pow(sigma, 6);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);
                    const double u3        = pow(u, 3);

                    const double phi_bar_4 = this->phi_bar_4(omega_1, omega_2);
                    const double phi_double_bar_4 = this->phi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(pow(m_B,-1)*pow(sigmabar,-5.0)*(m_B*(1 - sigma)*
                                       (-(2.0*m_B5*sigma5*t*tbar*(-1 + u)) +
                                       m_B3*sigma4*(q2 + 2.0*m_B*t*tbar*(m_v + 4.0*m_B*(-1 + u)) +
                                       q2*(-2.0 - 4.0*t2 + 4.0*t)*u) +
                                       m_v2*(2.0*m_B*m_v2*t*tbar - 2.0*m_B*(m_B2 - q2)*t*tbar*(-1 + u) +
                                       m_v*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u + 2.0*m_B2*t*tbar*(-2.0 + u))) +
                                       m_B2*sigma3*(-(6.0*m_B2*m_v*t*tbar) - 12.0*m_B3*t*tbar*(-1 + u) +
                                       m_v*q2*(1 + (-2.0 - 6.0*t2 + 6.0*t)*u) +
                                       2.0*m_B*(m_v2*t*(-tbar + 2.0*u - 2.0*t*u) +
                                       q2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)))) +
                                       m_B*sigma2*(6.0*m_B3*m_v*t*tbar + 8.0*m_B4*t*tbar*(-1 + u) +
                                       2.0*m_B*(-(m_v3*t*tbar) + m_v*q2*(-1 + t2 - t + 2.0*u + 5.0*t2*u - 5.0*t*u)) +
                                       q2*(2.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       m_B2*(-(2.0*m_v2*t*tbar*(-3.0 + 5.0*u)) +
                                       q2*(1 - 2.0*u + 4.0*t2*(-2.0 + u) - 4.0*t*(-2.0 + u)))) +
                                       sigma*(-(2.0*m_B4*m_v*t*tbar) - 2.0*m_B5*t*tbar*(-1 + u) +
                                       2.0*m_B3*t*tbar*(2.0*q2*(-1 + u) + m_v2*(-3.0 + 4.0*u)) +
                                       m_B2*m_v*(-(2.0*m_v2*t*tbar*(-3.0 + u)) +
                                       q2*(1 - 2.0*u - 2.0*t2*(1 + 2.0*u) + t*(2.0 + 4.0*u))) +
                                       m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(q2,2.0) +
                                       m_B*(-(2.0*m_v4*t*tbar*u) - m_v2*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) -
                                       2.0*t*tbar*(-1 + u)*pow(q2,2.0)))) +
                                       2.0*omega*u*(m_B5*sigma6*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_v2*(4.0*m_B3*t*tbar*(-1 + u) +
                                       m_B2*m_v*(-1 + 8.0*t2*(-1 + u) - 8.0*t*(-1 + u) + 2.0*u) -
                                       m_v*(m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       2.0*m_B*t*tbar*(-(2.0*q2*(-1 + u)) + m_v2*(-2.0 + u))) +
                                       m_B3*sigma4*(4.0*q2*t*(-tbar + u - t*u) -
                                       2.0*m_B*m_v*(-2.0 + 4.0*u + t*(5.0 - 11.0*u) + t2*(-5.0 + 11.0*u)) +
                                       2.0*m_B2*(-3.0 + t*(14.0 - 26*u) + 6.0*u + 2.0*t2*(-7.0 + 13.0*u)) +
                                       m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) -
                                       m_B4*sigma5*(m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       4.0*m_B*(-1 + 2.0*u + t*(3.0 - 7.0*u) + t2*(-3.0 + 7.0*u))) +
                                       m_B2*sigma3*(-(2.0*m_v*q2*t*tbar*(1 + u)) -
                                       4.0*m_B3*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       6.0*m_B2*m_v*(-1 + 2.0*u + t*(3.0 - 5.0*u) + t2*(-3.0 + 5.0*u)) +
                                       m_v3*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_B*(-(16.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-3.0 + t*(10.0 - 24*u) + 6.0*u + 2.0*t2*(-5.0 + 12.0*u)))) +
                                       m_B*sigma2*(m_B4*(-1 + t*(18.0 - 22*u) + 2.0*u + 2.0*t2*(-9.0 + 11.0*u)) -
                                       2.0*m_B3*m_v*(-2.0 + 4.0*u + t*(7.0 - 9.0*u) + t2*(-7.0 + 9.0*u)) +
                                       m_B2*(20.0*q2*t*tbar*(-1 + u) +
                                       m_v2*(3.0 + t2*(18.0 - 34*u) - 6.0*u + 2.0*t*(-9.0 + 17.0*u))) +
                                       q2*(q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) +
                                       m_B*(2.0*m_v*q2*t*tbar*(3.0 + u) +
                                       m_v3*(-3.0 + 6.0*u + 2.0*t2*(-5.0 + 9.0*u) - 2.0*t*(-5.0 + 9.0*u)))) -
                                       sigma*(-(4.0*m_B5*t*tbar*(-1 + u)) +
                                       m_B4*m_v*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u) -
                                       m_v*(m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_B3*(8.0*q2*t*(-tbar + u - t*u) +
                                       m_v2*(1 - 2.0*u + 2.0*t*(-7.0 + 10.0*u) + t2*(14.0 - 20.0*u))) +
                                       m_B2*(4.0*m_v*q2*t*tbar + m_v3*
                                       (-3.0 + 6.0*u + 4.0*t2*(-4.0 + 5.0*u) - 4.0*t*(-4.0 + 5.0*u))) +
                                       m_B*(-(m_v2*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)) +
                                       m_v4*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u)) -
                                       4.0*t*tbar*(-1 + u)*pow(q2,2.0)))) +
                                       4.0*u2*pow(omega,2.0)*(m_B4*sigma5*
                                       (-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       m_B3*sigma4*(m_B*(-3.0 + 6.0*u + t*(8.0 - 20.0*u) + 4.0*t2*(-2.0 + 5.0*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       m_v2*(-(2.0*m_B2*t*tbar*(-1 + u)) - 2.0*(m_v2 - q2)*t*tbar*(-1 + u) +
                                       m_B*m_v*(1 - 2.0*u + t2*(4.0 - 6.0*u) + t*(-4.0 + 6.0*u))) +
                                       m_B*sigma2*(2.0*m_v3*t*tbar*(-1 + u) - m_v*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_B3*(1 - 2.0*u + t2*(8.0 - 12.0*u) + 4.0*t*(-2.0 + 3.0*u)) +
                                       m_B2*m_v*(-3.0 + 6.0*u + t*(6.0 - 14.0*u) + 2.0*t2*(-3.0 + 7.0*u)) +
                                       m_B*(q2*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)))) +
                                       sigma*(-(2.0*m_B4*t*tbar*(-1 + u)) + 2.0*(m_v2 - q2)*q2*t*tbar*(-1 + u) -
                                       m_B3*m_v*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_B2*(4.0*q2*t*(-tbar + u - t*u) +
                                       m_v2*(1 - 2.0*u + t2*(6.0 - 10.0*u) + 2.0*t*(-3.0 + 5.0*u))) +
                                       m_B*m_v*(q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_v2*(-1 + 2.0*u + t*(6.0 - 8.0*u) + t2*(-6.0 + 8.0*u)))) +
                                       m_B2*sigma3*(2.0*m_v2*t*(-tbar + u - t*u) +
                                       m_B*m_v*(3.0 - 6.0*u + t2*(6.0 - 16.0*u) + 2.0*t*(-3.0 + 8.0*u)) +
                                       3.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0) - q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))));

                    const double C_3d1 = pow(m_B,-2.0)*pow(sigmabar,-6.0)*(2.0*m_B6*t*tbar*(-1 + u) +
                                         2.0*m_B5*t*tbar*(m_v - 4.0*omega*(-1 + u)*u) +
                                         12.0*m_B5*sigma5*(-(m_B*t*tbar*(-1 + u)) +
                                         omega*u*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) -
                                         2.0*m_B5*sigma6*(-(m_B*t*tbar*(-1 + u)) +
                                         omega*u*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                         m_B3*(16.0*omega*q2*t*tbar*(-1 + u)*u - 2.0*m_v2*omega*u*(-1 + 6.0*t2 - 6.0*t + 2.0*u) -
                                         2.0*m_v3*t*tbar*(-5.0 + 3.0*u) +
                                         m_v*q2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u)) +
                                         4.0*m_v*u2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)*pow(omega,2.0)) +
                                         2.0*m_B4*(-(m_v2*t*tbar) + m_v*omega*u*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u) -
                                         2.0*t*tbar*(-1 + u)*(q2 - 2.0*u2*pow(omega,2.0))) +
                                         m_B2*(8.0*m_v*omega*q2*t*tbar*u + 2.0*m_v4*t*tbar*(-4.0 + u) +
                                         4.0*m_v3*omega*u*(1 - 2.0*u - 2.0*t2*(-6.0 + 5.0*u) + 2.0*t*(-6.0 + 5.0*u)) +
                                         2.0*q2*t*tbar*(-1 + u)*(q2 - 8.0*u2*pow(omega,2.0)) +
                                         m_v2*(q2*(-1 + 2.0*u + 2.0*t2*(-5.0 + 6.0*u) - 2.0*t*(-5.0 + 6.0*u)) +
                                         4.0*u2*(-1 + 4.0*t2 - 4.0*t + 2.0*u)*pow(omega,2.0))) +
                                         2.0*omega*u*(4.0*omega*t*u*(-tbar + u - t*u) + m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))*
                                         (5.0*m_v4 - 6.0*m_v2*q2 + pow(q2,2.0)) +
                                         m_B*(2.0*m_v2*omega*q2*u*(1 + t2*(22 - 24*u) - 2.0*u + t*(-22 + 24*u)) +
                                         2.0*m_v4*omega*u*(-1 + 2.0*u + 2.0*t2*(-9.0 + 7.0*u) - 2.0*t*(-9.0 + 7.0*u)) +
                                         m_v*q2*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u -
                                         4.0*u2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)*pow(omega,2.0)) +
                                         4.0*m_v3*(q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                         2.0*u2*(-2.0 + 4.0*u + t*(7.0 - 11.0*u) + t2*(-7.0 + 11.0*u))*pow(omega,2.0)) -
                                         8.0*omega*t*tbar*(-1 + u)*u*pow(q2,2.0)) +
                                         2.0*sigma*(-(6.0*m_B6*t*tbar*(-1 + u)) -
                                         4.0*omega*(m_v2 - q2)*q2*u*(4.0*omega*t*u*(-tbar + u - t*u) +
                                         m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                         2.0*m_B5*(-(2.0*m_v*t*tbar) + omega*u*
                                         (1 - 2.0*u - 2.0*t2*(-5.0 + 7.0*u) + 2.0*t*(-5.0 + 7.0*u))) +
                                         m_B2*(-(4.0*m_v*omega*q2*t*tbar*(-1 + u)*u) + 2.0*m_v4*t*tbar*(2.0 + u) +
                                         2.0*m_v3*omega*u*(-3.0 + 22*t2*(-1 + u) - 22*t*(-1 + u) + 6.0*u) +
                                         q2*(q2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                         4.0*u2*(1 + (-2.0 - 4.0*t2 + 4.0*t)*u)*pow(omega,2.0)) +
                                         m_v2*(q2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)) +
                                         4.0*u2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u))*pow(omega,2.0)))
                                         + m_B3*(-(8.0*omega*q2*t*tbar*(-1 + u)*u) + 6.0*m_v3*t*tbar*(-2.0 + u) -
                                         2.0*m_v2*omega*u*(1 - 2.0*u + 2.0*t2*(-5.0 + 3.0*u) + t*(10.0 - 6.0*u)) +
                                         m_v*(q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u +
                                         4.0*u2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))*pow(omega,2.0))) +
                                         m_B*(4.0*m_v4*omega*u*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u)) +
                                         2.0*m_v2*omega*q2*u*(3.0 - 6.0*u + t2*(2.0 - 12.0*u) + 2.0*t*(-1 + 6.0*u)) +
                                         m_v*q2*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u -
                                         4.0*u2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)*pow(omega,2.0)) +
                                         2.0*m_v3*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u -
                                         4.0*u2*(-1 + 2.0*u + t*(5.0 - 7.0*u) + t2*(-5.0 + 7.0*u))*pow(omega,2.0)) +
                                         2.0*omega*u*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u))*pow(q2,2.0)) +
                                         m_B4*(-(2.0*m_v2*t*tbar*(-2.0 + u)) + q2*(-1 + 4.0*t2 - 4.0*t + 2.0*u) +
                                         4.0*m_v*omega*u*(-1 + 2.0*u + t*(3.0 - 5.0*u) + t2*(-3.0 + 5.0*u)) +
                                         4.0*u2*(-1 + 2.0*u)*pow(omega,2.0)*pow(1 - 2.0*t,2.0))) +
                                         2.0*m_B2*sigma3*(-(20.0*m_B4*t*tbar*(-1 + u)) +
                                         2.0*m_v2*omega*u*(-(4.0*omega*t*tbar*(-1 + u)*u) +
                                         m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                         4.0*m_B3*(-(m_v*t*tbar) - 4.0*omega*u*
                                         (-1 + 2.0*u + t*(3.0 - 7.0*u) + t2*(-3.0 + 7.0*u))) +
                                         m_B*(-(2.0*m_v3*t*tbar) + 16.0*omega*q2*t*tbar*(-1 + u)*u -
                                         2.0*m_v2*omega*(1 + 6.0*t2 - 6.0*t)*u*(-1 + 2.0*u) +
                                         m_v*q2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u)) +
                                         4.0*m_v*u2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)*pow(omega,2.0)) +
                                         m_B2*(-(2.0*m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                         q2*(-3.0 + 6.0*u + t*(4.0 - 16.0*u) + 4.0*t2*(-1 + 4.0*u)) +
                                         4.0*m_v*omega*u*(-tbar + 2.0*u - 7.0*t*u + t2*(-1 + 7.0*u)) +
                                         4.0*u2*(-3.0 + 6.0*u + t*(4.0 - 16.0*u) + 4.0*t2*(-1 + 4.0*u))*pow(omega,2.0)) +
                                         (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(q2,2.0) +
                                         2.0*q2*u*(m_v2*t*tbar + 2.0*m_v*omega*t*tbar*(1 + u) +
                                         2.0*u*(-1 + 2.0*u)*pow(omega,2.0)*pow(1 - 2.0*t,2.0))) +
                                         m_B3*sigma4*(30*m_B3*t*tbar*(-1 + u) - 8.0*omega*q2*t*tbar*(-1 + u)*u +
                                         2.0*m_v2*omega*u*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                         2.0*m_B2*(-(m_v*t*tbar) + 2.0*omega*u*
                                         (7.0 + t2*(16.0 - 44*u) - 14.0*u + 4.0*t*(-4.0 + 11.0*u))) +
                                         m_v*(q2 + q2*(-2.0 - 6.0*t2 + 6.0*t)*u +
                                         4.0*u2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))*pow(omega,2.0)) -
                                         2.0*m_B*(-(m_v2*t*tbar*(-1 + 2.0*u)) +
                                         q2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                         4.0*u2*(-tbar + 2.0*u - 5.0*t*u + t2*(-1 + 5.0*u))*pow(omega,2.0) +
                                         m_v*omega*u*(-1 + 2.0*u*pow(1 - 2.0*t,2.0)))) -
                                         3.0*m_B*sigma2*(-(10.0*m_B5*t*tbar*(-1 + u)) +
                                         2.0*m_v2*omega*q2*u*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u)) +
                                         2.0*m_B4*(-(2.0*m_v*t*tbar) + omega*u*
                                         (3.0 + t2*(14.0 - 26*u) - 6.0*u + 2.0*t*(-7.0 + 13.0*u))) +
                                         8.0*m_v3*u2*t*tbar*(-1 + u)*pow(omega,2.0) +
                                         2.0*m_B2*t*tbar*(m_v*q2*(-1 + u) + 4.0*omega*q2*(-1 + u)*u + m_v3*(-3.0 + u) +
                                         2.0*m_v2*omega*u*(-4.0 + 5.0*u) + 4.0*m_v*u3*pow(omega,2.0)) +
                                         m_v*q2*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u -
                                         4.0*u2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)*pow(omega,2.0)) +
                                         2.0*m_B3*(-q2 - 2.0*m_v2*t*tbar*(-1 + u) + q2*(2.0 + 4.0*t2 - 4.0*t)*u +
                                         2.0*m_v*omega*u*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                         4.0*u2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))*pow(omega,2.0)) +
                                         2.0*omega*u*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u))*pow(q2,2.0) +
                                         m_B*(2.0*m_v4*t*tbar*u + 8.0*m_v*omega*q2*t*tbar*u +
                                         4.0*m_v3*omega*u*(-1 + 2.0*u + t*(4.0 - 6.0*u) + t2*(-4.0 + 6.0*u)) +
                                         2.0*q2*t*tbar*(-1 + u)*(q2 - 8.0*u2*pow(omega,2.0)) +
                                         m_v2*(-1 + 2.0*u)*(q2*(1 + 2.0*t2 - 2.0*t) +
                                         4.0*u2*pow(omega,2.0)*pow(1 - 2.0*t,2.0)))));

                    return C_3 * phi_bar_4 + C_3d1 * phi_double_bar_4;
                }

                inline
                double I1_V23_psi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

                    const double C_1 = (pow(m_B,-2.0)*(-(m_B*sigmabar*(-(2.0*m_B2*sigma2*t*tbar) -
                                       2.0*t*tbar*(-m_v2 + q2 + m_B2*(-1 + u) - q2*u) +
                                       sigma*(2.0*m_B2*t*tbar*u + q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u)))) +
                                       2.0*omega*u*(-(2.0*(m_B2 - q2)*t*tbar*(-1 + u)) +
                                       m_v2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_B2*sigma3*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_B2*sigma2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       sigma*(m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) -
                                       2.0*(-(m_B2*t*tbar*(-1 + u)) + q2*(1 + t2 - t - 2.0*u - 5.0*t2*u + 5.0*t*u)))) -
                                       4.0*m_B*u2*sigma*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(omega,2.0))*
                                       pow(sigmabar,-3.0))/2.0;

                    return C_1 * psi_bar_4;
                }

                inline
                double I2_V23_psi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double sigma2    = pow(sigma, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);

                    const double C_2 = -((-m_v2 + m_B2*sigma2 + (-m_B2 + q2)*sigma)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_B,-2.0)*
                                       (2.0*m_B2*omega*sigmabar2*u + 2.0*omega*(-m_v2 + q2)*u -
                                       m_B*sigmabar*(q2 + 4.0*u2*pow(omega,2.0)))*pow(sigmabar,-4.0))/2.0;

                    return C_2 * psi_bar_4;
                }

                inline
                double I1_V23_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to psi_double_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double t2        = pow(t, 2);

                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_1 = -(m_v*(m_B*sigmabar*t*tbar + omega*u*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u) +
                                       sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))))*pow(m_B,-2.0)*
                                       pow(sigmabar,-3.0));

                    return C_1 * psi_double_bar_4;
                }

                inline
                double I2_V23_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_double_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigma2    = pow(sigma, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_2 = -(m_v*pow(m_B,-2.0)*(-(m_B*sigmabar*(2.0*t*tbar*(-m_v2 - m_B2*sigmabar*(1 + sigma - u)) +
                                       q2*(1 + sigma - 2.0*u + sigma*(-2.0 - 6.0*t2 + 6.0*t)*u - 2.0*t2*(1 + u) +
                                       2.0*t*(1 + u)))) + 2.0*omega*u*
                                       (m_B2*sigmabar*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u) +
                                       sigma*(-2.0 + (4.0 + 8.0*t2 - 8.0*t)*u) +
                                       sigma2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       m_v2*(-3.0 + 6.0*u + 2.0*t2*(1 + 5.0*u) - 2.0*t*(1 + 5.0*u) +
                                       sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u) +
                                       sigma*(3.0 - 6.0*u + t2*(2.0 - 14.0*u) + 2.0*t*(-1 + 7.0*u)))) +
                                       4.0*m_B*u2*(1 + sigma)*sigmabar*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(omega,2.0))*
                                       pow(sigmabar,-4.0))/2.0;

                    return C_2 * psi_double_bar_4;
                }

                inline
                double I3_V23_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigmabar2 = pow(sigmabar, 2);
                    const double sigma2    = pow(sigma, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(m_v*(-m_v2 + m_B2*sigma2 + (-m_B2 + q2)*sigma)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*
                                       pow(m_B,-2.0)*(2.0*m_B2*omega*sigmabar2*u + 2.0*omega*(-m_v2 + q2)*u -
                                       m_B*sigmabar*(q2 + 4.0*u2*pow(omega,2.0)))*pow(sigmabar,-5.0));

                    return C_3 * psi_double_bar_4;
                }

                inline
                double I3d1_V23_psi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to psi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigmabar2 = pow(sigmabar, 2), sigmabar3 = pow(sigmabar, 3);
                    const double sigma2    = pow(sigma, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double psi_bar_4 = this->psi_bar_4(omega_1, omega_2);
                    const double psi_double_bar_4 = this->psi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(m_v*(-m_v2 + m_B2*sigma2 + (-m_B2 + q2)*sigma)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_B,-1)*
                                       (2.0*m_B2*omega*sigmabar2*u + 2.0*omega*(-m_v2 + q2)*u -
                                       m_B*sigmabar*(q2 + 4.0*u2*pow(omega,2.0)))*pow(sigmabar,-5.0));

                    const double C_3d1 = -(m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(m_B,-2.0)*
                                         (-(2.0*m_B4*omega*sigmabar3*(1 + sigma)*u) +
                                         2.0*omega*(m_v2 - q2)*(5.0*m_v2 - q2*(1 + 4.0*sigma))*u -
                                         4.0*m_B2*omega*(q2*sigma*(1 + sigma) + m_v2*(1 - 3.0*sigma))*sigmabar*u +
                                         m_B3*sigmabar2*(1 + 2.0*sigma)*(q2 + 4.0*u2*pow(omega,2.0)) -
                                         m_B*(-(4.0*m_v2) + q2 + 3.0*q2*sigma)*sigmabar*(q2 + 4.0*u2*pow(omega,2.0)))*
                                         pow(sigmabar,-6.0));

                    return C_3 * psi_bar_4 + C_3d1 * psi_double_bar_4;
                }

                inline
                double I1_V23_chi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2);
                    const double m_B2      = pow(m_B, 2);
                    const double sigma2    = pow(sigma, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*(m_B*(q2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u) +
                                       sigma*(-1 + 2.0*(1 + t2 - t)*u)) +
                                       2.0*t*tbar*(2.0*m_B*m_v + m_B2*u - 2.0*m_v2*u + m_B2*sigma2*(-3.0 + 4.0*u) +
                                       m_B*sigma*(-(2.0*m_v) + m_B*(3.0 - 5.0*u)))) +
                                       2.0*omega*u*(-(2.0*m_v2*t*tbar*u) +
                                       q2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       2.0*m_B*m_v*(-(t*tbar) + sigma*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                       m_B2*(-1 + 2.0*u + t*(2.0 - 10.0*u) + 2.0*t2*(-1 + 5.0*u) +
                                       3.0*sigma2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       2.0*sigma*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)))) +
                                       4.0*u2*(-(2.0*m_v*t*tbar) + m_B*
                                       (1 - 2.0*u + t2*(2.0 - 8.0*u) + t*(-2.0 + 8.0*u) +
                                       sigma*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))))*pow(omega,2.0))*
                                       pow(sigmabar,-2.0))/2.0;

                    return C_1 * chi_bar_4;
                }

                inline
                double I2_V23_chi_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);

                    const double C_2 = -(pow(m_B,-2.0)*pow(sigmabar,-3.0)*(-(4.0*u2*pow(omega,2.0)*
                                       (2.0*m_B3*sigma3*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_v*(2.0*m_B2*t*tbar + 2.0*(m_v2 - q2)*t*tbar +
                                       m_B*m_v*(1 - 2.0*u - 2.0*t2*(1 + 2.0*u) + t*(2.0 + 4.0*u))) +
                                       m_B*sigma*(-(4.0*m_B2*t*tbar*(-1 + u)) -
                                       2.0*m_B*m_v*(1 + 3.0*t2 - 3.0*t)*(-1 + 2.0*u) + 4.0*q2*t*(-tbar + u - t*u) +
                                       m_v2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))) -
                                       2.0*m_B2*sigma2*(m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_B*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))) -
                                       2.0*omega*u*(2.0*m_B4*sigma4*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_B2*sigma2*(-(m_v2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u)) +
                                       2.0*q2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)) -
                                       2.0*m_B*m_v*(-2.0 + 4.0*u + t*(7.0 - 12.0*u) + t2*(-7.0 + 12.0*u)) +
                                       m_B2*(-2.0 + 4.0*u + 4.0*t2*(-5.0 + 7.0*u) - 4.0*t*(-5.0 + 7.0*u))) +
                                       m_v*(-(4.0*m_B3*t*tbar) + 2.0*m_B*(m_v2 + 2.0*q2)*t*tbar +
                                       m_B2*m_v*(-1 + 6.0*t2 - 6.0*t + 2.0*u) -
                                       m_v*(m_v2 - q2)*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))) -
                                       2.0*m_B*sigma*(-(4.0*m_B3*t*tbar*(-1 + u)) +
                                       m_v*q2*(1 + t - 2.0*u + 6.0*t*u - t2*(1 + 6.0*u)) +
                                       m_v3*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_B2*m_v*(1 - 2.0*u + t2*(7.0 - 6.0*u) + t*(-7.0 + 6.0*u)) +
                                       m_B*(4.0*q2*t*(-tbar + u - t*u) +
                                       m_v2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u)))) -
                                       2.0*m_B3*sigma3*(m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       2.0*m_B*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0))) +
                                       m_B*(-(4.0*m_B4*sigma4*t*tbar*(-1 + u)) -
                                       2.0*m_B2*sigma3*(-(2.0*m_B*t*tbar*(m_v + 3.0*m_B*(-1 + u))) +
                                       q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       2.0*m_B*sigma2*(m_v*q2 - 5.0*m_B2*m_v*t*tbar - 6.0*m_B3*t*tbar*(-1 + u) +
                                       m_v*q2*(-2.0 - 6.0*t2 + 6.0*t)*u + m_B*m_v2*t*(-tbar + 2.0*u - 2.0*t*u) +
                                       m_B*q2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u))) +
                                       m_v*(4.0*m_B*m_v2*t*tbar - 2.0*m_v3*t*tbar - 2.0*m_B*(m_B2 - q2)*t*tbar +
                                       4.0*m_B2*m_v*t*(-tbar + u - t*u) +
                                       m_v*q2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u))) +
                                       sigma*(8.0*m_B3*m_v*t*tbar + 4.0*m_B4*t*tbar*(-1 + u) -
                                       2.0*m_B2*t*tbar*(2.0*q2*(-1 + u) + m_v2*(-3.0 + 4.0*u)) +
                                       2.0*m_B*(-(2.0*m_v3*t*tbar) +
                                       m_v*q2*(-1 + t2 - t + 2.0*u + 6.0*t2*u - 6.0*t*u)) +
                                       m_v2*q2*(1 - 2.0*u*pow(1 - 2.0*t,2.0))))))/2.0;

                    return C_2 * chi_bar_4;
                }

                inline
                double I1_V23_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & /*q2*/) const
                {
                    // contribution proportional to chi_double_bar_4 with a single pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double sigma2    = pow(sigma, 2);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_1 = -(pow(m_B,-2.0)*(m_B*sigmabar*t*tbar*(m_v + 2.0*m_B*(1 + sigma - u)) +
                                       omega*u*(m_v*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)) +
                                       2.0*m_B*sigma2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)) +
                                       sigma*(2.0*m_B*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)))) +
                                       8.0*u2*(1 + sigma)*t*tbar*(-1 + u)*pow(omega,2.0))*pow(sigmabar,-3.0));

                    return C_1 * chi_double_bar_4;
                }

                inline
                double I2_V23_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_double_bar_4 with a double pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_2 = (pow(m_B,-2.0)*pow(sigmabar,-4.0)*(-(2.0*omega*u*
                                       (-(8.0*m_B3*t*tbar*(-1 + u)) + 2.0*m_B3*sigma4*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_v*q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) -
                                       8.0*m_B*t*tbar*(q2 - q2*u + m_v2*(-2.0 + u)) +
                                       m_v3*(-3.0 + 6.0*u + 2.0*t2*(1 + 5.0*u) - 2.0*t*(1 + 5.0*u)) +
                                       m_B2*m_v*(1 - 2.0*u + 2.0*t*(-5.0 + 3.0*u) + t2*(10.0 - 6.0*u)) +
                                       m_B2*sigma3*(m_B*(6.0 + t2*(4.0 - 24*u) - 12.0*u + 4.0*t*(-1 + 6.0*u)) +
                                       3.0*m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) -
                                       sigma*(2.0*m_B*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       2.0*m_B3*(1 + 6.0*t2 - 6.0*t)*(-1 + 2.0*u) -
                                       2.0*m_B*m_v2*(-3.0 + 6.0*u + 2.0*t2*(1 + 6.0*u) - 2.0*t*(1 + 6.0*u)) +
                                       m_v3*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       m_v*q2*(-3.0 + 6.0*u + t*(2.0 - 14.0*u) + 2.0*t2*(-1 + 7.0*u)) +
                                       m_B2*m_v*(3.0 - 6.0*u - 2.0*t2*(-9.0 + 7.0*u) + 2.0*t*(-9.0 + 7.0*u))) +
                                       m_B*sigma2*(m_B*m_v*(5.0 + t2*(14.0 - 26*u) - 10.0*u + 2.0*t*(-7.0 + 13.0*u)) +
                                       2.0*m_B2*(-3.0 + 6.0*u + t*(2.0 - 16.0*u) + 2.0*t2*(-1 + 8.0*u)) +
                                       2.0*(q2*(3.0 - 6.0*u + t2*(2.0 - 16.0*u) + 2.0*t*(-1 + 8.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(2.0 - 8.0*u) + t2*(-2.0 + 8.0*u)))))) +
                                       4.0*u2*(-(4.0*m_B2*t*tbar*(-1 + u)) - 4.0*(3.0*m_v2 - q2)*t*tbar*(-1 + u) -
                                       4.0*(m_v2 - 3.0*q2)*sigma*t*tbar*(-1 + u) +
                                       m_B*m_v*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u) -
                                       2.0*m_B2*sigma*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)) -
                                       m_B*sigma2*(-(20.0*m_B*t*tbar*(-1 + u)) +
                                       m_v*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u)) +
                                       2.0*m_B2*sigma3*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)))*
                                       pow(omega,2.0) - m_B*sigmabar*(2.0*m_v3*t*tbar + 4.0*m_B*(m_B2 - q2)*t*tbar*(-1 + u) -
                                       4.0*m_B*m_v2*t*tbar*(1 + u) - 2.0*m_B2*m_v*t*tbar*(-3.0 + u) -
                                       4.0*m_B3*sigma3*t*tbar*(-1 + 2.0*u) +
                                       m_v*q2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)) +
                                       sigma*(2.0*m_B2*m_v*t*tbar*(-6.0 + u) - 4.0*m_B3*t*tbar*(-3.0 + 4.0*u) +
                                       m_v*q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                       m_B*(4.0*m_v2*t*tbar*(1 + 2.0*u) + q2*(-2.0 + 4.0*t2 - 4.0*t + 4.0*u))) -
                                       2.0*m_B*sigma2*(-(m_B*t*tbar*(3.0*m_v + 2.0*m_B*(-3.0 + 5.0*u))) +
                                       q2*(1 - 2.0*u*pow(1 - 2.0*t,2.0))))))/2.0;

                    return C_2 * chi_double_bar_4;
                }

                inline
                double I3_V23_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4), m_B5 = pow(m_B, 5);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5), sigma6 = pow(sigma, 6);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(pow(m_B,-2.0)*pow(sigmabar,-5.0)*(-(m_B*sigmabar*
                                       (-(4.0*m_B5*sigma5*t*tbar*(-1 + u)) -
                                       2.0*m_B3*sigma4*(-(2.0*m_B*t*tbar*(m_v + 4.0*m_B*(-1 + u))) +
                                       q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       m_v2*(4.0*m_B*m_v2*t*tbar - 4.0*m_B*(m_B2 - q2)*t*tbar*(-1 + u) +
                                       m_v*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u + 4.0*m_B2*t*tbar*(-2.0 + u))) +
                                       2.0*m_B2*sigma3*(-(6.0*m_B2*m_v*t*tbar) - 12.0*m_B3*t*tbar*(-1 + u) +
                                       m_v*q2*(1 + (-2.0 - 6.0*t2 + 6.0*t)*u) +
                                       2.0*m_B*(m_v2*t*(-tbar + 2.0*u - 2.0*t*u) +
                                       q2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)))) +
                                       m_B*sigma2*(12.0*m_B3*m_v*t*tbar + 16.0*m_B4*t*tbar*(-1 + u) +
                                       2.0*q2*(2.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       2.0*m_B2*(-(2.0*m_v2*t*tbar*(-3.0 + 5.0*u)) +
                                       q2*(1 - 2.0*u + 4.0*t2*(-2.0 + u) - 4.0*t*(-2.0 + u))) +
                                       m_B*(-(4.0*m_v3*t*tbar) +
                                       m_v*q2*(-3.0 + 6.0*u + 4.0*t2*(1 + 4.0*u) - 4.0*t*(1 + 4.0*u)))) +
                                       sigma*(-(4.0*m_B4*m_v*t*tbar) - 4.0*m_B5*t*tbar*(-1 + u) +
                                       4.0*m_B3*t*tbar*(2.0*q2*(-1 + u) + m_v2*(-3.0 + 4.0*u)) +
                                       m_B2*m_v*(-(4.0*m_v2*t*tbar*(-3.0 + u)) +
                                       q2*(1 - 2.0*u - 4.0*t2*(1 + u) + 4.0*t*(1 + u))) +
                                       m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(q2,2.0) +
                                       2.0*m_B*(-(2.0*m_v4*t*tbar*u) - m_v2*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) -
                                       2.0*t*tbar*(-1 + u)*pow(q2,2.0))))) -
                                       2.0*omega*u*(2.0*m_B5*sigma6*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_v2*(8.0*m_B3*t*tbar*(-1 + u) - m_v*(m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       4.0*m_B*t*tbar*(-(2.0*q2*(-1 + u)) + m_v2*(-2.0 + u)) +
                                       m_B2*m_v*(-1 + 2.0*u + 4.0*t2*(-4.0 + 3.0*u) - 4.0*t*(-4.0 + 3.0*u))) +
                                       m_B3*sigma4*(-(m_B*m_v*(7.0 + 20.0*t2 - 20.0*t)*(-1 + 2.0*u)) +
                                       8.0*q2*t*(-tbar + u - t*u) + 4.0*m_B2*
                                       (-3.0 + t*(14.0 - 26*u) + 6.0*u + 2.0*t2*(-7.0 + 13.0*u)) +
                                       2.0*m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) -
                                       2.0*m_B4*sigma5*(m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       4.0*m_B*(-1 + 2.0*u + t*(3.0 - 7.0*u) + t2*(-3.0 + 7.0*u))) +
                                       m_B*sigma2*(m_B4*(-2.0 + t*(36 - 44*u) + 4.0*u + 4.0*t2*(-9.0 + 11.0*u)) +
                                       m_B3*m_v*(5.0 - 10.0*u - 4.0*t2*(-7.0 + 6.0*u) + 4.0*t*(-7.0 + 6.0*u)) +
                                       m_B2*(40*q2*t*tbar*(-1 + u) + m_v2*
                                       (6.0 + t2*(36 - 68*u) - 12.0*u + 4.0*t*(-9.0 + 17.0*u))) +
                                       2.0*q2*(q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) +
                                       m_B*(m_v*q2*(1 - 2.0*u - 4.0*t2*(3.0 + 2.0*u) + 4.0*t*(3.0 + 2.0*u)) +
                                       4.0*m_v3*(-1 + 2.0*u + t*(5.0 - 7.0*u) + t2*(-5.0 + 7.0*u)))) +
                                       m_B2*sigma3*(-(8.0*m_B3*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u))) +
                                       3.0*m_B2*m_v*(-3.0 + 6.0*u + 4.0*t2*(-3.0 + 4.0*u) - 4.0*t*(-3.0 + 4.0*u)) +
                                       2.0*m_B*(-(16.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-3.0 + t*(10.0 - 24*u) + 6.0*u + 2.0*t2*(-5.0 + 12.0*u))) -
                                       m_v*(2.0*m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       q2*(1 - 2.0*u - 4.0*t2*(1 + 2.0*u) + t*(4.0 + 8.0*u)))) +
                                       sigma*(8.0*m_B5*t*tbar*(-1 + u) + m_v*(m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_B4*m_v*(-1 + 2.0*u + 4.0*t2*(-2.0 + u) - 4.0*t*(-2.0 + u)) +
                                       2.0*m_B3*(-(8.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-1 + 2.0*u + 2.0*t2*(-7.0 + 10.0*u) + t*(14.0 - 20.0*u))) +
                                       m_B2*(-(8.0*m_v*q2*t*tbar) +
                                       m_v3*(3.0 - 6.0*u - 4.0*t2*(-8.0 + 7.0*u) + 4.0*t*(-8.0 + 7.0*u))) +
                                       m_B*(2.0*m_v2*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_v4*(2.0 - 4.0*u - 4.0*t2*(1 + 2.0*u) + t*(4.0 + 8.0*u)) +
                                       8.0*t*(-tbar + u - t*u)*pow(q2,2.0)))) -
                                       4.0*u2*pow(omega,2.0)*(m_v2*(-(4.0*m_B2*t*tbar*(-1 + u)) -
                                       4.0*(m_v2 - q2)*t*tbar*(-1 + u) +
                                       m_B*m_v*(1 - 8.0*t2*(-1 + u) + 8.0*t*(-1 + u) - 2.0*u)) +
                                       2.0*m_B4*sigma5*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       2.0*m_B3*sigma4*(m_B*(-3.0 + 6.0*u + t*(8.0 - 20.0*u) + 4.0*t2*(-2.0 + 5.0*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       sigma*(-(4.0*m_B4*t*tbar*(-1 + u)) + 4.0*(m_v2 - q2)*q2*t*tbar*(-1 + u) +
                                       m_B3*m_v*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u) +
                                       m_B*m_v*(m_v2*(-1 + 12.0*t2*(-1 + u) - 12.0*t*(-1 + u) + 2.0*u) +
                                       q2*(-1 + 4.0*t2*(-1 + u) - 4.0*t*(-1 + u) + 2.0*u)) -
                                       2.0*m_B2*(-(4.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)))) +
                                       m_B*sigma2*(m_v*q2*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u) +
                                       4.0*m_v3*t*(-tbar + u - t*u) -
                                       2.0*m_B3*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       4.0*m_B2*m_v*(-1 + 2.0*u + t*(3.0 - 5.0*u) + t2*(-3.0 + 5.0*u)) +
                                       2.0*m_B*(q2*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)))) +
                                       m_B2*sigma3*(4.0*m_v2*t*(-tbar + u - t*u) +
                                       m_B*m_v*(5.0 - 10.0*u - 4.0*t2*(-3.0 + 7.0*u) + 4.0*t*(-3.0 + 7.0*u)) +
                                       6.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0) - 2.0*q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))));

                    return C_3 * chi_double_bar_4;
                }

                inline
                double I3d1_V23_chi_double_bar_4(const double & sigma, const double & omega, const double & t, const double & u, const double & q2) const
                {
                    // contribution proportional to chi_double_bar_4 with a triple pole in k2
                    const double m_v      = this->m_v();
                    const double omega_1  = sigma * m_B;
                    const double omega_2  = 2.0 * omega;
                    const double sigmabar = 1.0 - sigma;
                    const double tbar     = 1.0 - t;

                    const double m_v2      = pow(m_v, 2), m_v3 = pow(m_v, 3), m_v4 = pow(m_v, 4);
                    const double m_B2      = pow(m_B, 2), m_B3 = pow(m_B, 3), m_B4 = pow(m_B, 4), m_B5 = pow(m_B, 5), m_B6 = pow(m_B, 6);
                    const double sigma2    = pow(sigma, 2), sigma3 = pow(sigma, 3), sigma4 = pow(sigma, 4), sigma5 = pow(sigma, 5), sigma6 = pow(sigma, 6);
                    const double t2        = pow(t, 2);
                    const double u2        = pow(u, 2);

                    const double chi_bar_4 = this->chi_bar_4(omega_1, omega_2);
                    const double chi_double_bar_4 = this->chi_double_bar_4(omega_1, omega_2);

                    const double C_3 = -(pow(sigmabar,-5.0)*(-(sigmabar*(-(4.0*m_B5*sigma5*t*tbar*(-1 + u)) -
                                       2.0*m_B3*sigma4*(-(2.0*m_B*t*tbar*(m_v + 4.0*m_B*(-1 + u))) +
                                       q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       m_v2*(4.0*m_B*m_v2*t*tbar - 4.0*m_B*(m_B2 - q2)*t*tbar*(-1 + u) +
                                       m_v*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u + 4.0*m_B2*t*tbar*(-2.0 + u))) +
                                       2.0*m_B2*sigma3*(-(6.0*m_B2*m_v*t*tbar) - 12.0*m_B3*t*tbar*(-1 + u) +
                                       m_v*q2*(1 + (-2.0 - 6.0*t2 + 6.0*t)*u) +
                                       2.0*m_B*(m_v2*t*(-tbar + 2.0*u - 2.0*t*u) +
                                       q2*(-1 + 2.0*u + 2.0*t2*(1 + u) - 2.0*t*(1 + u)))) +
                                       m_B*sigma2*(12.0*m_B3*m_v*t*tbar + 16.0*m_B4*t*tbar*(-1 + u) +
                                       2.0*q2*(2.0*m_v2*t*tbar*u + q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                       2.0*m_B2*(-(2.0*m_v2*t*tbar*(-3.0 + 5.0*u)) +
                                       q2*(1 - 2.0*u + 4.0*t2*(-2.0 + u) - 4.0*t*(-2.0 + u))) +
                                       m_B*(-(4.0*m_v3*t*tbar) +
                                       m_v*q2*(-3.0 + 6.0*u + 4.0*t2*(1 + 4.0*u) - 4.0*t*(1 + 4.0*u)))) +
                                       sigma*(-(4.0*m_B4*m_v*t*tbar) - 4.0*m_B5*t*tbar*(-1 + u) +
                                       4.0*m_B3*t*tbar*(2.0*q2*(-1 + u) + m_v2*(-3.0 + 4.0*u)) +
                                       m_B2*m_v*(-(4.0*m_v2*t*tbar*(-3.0 + u)) +
                                       q2*(1 - 2.0*u - 4.0*t2*(1 + u) + 4.0*t*(1 + u))) +
                                       m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(q2,2.0) +
                                       2.0*m_B*(-(2.0*m_v4*t*tbar*u) - m_v2*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) -
                                       2.0*t*tbar*(-1 + u)*pow(q2,2.0))))) -
                                       2.0*omega*u*pow(m_B,-1)*(2.0*m_B5*sigma6*
                                       (-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       m_v2*(8.0*m_B3*t*tbar*(-1 + u) - m_v*(m_v2 - q2)*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       4.0*m_B*t*tbar*(-(2.0*q2*(-1 + u)) + m_v2*(-2.0 + u)) +
                                       m_B2*m_v*(-1 + 2.0*u + 4.0*t2*(-4.0 + 3.0*u) - 4.0*t*(-4.0 + 3.0*u))) +
                                       m_B3*sigma4*(-(m_B*m_v*(7.0 + 20.0*t2 - 20.0*t)*(-1 + 2.0*u)) +
                                       8.0*q2*t*(-tbar + u - t*u) + 4.0*m_B2*
                                       (-3.0 + t*(14.0 - 26*u) + 6.0*u + 2.0*t2*(-7.0 + 13.0*u)) +
                                       2.0*m_v2*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) -
                                       2.0*m_B4*sigma5*(m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u)) +
                                       4.0*m_B*(-1 + 2.0*u + t*(3.0 - 7.0*u) + t2*(-3.0 + 7.0*u))) +
                                       m_B*sigma2*(m_B4*(-2.0 + t*(36 - 44*u) + 4.0*u + 4.0*t2*(-9.0 + 11.0*u)) +
                                       m_B3*m_v*(5.0 - 10.0*u - 4.0*t2*(-7.0 + 6.0*u) + 4.0*t*(-7.0 + 6.0*u)) +
                                       m_B2*(40*q2*t*tbar*(-1 + u) + m_v2*
                                       (6.0 + t2*(36 - 68*u) - 12.0*u + 4.0*t*(-9.0 + 17.0*u))) +
                                       2.0*q2*(q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                       m_v2*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u))) +
                                       m_B*(m_v*q2*(1 - 2.0*u - 4.0*t2*(3.0 + 2.0*u) + 4.0*t*(3.0 + 2.0*u)) +
                                       4.0*m_v3*(-1 + 2.0*u + t*(5.0 - 7.0*u) + t2*(-5.0 + 7.0*u)))) +
                                       m_B2*sigma3*(-(8.0*m_B3*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u))) +
                                       3.0*m_B2*m_v*(-3.0 + 6.0*u + 4.0*t2*(-3.0 + 4.0*u) - 4.0*t*(-3.0 + 4.0*u)) +
                                       2.0*m_B*(-(16.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-3.0 + t*(10.0 - 24*u) + 6.0*u + 2.0*t2*(-5.0 + 12.0*u))) -
                                       m_v*(2.0*m_v2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                       q2*(1 - 2.0*u - 4.0*t2*(1 + 2.0*u) + t*(4.0 + 8.0*u)))) +
                                       sigma*(8.0*m_B5*t*tbar*(-1 + u) + m_v*(m_v2 - q2)*q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                       m_B4*m_v*(-1 + 2.0*u + 4.0*t2*(-2.0 + u) - 4.0*t*(-2.0 + u)) +
                                       2.0*m_B3*(-(8.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-1 + 2.0*u + 2.0*t2*(-7.0 + 10.0*u) + t*(14.0 - 20.0*u))) +
                                       m_B2*(-(8.0*m_v*q2*t*tbar) +
                                       m_v3*(3.0 - 6.0*u - 4.0*t2*(-8.0 + 7.0*u) + 4.0*t*(-8.0 + 7.0*u))) +
                                       m_B*(2.0*m_v2*q2*(1 + 2.0*t2 - 2.0*t)*(-1 + 2.0*u) +
                                       m_v4*(2.0 - 4.0*u - 4.0*t2*(1 + 2.0*u) + t*(4.0 + 8.0*u)) +
                                       8.0*t*(-tbar + u - t*u)*pow(q2,2.0)))) -
                                       4.0*u2*pow(m_B,-1)*pow(omega,2.0)*
                                       (m_v2*(-(4.0*m_B2*t*tbar*(-1 + u)) - 4.0*(m_v2 - q2)*t*tbar*(-1 + u) +
                                       m_B*m_v*(1 - 8.0*t2*(-1 + u) + 8.0*t*(-1 + u) - 2.0*u)) +
                                       2.0*m_B4*sigma5*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) -
                                       2.0*m_B3*sigma4*(m_B*(-3.0 + 6.0*u + t*(8.0 - 20.0*u) + 4.0*t2*(-2.0 + 5.0*u)) +
                                       m_v*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                       sigma*(-(4.0*m_B4*t*tbar*(-1 + u)) + 4.0*(m_v2 - q2)*q2*t*tbar*(-1 + u) +
                                       m_B3*m_v*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u) +
                                       m_B*m_v*(m_v2*(-1 + 12.0*t2*(-1 + u) - 12.0*t*(-1 + u) + 2.0*u) +
                                       q2*(-1 + 4.0*t2*(-1 + u) - 4.0*t*(-1 + u) + 2.0*u)) -
                                       2.0*m_B2*(-(4.0*q2*t*tbar*(-1 + u)) +
                                       m_v2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)))) +
                                       m_B*sigma2*(m_v*q2*(1 - 4.0*t2*(-1 + u) + 4.0*t*(-1 + u) - 2.0*u) +
                                       4.0*m_v3*t*(-tbar + u - t*u) -
                                       2.0*m_B3*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       4.0*m_B2*m_v*(-1 + 2.0*u + t*(3.0 - 5.0*u) + t2*(-3.0 + 5.0*u)) +
                                       2.0*m_B*(q2*(-1 + 2.0*u + t*(8.0 - 12.0*u) + 4.0*t2*(-2.0 + 3.0*u)) +
                                       m_v2*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u)))) +
                                       m_B2*sigma3*(4.0*m_v2*t*(-tbar + u - t*u) +
                                       m_B*m_v*(5.0 - 10.0*u - 4.0*t2*(-3.0 + 7.0*u) + 4.0*t*(-3.0 + 7.0*u)) +
                                       6.0*m_B2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0) - 2.0*q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)))));

                    const double C_3d1 = pow(m_B,-2.0)*pow(sigmabar,-6.0)*(-(4.0*m_B6*t*tbar*(-1 + u)) -
                                         4.0*m_B5*t*tbar*(m_v - 4.0*omega*(-1 + u)*u) -
                                         24*m_B5*sigma5*(-(m_B*t*tbar*(-1 + u)) +
                                         omega*u*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                         4.0*m_B5*sigma6*(-(m_B*t*tbar*(-1 + u)) +
                                         omega*u*(1 - 2.0*u + t2*(2.0 - 6.0*u) + t*(-2.0 + 6.0*u))) +
                                         m_B3*(-32*omega*q2*t*tbar*(-1 + u)*u + 4.0*m_v2*omega*u*(-1 + 6.0*t2 - 6.0*t + 2.0*u) +
                                         4.0*m_v3*t*tbar*(-5.0 + 3.0*u) + m_v*q2*(1 - 2.0*u - 4.0*t2*(1 + u) + 4.0*t*(1 + u)) -
                                         4.0*m_v*u2*(-1 + 4.0*t2*(-1 + u) - 4.0*t*(-1 + u) + 2.0*u)*pow(omega,2.0)) -
                                         2.0*m_B4*(-(2.0*m_v2*t*tbar) + m_v*omega*u*
                                         (1 - 2.0*u - 4.0*t2*(-2.0 + u) + 4.0*t*(-2.0 + u)) -
                                         4.0*t*tbar*(-1 + u)*(q2 - 2.0*u2*pow(omega,2.0))) +
                                         2.0*m_B2*(-(8.0*m_v*omega*q2*t*tbar*u) - 2.0*m_v4*t*tbar*(-4.0 + u) +
                                         2.0*m_v3*omega*u*(-1 + 2.0*u + 8.0*t2*(-3.0 + 2.0*u) - 8.0*t*(-3.0 + 2.0*u)) -
                                         2.0*q2*t*tbar*(-1 + u)*(q2 - 8.0*u2*pow(omega,2.0)) +
                                         m_v2*(q2*(1 - 2.0*u - 2.0*t2*(-5.0 + 6.0*u) + 2.0*t*(-5.0 + 6.0*u)) +
                                         4.0*u2*(1 - 4.0*t2 + 4.0*t - 2.0*u)*pow(omega,2.0))) +
                                         2.0*m_B3*sigma4*(-30*m_B3*t*tbar*(-1 + u) + 8.0*omega*q2*t*u*(-tbar + u - t*u) -
                                         2.0*m_v2*omega*u*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                         2.0*m_B2*(-(m_v*t*tbar) + 2.0*omega*u*
                                         (7.0 + t2*(16.0 - 44*u) - 14.0*u + 4.0*t*(-4.0 + 11.0*u))) +
                                         m_B*(-(2.0*m_v2*t*tbar*(-1 + 2.0*u)) +
                                         m_v*omega*u*(-3.0 + (6.0 + 20.0*t2 - 20.0*t)*u) +
                                         2.0*q2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                         8.0*u2*(-tbar + 2.0*u - 5.0*t*u + t2*(-1 + 5.0*u))*pow(omega,2.0)) +
                                         m_v*(q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                         4.0*u2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))*pow(omega,2.0))) -
                                         2.0*omega*u*(8.0*omega*t*u*(-tbar + u - t*u) + m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u))*
                                         (5.0*m_v4 - 6.0*m_v2*q2 + pow(q2,2.0)) +
                                         m_B*(4.0*m_v2*omega*q2*u*(-1 + t*(22 - 24*u) + 2.0*u + t2*(-22 + 24*u)) +
                                         4.0*m_v4*omega*u*(1 - 2.0*u - 2.0*t2*(-9.0 + 7.0*u) + 2.0*t*(-9.0 + 7.0*u)) +
                                         m_v*q2*(q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                         4.0*u2*(-1 + 4.0*t2*(-1 + u) - 4.0*t*(-1 + u) + 2.0*u)*pow(omega,2.0)) -
                                         4.0*m_v3*(q2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u) +
                                         4.0*u2*(-1 + 7.0*t2*(-1 + u) - 7.0*t*(-1 + u) + 2.0*u)*pow(omega,2.0)) +
                                         16.0*omega*t*tbar*(-1 + u)*u*pow(q2,2.0)) -
                                         2.0*sigma*(-(12.0*m_B6*t*tbar*(-1 + u)) -
                                         4.0*omega*(m_v2 - q2)*q2*u*(8.0*omega*t*u*(-tbar + u - t*u) +
                                         m_v*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)) +
                                         4.0*m_B5*(-(2.0*m_v*t*tbar) + omega*u*
                                         (1 - 2.0*u - 2.0*t2*(-5.0 + 7.0*u) + 2.0*t*(-5.0 + 7.0*u))) -
                                         2.0*m_B3*(-(6.0*m_v3*t*tbar*(-2.0 + u)) + 8.0*omega*q2*t*u*(-tbar + u - t*u) +
                                         2.0*m_v2*omega*u*(1 - 2.0*u + 2.0*t2*(-5.0 + 3.0*u) + t*(10.0 - 6.0*u)) +
                                         m_v*(q2*(-1 + (2.0 + 6.0*t2 - 6.0*t)*u) +
                                         4.0*u2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))*pow(omega,2.0))) +
                                         m_B*(8.0*m_v4*omega*u*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u)) -
                                         4.0*m_v2*omega*q2*u*(-3.0 + 6.0*u + t*(2.0 - 12.0*u) + 2.0*t2*(-1 + 6.0*u)) +
                                         2.0*m_v3*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u +
                                         4.0*u2*(1 - 10.0*t2*(-1 + u) + 10.0*t*(-1 + u) - 2.0*u)*pow(omega,2.0)) +
                                         m_v*q2*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u -
                                         4.0*u2*(-1 + 4.0*t2*(-1 + u) - 4.0*t*(-1 + u) + 2.0*u)*pow(omega,2.0)) +
                                         4.0*omega*u*(-1 + 2.0*u + t*(6.0 - 10.0*u) + 2.0*t2*(-3.0 + 5.0*u))*pow(q2,2.0)) -
                                         2.0*m_B2*(-(2.0*m_v4*t*tbar*(2.0 + u)) +
                                         2.0*m_v3*omega*u*(1 - 2.0*u + t2*(22 - 14.0*u) + 2.0*t*(-11.0 + 7.0*u)) +
                                         q2*(q2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                         4.0*u2*(-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(omega,2.0)) +
                                         m_v2*(q2*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u)) +
                                         4.0*u2*(1 - 2.0*u + t2*(6.0 - 10.0*u) + 2.0*t*(-3.0 + 5.0*u))*pow(omega,2.0))\
                                         - m_v*omega*q2*u*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0)) +
                                         2.0*m_B4*(-(2.0*m_v2*t*tbar*(-2.0 + u)) + q2*(-1 + 4.0*t2 - 4.0*t + 2.0*u) +
                                         m_v*omega*u*(-3.0 + 6.0*u + 4.0*t2*(-3.0 + 4.0*u) - 4.0*t*(-3.0 + 4.0*u)) +
                                         4.0*u2*(-1 + 2.0*u)*pow(omega,2.0)*pow(1 - 2.0*t,2.0))) +
                                         3.0*m_B*sigma2*(-(20.0*m_B5*t*tbar*(-1 + u)) +
                                         4.0*m_v2*omega*q2*u*(-1 + 2.0*u - 2.0*t*(1 + 2.0*u) + t2*(2.0 + 4.0*u)) +
                                         4.0*m_B4*(-(2.0*m_v*t*tbar) + omega*u*
                                         (3.0 + t2*(14.0 - 26*u) - 6.0*u + 2.0*t*(-7.0 + 13.0*u))) +
                                         16.0*m_v3*u2*t*tbar*(-1 + u)*pow(omega,2.0) +
                                         m_v*q2*(q2 + q2*(-2.0 - 4.0*t2 + 4.0*t)*u -
                                         4.0*u2*(-1 + 4.0*t2*(-1 + u) - 4.0*t*(-1 + u) + 2.0*u)*pow(omega,2.0)) +
                                         4.0*m_B3*(-q2 - 2.0*m_v2*t*tbar*(-1 + u) + q2*(2.0 + 4.0*t2 - 4.0*t)*u +
                                         2.0*m_v*omega*u*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u)) +
                                         4.0*u2*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))*pow(omega,2.0)) +
                                         4.0*omega*u*(1 - 2.0*u - 2.0*t2*(1 + u) + 2.0*t*(1 + u))*pow(q2,2.0) +
                                         m_B2*(16.0*omega*q2*t*tbar*(-1 + u)*u + 4.0*m_v3*t*tbar*(-3.0 + u) +
                                         8.0*m_v2*omega*t*tbar*u*(-4.0 + 5.0*u) - m_v*q2*(-1 + 2.0*u)*pow(1 - 2.0*t,2.0) +
                                         4.0*m_v*u2*pow(omega,2.0)*(1 - 2.0*u*pow(1 - 2.0*t,2.0))) -
                                         2.0*m_B*(-(2.0*m_v4*t*tbar*u) - 8.0*m_v*omega*q2*t*tbar*u -
                                         2.0*m_v3*omega*u*(-1 + 8.0*t2*(-1 + u) - 8.0*t*(-1 + u) + 2.0*u) -
                                         2.0*q2*t*tbar*(-1 + u)*(q2 - 8.0*u2*pow(omega,2.0)) -
                                         m_v2*(-1 + 2.0*u)*(q2*(1 + 2.0*t2 - 2.0*t) +
                                         4.0*u2*pow(omega,2.0)*pow(1 - 2.0*t,2.0)))) -
                                         2.0*m_B2*sigma3*(-40*m_B4*t*tbar*(-1 + u) +
                                         8.0*m_B3*(-(m_v*t*tbar) - 4.0*omega*u*
                                         (-1 + 2.0*u + t*(3.0 - 7.0*u) + t2*(-3.0 + 7.0*u))) +
                                         m_B*(-(4.0*m_v3*t*tbar) + 32*omega*q2*t*tbar*(-1 + u)*u -
                                         4.0*m_v2*omega*(1 + 6.0*t2 - 6.0*t)*u*(-1 + 2.0*u) +
                                         m_v*q2*(-1 + 2.0*u + 4.0*t2*(1 + u) - 4.0*t*(1 + u)) +
                                         4.0*m_v*u2*(-1 + 4.0*t2*(-1 + u) - 4.0*t*(-1 + u) + 2.0*u)*pow(omega,2.0)) +
                                         2.0*m_B2*(-(2.0*m_v2*t*tbar*(-2.0 + 3.0*u)) +
                                         q2*(-3.0 + 6.0*u + t*(4.0 - 16.0*u) + 4.0*t2*(-1 + 4.0*u)) +
                                         m_v*omega*u*(-5.0 + t*(4.0 - 32*u) + 10.0*u + 4.0*t2*(-1 + 8.0*u)) +
                                         4.0*u2*(-3.0 + 6.0*u + t*(4.0 - 16.0*u) + 4.0*t2*(-1 + 4.0*u))*pow(omega,2.0)) +
                                         2.0*(2.0*m_v2*omega*u*(-(4.0*omega*t*tbar*(-1 + u)*u) +
                                         m_v*(-1 + 2.0*u + t*(2.0 - 6.0*u) + t2*(-2.0 + 6.0*u))) +
                                         (-1 + (2.0 + 4.0*t2 - 4.0*t)*u)*pow(q2,2.0) +
                                         q2*u*(2.0*m_v2*t*tbar + m_v*omega*
                                         (1 - 2.0*u - 4.0*t2*(1 + 2.0*u) + t*(4.0 + 8.0*u)) +
                                         4.0*u*(-1 + 2.0*u)*pow(omega,2.0)*pow(1 - 2.0*t,2.0)))));

                    return C_3 * chi_bar_4 + C_3d1 * chi_double_bar_4;
                }

                double integrand_V23(const std::array<double, 4> & args, const double & q2) const
                {
                    const double sigma = args[0];
                    const double x     = args[1];
                    const double t     = args[2];
                    const double u     = args[3];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = - f_B() * m_B / (8.0  * M_PI * M_PI * f_V() * m_V());
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_V2 = pow(m_V(), 2);
                    const double M4   = pow(M2, 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_V2) / M2());

                    const double I1 = I1_V23_phi_3(sigma, omega, t, u, q2)             + I1_V23_phi_4(sigma, omega, t, u, q2)
                                    + I1_V23_phi_bar_3(sigma, omega, t, u, q2)         + I1_V23_phi_bar_4(sigma, omega, t, u, q2)
                                    + I1_V23_phi_double_bar_3(sigma, omega, t, u, q2)  + I1_V23_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I1_V23_psi_bar_4(sigma, omega, t, u, q2)         + I1_V23_chi_bar_4(sigma, omega, t, u, q2)
                                    + I1_V23_psi_double_bar_4(sigma, omega, t, u, q2)  + I1_V23_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I2 = I2_V23_phi_bar_3(sigma, omega, t, u, q2)         + I2_V23_phi_bar_4(sigma, omega, t, u, q2)
                                    + I2_V23_phi_double_bar_3(sigma, omega, t, u, q2)  + I2_V23_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I2_V23_psi_bar_4(sigma, omega, t, u, q2)         + I2_V23_chi_bar_4(sigma, omega, t, u, q2)
                                    + I2_V23_psi_double_bar_4(sigma, omega, t, u, q2)  + I2_V23_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3 = I3_V23_phi_double_bar_3(sigma, omega, t, u, q2)  + I3_V23_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I3_V23_psi_double_bar_4(sigma, omega, t, u, q2)  + I3_V23_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result = 0.0;
                    result += - I1;
                    result +=   I2 / M2;
                    result += - I3 / (2.0 * M4);
                    result *=   prefactor * exp;

                    return result;
                }

                double surface_V23(const std::array<double, 3> & args, const double & sigma, const double & q2) const
                {
                    const double x     = args[0];
                    const double t     = args[1];
                    const double u     = args[2];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = - f_B() * m_B / (8.0  * M_PI * M_PI * f_V() * m_V());
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_B2 = pow(m_B(), 2);
                    const double m_v2 = pow(m_v(), 2);
                    const double m_V2 = pow(m_V(), 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_V2) / M2());

                    const double sigmabar = 1.0 - sigma, sigmabar2 = pow(sigmabar, 2);
                    const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
                    const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;

                    const double I2   = I2_V23_phi_bar_3(sigma, omega, t, u, q2)           + I2_V23_phi_bar_4(sigma, omega, t, u, q2)
                                      + I2_V23_phi_double_bar_3(sigma, omega, t, u, q2)    + I2_V23_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I2_V23_psi_bar_4(sigma, omega, t, u, q2)           + I2_V23_chi_bar_4(sigma, omega, t, u, q2)
                                      + I2_V23_psi_double_bar_4(sigma, omega, t, u, q2)    + I2_V23_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3   = I3_V23_phi_double_bar_3(sigma, omega, t, u, q2)    + I3_V23_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3_V23_psi_double_bar_4(sigma, omega, t, u, q2)    + I3_V23_chi_double_bar_4(sigma, omega, t, u, q2);

                    const double I3d1 = I3d1_V23_phi_double_bar_3(sigma, omega, t, u, q2)  + I3d1_V23_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3d1_V23_psi_double_bar_4(sigma, omega, t, u, q2)  + I3d1_V23_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result = 0.0;
                    result +=  eta * I2 / m_B2;
                    result += -0.5 * eta / m_B2 * (I3 / M2() + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
                    result *=  prefactor * exp;

                    return result;
                }

                double integrand_V23_m1(const std::array<double, 4> & args, const double & q2) const
                {
                    const double sigma = args[0];
                    const double x     = args[1];
                    const double t     = args[2];
                    const double u     = args[3];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = - f_B() * m_B / (8.0  * M_PI * M_PI * f_V() * m_V());
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_V2 = pow(m_V(), 2);
                    const double M4   = pow(M2, 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_V2) / M2());

                    const double I1 = I1_V23_phi_3(sigma, omega, t, u, q2)             + I1_V23_phi_4(sigma, omega, t, u, q2)
                                    + I1_V23_phi_bar_3(sigma, omega, t, u, q2)         + I1_V23_phi_bar_4(sigma, omega, t, u, q2)
                                    + I1_V23_phi_double_bar_3(sigma, omega, t, u, q2)  + I1_V23_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I1_V23_psi_bar_4(sigma, omega, t, u, q2)         + I1_V23_chi_bar_4(sigma, omega, t, u, q2)
                                    + I1_V23_psi_double_bar_4(sigma, omega, t, u, q2)  + I1_V23_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I2 = I2_V23_phi_bar_3(sigma, omega, t, u, q2)         + I2_V23_phi_bar_4(sigma, omega, t, u, q2)
                                    + I2_V23_phi_double_bar_3(sigma, omega, t, u, q2)  + I2_V23_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I2_V23_psi_bar_4(sigma, omega, t, u, q2)         + I2_V23_chi_bar_4(sigma, omega, t, u, q2)
                                    + I2_V23_psi_double_bar_4(sigma, omega, t, u, q2)  + I2_V23_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3 = I3_V23_phi_double_bar_3(sigma, omega, t, u, q2)  + I3_V23_phi_double_bar_4(sigma, omega, t, u, q2)
                                    + I3_V23_psi_double_bar_4(sigma, omega, t, u, q2)  + I3_V23_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result1 =    0.0;
                           result1 += - I1;
                           result1 +=   I2 / M2;
                           result1 += - I3 / (2.0 * M4);
                           result1 *=   exp * s(sigma, q2);

                    double result2 =    0.0;
                           result2 += - I2;
                           result2 +=   I3 / M2;
                           result2 *=   exp;

                    return prefactor * (result1 + result2);
                }

                double surface_V23_m1(const std::array<double, 3> & args, const double & sigma, const double & q2) const
                {
                    const double x     = args[0];
                    const double t     = args[1];
                    const double u     = args[2];

                    const double omega = x / (1.0 - x);

                    const double prefactordx = 1.0 / ((1.0 - x) * (1.0 - x));
                    const double prefactorq2 = 1.0 / (m_c() * m_c() - t * (1.0 -t) * (q2 - 2.0 * m_B * u * omega));
                    const double prefactorHS = - f_B() * m_B / (8.0  * M_PI * M_PI * f_V() * m_V());
                    const double prefactor   = prefactordx * prefactorq2 * prefactorHS;

                    const double m_B2 = pow(m_B(), 2);
                    const double m_v2 = pow(m_v(), 2);
                    const double m_V2 = pow(m_V(), 2);
                    const double exp  = std::exp((-s(sigma, q2) + m_V2) / M2());

                    const double sigmabar = 1.0 - sigma, sigmabar2 = pow(sigmabar, 2);
                    const double eta      = 1.0 / (1.0 + (m_v2 - q2) / (sigmabar2 * m_B2));
                    const double etad1    = 2.0 * (eta - 1.0) * eta / sigmabar;

                    const double I2   = I2_V23_phi_bar_3(sigma, omega, t, u, q2)           + I2_V23_phi_bar_4(sigma, omega, t, u, q2)
                                      + I2_V23_phi_double_bar_3(sigma, omega, t, u, q2)    + I2_V23_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I2_V23_psi_bar_4(sigma, omega, t, u, q2)           + I2_V23_chi_bar_4(sigma, omega, t, u, q2)
                                      + I2_V23_psi_double_bar_4(sigma, omega, t, u, q2)    + I2_V23_chi_double_bar_4(sigma, omega, t, u, q2);
                    const double I3   = I3_V23_phi_double_bar_3(sigma, omega, t, u, q2)    + I3_V23_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3_V23_psi_double_bar_4(sigma, omega, t, u, q2)    + I3_V23_chi_double_bar_4(sigma, omega, t, u, q2);

                    const double I3d1 = I3d1_V23_phi_double_bar_3(sigma, omega, t, u, q2)  + I3d1_V23_phi_double_bar_4(sigma, omega, t, u, q2)
                                      + I3d1_V23_psi_double_bar_4(sigma, omega, t, u, q2)  + I3d1_V23_chi_double_bar_4(sigma, omega, t, u, q2);

                    double result1 =   0.0;
                           result1 +=  eta * I2 / m_B2;
                           result1 += -0.5 * eta / m_B2 * (I3 / M2 + eta / m_B2 * I3d1 + I3 * etad1 / m_B2);
                           result1 *=  exp * s(sigma, q2);

                    double result2 =  0.0;
                           result2 += 0.5 * eta * I3 / m_B2;
                           result2 *= exp;

                    return prefactor * (result1 + result2);
                }

                double V23(const double & q2) const
                {
                    const double sigma_0 = this->sigma_0(q2, s0_0_V23(), s0_1_V23());

                    double value_integral = 0.0;
                    double value_surface  = 0.0;

                    const std::function<double (const std::array<double, 4> &)> integrand = std::bind(&LCSR::integrand_V23, this, std::placeholders::_1, q2);
                    const std::function<double (const std::array<double, 3> &)> surface = std::bind(&LCSR::surface_V23, this, std::placeholders::_1, sigma_0, q2);

                    value_integral = integrate(integrand, { 0.0, 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0, 1.0 }, cubature::Config());
                    value_surface  = integrate(surface, { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, cubature::Config());

                    return value_integral + value_surface;
                }

                double normalized_first_moment_V23(const double & q2) const
                {
                    const double sigma_0 = this->sigma_0(q2, s0_0_V23(), s0_1_V23());

                    double value_integral = 0.0;
                    double value_surface  = 0.0;

                    const std::function<double (const std::array<double, 4> &)> integrand = std::bind(&LCSR::integrand_V23, this, std::placeholders::_1, q2);
                    const std::function<double (const std::array<double, 3> &)> surface = std::bind(&LCSR::surface_V23, this, std::placeholders::_1, sigma_0, q2);

                    value_integral = integrate(integrand, { 0.0, 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0, 1.0 }, cubature::Config());
                    value_surface  = integrate(surface, { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, cubature::Config());

                    const double denominator = value_integral + value_surface;

                    double value_integral_m1 = 0.0;
                    double value_surface_m1  = 0.0;

                    const std::function<double (const std::array<double, 4> &)> integrand_m1 = std::bind(&LCSR::integrand_V23_m1, this, std::placeholders::_1, q2);
                    const std::function<double (const std::array<double, 3> &)> surface_m1 = std::bind(&LCSR::surface_V23_m1, this, std::placeholders::_1, sigma_0, q2);

                    value_integral_m1 = integrate(integrand_m1, { 0.0, 0.0, 0.0, 0.0 }, { sigma_0, 1.0, 1.0, 1.0 }, cubature::Config());
                    value_surface_m1  = integrate(surface_m1, { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, cubature::Config());

                    const double numerator = value_integral_m1 + value_surface_m1;

                    return numerator / denominator;
                }

                /* Diagnostics */

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    /* paramaters */
                    results.add({ this->m_v(),                                                     "m_v(mu) in the MSbar scheme"                     });
                    results.add({ this->m_c(),                                                     "m_c(mu) in the MSbar scheme"                     });
                    results.add({ this->f_V(),                                                     "final state decay constant"                      });

                    /* V1 */

                    /* I phi_3(sigma, omega, t, u, q2) */
                    results.add({ this->I1_V1_phi_3(0.02, 0.4, 0.3, 0.7, 2.0),                      "I1_V1_phi_3(0.02, 0.4, 0.3, 0.7, 2.0)"          });
                    results.add({ this->I1_V1_phi_3(0.01, 0.7, 0.9, 0.1, 3.0),                      "I1_V1_phi_3(0.01, 0.7, 0.9, 0.1, 3.0)"          });

                    results.add({ this->I1_V1_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),                  "I1_V1_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I1_V1_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),                  "I1_V1_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)"      });
                    results.add({ this->I2_V1_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),                  "I2_V1_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I2_V1_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),                  "I2_V1_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)"      });

                    results.add({ this->I1_V1_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),           "I1_V1_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I1_V1_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),           "I1_V1_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I2_V1_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),           "I2_V1_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I2_V1_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),           "I2_V1_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3_V1_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),           "I3_V1_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I3_V1_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),           "I3_V1_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3d1_V1_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),         "I3d1_V1_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)" });
                    results.add({ this->I3d1_V1_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),         "I3d1_V1_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)" });

                    /* I phi_4(sigma, omega, t, u, q2) */
                    results.add({ this->I1_V1_phi_4(0.02, 0.4, 0.3, 0.7, 2.0),                      "I1_V1_phi_4(0.02, 0.4, 0.3, 0.7, 2.0)"          });
                    results.add({ this->I1_V1_phi_4(0.01, 0.7, 0.9, 0.1, 3.0),                      "I1_V1_phi_4(0.01, 0.7, 0.9, 0.1, 3.0)"          });

                    results.add({ this->I1_V1_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I1_V1_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I1_V1_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I1_V1_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });
                    results.add({ this->I2_V1_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I2_V1_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I2_V1_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I2_V1_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });

                    results.add({ this->I1_V1_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I1_V1_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I1_V1_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I1_V1_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I2_V1_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I2_V1_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I2_V1_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I2_V1_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3_V1_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I3_V1_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I3_V1_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I3_V1_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3d1_V1_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),         "I3d1_V1_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)" });
                    results.add({ this->I3d1_V1_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),         "I3d1_V1_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)" });

                    /* I psi_4(sigma, omega, t, u, q2) */
                    results.add({ this->I1_V1_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I1_V1_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I1_V1_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I1_V1_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });
                    results.add({ this->I2_V1_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I2_V1_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I2_V1_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I2_V1_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });

                    results.add({ this->I1_V1_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I1_V1_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I1_V1_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I1_V1_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I2_V1_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I2_V1_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I2_V1_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I2_V1_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3_V1_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I3_V1_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I3_V1_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I3_V1_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3d1_V1_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),         "I3d1_V1_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)" });
                    results.add({ this->I3d1_V1_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),         "I3d1_V1_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)" });

                    /* I chi_4(sigma, omega, t, u, q2) */
                    results.add({ this->I1_V1_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I1_V1_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I1_V1_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I1_V1_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });
                    results.add({ this->I2_V1_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I2_V1_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I2_V1_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I2_V1_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });

                    results.add({ this->I1_V1_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I1_V1_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I1_V1_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I1_V1_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I2_V1_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I2_V1_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I2_V1_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I2_V1_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3_V1_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I3_V1_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I3_V1_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I3_V1_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3d1_V1_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),         "I3d1_V1_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)" });
                    results.add({ this->I3d1_V1_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),         "I3d1_V1_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)" });

                    /* V2 */

                    /* I phi_3(sigma, omega, t, u, q2) */
                    results.add({ this->I1_V2_phi_3(0.02, 0.4, 0.3, 0.7, 2.0),                      "I1_V2_phi_3(0.02, 0.4, 0.3, 0.7, 2.0)"          });
                    results.add({ this->I1_V2_phi_3(0.01, 0.7, 0.9, 0.1, 3.0),                      "I1_V2_phi_3(0.01, 0.7, 0.9, 0.1, 3.0)"          });

                    results.add({ this->I1_V2_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),                  "I1_V2_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I1_V2_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),                  "I1_V2_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)"      });
                    results.add({ this->I2_V2_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),                  "I2_V2_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I2_V2_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),                  "I2_V2_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)"      });

                    results.add({ this->I1_V2_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),           "I1_V2_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I1_V2_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),           "I1_V2_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I2_V2_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),           "I2_V2_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I2_V2_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),           "I2_V2_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3_V2_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),           "I3_V2_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I3_V2_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),           "I3_V2_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3d1_V2_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),         "I3d1_V2_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)" });
                    results.add({ this->I3d1_V2_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),         "I3d1_V2_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)" });

                    /* I phi_4(sigma, omega, t, u, q2) */
                    results.add({ this->I1_V2_phi_4(0.02, 0.4, 0.3, 0.7, 2.0),                      "I1_V2_phi_4(0.02, 0.4, 0.3, 0.7, 2.0)"          });
                    results.add({ this->I1_V2_phi_4(0.01, 0.7, 0.9, 0.1, 3.0),                      "I1_V2_phi_4(0.01, 0.7, 0.9, 0.1, 3.0)"          });

                    results.add({ this->I1_V2_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I1_V2_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I1_V2_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I1_V2_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });
                    results.add({ this->I2_V2_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I2_V2_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I2_V2_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I2_V2_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });

                    results.add({ this->I1_V2_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I1_V2_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I1_V2_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I1_V2_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I2_V2_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I2_V2_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I2_V2_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I2_V2_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3_V2_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I3_V2_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I3_V2_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I3_V2_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3d1_V2_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),         "I3d1_V2_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)" });
                    results.add({ this->I3d1_V2_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),         "I3d1_V2_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)" });

                    /* I psi_4(sigma, omega, t, u, q2) */
                    results.add({ this->I1_V2_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I1_V2_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I1_V2_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I1_V2_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });
                    results.add({ this->I2_V2_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I2_V2_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I2_V2_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I2_V2_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });

                    results.add({ this->I1_V2_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I1_V2_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I1_V2_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I1_V2_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I2_V2_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I2_V2_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I2_V2_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I2_V2_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3_V2_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I3_V2_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I3_V2_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I3_V2_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3d1_V2_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),         "I3d1_V2_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)" });
                    results.add({ this->I3d1_V2_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),         "I3d1_V2_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)" });

                    /* I chi_4(sigma, omega, t, u, q2) */
                    results.add({ this->I1_V2_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I1_V2_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I1_V2_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I1_V2_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });
                    results.add({ this->I2_V2_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                  "I2_V2_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"      });
                    results.add({ this->I2_V2_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                  "I2_V2_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"      });

                    results.add({ this->I1_V2_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I1_V2_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I1_V2_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I1_V2_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I2_V2_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I2_V2_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I2_V2_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I2_V2_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3_V2_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),           "I3_V2_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"   });
                    results.add({ this->I3_V2_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),           "I3_V2_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"   });
                    results.add({ this->I3d1_V2_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),         "I3d1_V2_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)" });
                    results.add({ this->I3d1_V2_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),         "I3d1_V2_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)" });

                    /* V23 */

                    /* I phi_3(sigma, omega, t, u, q2) */
                    results.add({ this->I1_V23_phi_3(0.02, 0.4, 0.3, 0.7, 2.0),                     "I1_V23_phi_3(0.02, 0.4, 0.3, 0.7, 2.0)"         });
                    results.add({ this->I1_V23_phi_3(0.01, 0.7, 0.9, 0.1, 3.0),                     "I1_V23_phi_3(0.01, 0.7, 0.9, 0.1, 3.0)"         });

                    results.add({ this->I1_V23_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),                 "I1_V23_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)"     });
                    results.add({ this->I1_V23_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),                 "I1_V23_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)"     });
                    results.add({ this->I2_V23_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),                 "I2_V23_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)"     });
                    results.add({ this->I2_V23_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),                 "I2_V23_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)"     });

                    results.add({ this->I1_V23_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),          "I1_V23_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)"  });
                    results.add({ this->I1_V23_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),          "I1_V23_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)"  });
                    results.add({ this->I2_V23_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),          "I2_V23_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)"  });
                    results.add({ this->I2_V23_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),          "I2_V23_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)"  });
                    results.add({ this->I3_V23_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),          "I3_V23_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)"  });
                    results.add({ this->I3_V23_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),          "I3_V23_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)"  });
                    results.add({ this->I3d1_V23_phi_double_bar_3(0.02, 0.4, 0.3, 0.7, 2.0),        "I3d1_V23_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)"});
                    results.add({ this->I3d1_V23_phi_double_bar_3(0.01, 0.7, 0.9, 0.1, 3.0),        "I3d1_V23_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)"});

                    /* I phi_4(sigma, omega, t, u, q2) */
                    results.add({ this->I1_V23_phi_4(0.02, 0.4, 0.3, 0.7, 2.0),                     "I1_V23_phi_4(0.02, 0.4, 0.3, 0.7, 2.0)"         });
                    results.add({ this->I1_V23_phi_4(0.01, 0.7, 0.9, 0.1, 3.0),                     "I1_V23_phi_4(0.01, 0.7, 0.9, 0.1, 3.0)"         });

                    results.add({ this->I1_V23_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                 "I1_V23_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"     });
                    results.add({ this->I1_V23_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                 "I1_V23_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"     });
                    results.add({ this->I2_V23_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                 "I2_V23_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"     });
                    results.add({ this->I2_V23_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                 "I2_V23_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"     });

                    results.add({ this->I1_V23_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),          "I1_V23_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"  });
                    results.add({ this->I1_V23_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),          "I1_V23_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"  });
                    results.add({ this->I2_V23_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),          "I2_V23_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"  });
                    results.add({ this->I2_V23_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),          "I2_V23_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"  });
                    results.add({ this->I3_V23_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),          "I3_V23_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"  });
                    results.add({ this->I3_V23_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),          "I3_V23_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"  });
                    results.add({ this->I3d1_V23_phi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),        "I3d1_V23_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"});
                    results.add({ this->I3d1_V23_phi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),        "I3d1_V23_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"});

                    /* I psi_4(sigma, omega, t, u, q2) */
                    results.add({ this->I1_V23_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                 "I1_V23_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"     });
                    results.add({ this->I1_V23_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                 "I1_V23_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"     });
                    results.add({ this->I2_V23_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                 "I2_V23_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"     });
                    results.add({ this->I2_V23_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                 "I2_V23_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"     });

                    results.add({ this->I1_V23_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),          "I1_V23_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"  });
                    results.add({ this->I1_V23_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),          "I1_V23_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"  });
                    results.add({ this->I2_V23_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),          "I2_V23_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"  });
                    results.add({ this->I2_V23_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),          "I2_V23_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"  });
                    results.add({ this->I3_V23_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),          "I3_V23_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"  });
                    results.add({ this->I3_V23_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),          "I3_V23_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"  });
                    results.add({ this->I3d1_V23_psi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),        "I3d1_V23_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"});
                    results.add({ this->I3d1_V23_psi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),        "I3d1_V23_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"});

                    /* I chi_4(sigma, omega, t, u, q2) */
                    results.add({ this->I1_V23_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                 "I1_V23_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"     });
                    results.add({ this->I1_V23_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                 "I1_V23_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"     });
                    results.add({ this->I2_V23_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),                 "I2_V23_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)"     });
                    results.add({ this->I2_V23_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),                 "I2_V23_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)"     });

                    results.add({ this->I1_V23_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),          "I1_V23_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"  });
                    results.add({ this->I1_V23_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),          "I1_V23_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"  });
                    results.add({ this->I2_V23_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),          "I2_V23_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"  });
                    results.add({ this->I2_V23_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),          "I2_V23_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"  });
                    results.add({ this->I3_V23_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),          "I3_V23_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"  });
                    results.add({ this->I3_V23_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),          "I3_V23_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"  });
                    results.add({ this->I3d1_V23_chi_double_bar_4(0.02, 0.4, 0.3, 0.7, 2.0),        "I3d1_V23_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)"});
                    results.add({ this->I3d1_V23_chi_double_bar_4(0.01, 0.7, 0.9, 0.1, 3.0),        "I3d1_V23_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)"});

                return results;
                }

                virtual complex<double> H_perp(const double & q2) const
                {
                    const double m_B     = this-> m_B();
                    const double m_B2    = pow(m_B, 2);
                    const double m_B3    = pow(m_B, 3);
                    const double m_V     = this-> m_V();
                    const double m_V2    = pow(m_V, 2);
                    const double Q_c     = 2.0 / 3.0; //Charm quark charge
                    auto         wc      = model->wilson_coefficients_b_to_s(mu_ren(), "mu",false);
                    auto         C_1_EOS = wc.c1();
                    auto         C_2_EOS = wc.c2();
                    auto         C_1_AK  = C_2_EOS - 1.0 / 6.0 * C_1_EOS;

                    return (2.0 * C_1_AK) * sqrt(eos::lambda(m_B2 , m_V2, q2)) * Q_c * this->V1(q2) / (sqrt(2.0) * m_B3);
                }

                virtual complex<double> H_para(const double & q2) const
                {
                    const double m_B     = this-> m_B();
                    const double m_B2    = pow(m_B, 2);
                    const double m_B3    = pow(m_B, 3);
                    const double m_V     = this-> m_V();
                    const double m_V2    = pow(m_V, 2);
                    const double Q_c     = 2.0 / 3.0; //Charm quark charge
                    auto         wc      = model->wilson_coefficients_b_to_s(mu_ren(), "mu",false);
                    auto         C_1_EOS = wc.c1();
                    auto         C_2_EOS = wc.c2();
                    auto         C_1_AK  = C_2_EOS - 1.0 / 6.0 * C_1_EOS;

                    return - (2.0 * C_1_AK) * sqrt(2.0) * Q_c * this->V2(q2) * (m_B2 - m_V2) / m_B3;
                }

                virtual complex<double> H_long(const double & q2) const
                {
                    const double m_B     = this-> m_B();
                    const double m_B2    = pow(m_B, 2);
                    const double m_B4    = pow(m_B, 4);
                    const double m_V     = this-> m_V();
                    const double m_V2    = pow(m_V, 2);
                    const double Q_c     = 2.0 / 3.0; //Charm quark charge
                    auto         wc      = model->wilson_coefficients_b_to_s(mu_ren(), "mu",false);
                    auto         C_1_EOS = wc.c1();
                    auto         C_2_EOS = wc.c2();
                    auto         C_1_AK  = C_2_EOS - 1.0 / 6.0 * C_1_EOS;

                    return - (2.0 * C_1_AK) * q2 * Q_c
                    * (
                        this->V2(q2) * (m_B2 - m_V2) * (m_B2 + 3.0 * m_V2 - q2) - eos::lambda(m_B2 , m_V2, q2) *
                        ((m_B2 - m_V2) / (m_B2 - m_V2 - q2) * (this->V2(q2) - this->V23(q2))) //V3
                    )
                    / (2.0 * m_B4 * m_V * (m_B2 - m_V2));
                }

                virtual complex<double> Hhat_perp(const double & q2) const
                {
                    return 0.0;
                }
                virtual complex<double> Hhat_para(const double & q2) const
                {
                    return 0.0;
                }
                virtual complex<double> Hhat_long(const double & q2) const
                {
                    return 0.0;
                }


                virtual complex<double> normalized_moment_V1(const double & q2) const
                {
                    return this->normalized_first_moment_V1(q2);
                }

                virtual complex<double> normalized_moment_V2(const double & q2) const
                {
                    return this->normalized_first_moment_V2(q2);
                }

                virtual complex<double> normalized_moment_V23(const double & q2) const
                {
                    return this->normalized_first_moment_V23(q2);
                }

                static NonlocalFormFactorPtr<nc::PToV> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nc::PToV>(new LCSR(p, o));
                }

        };

        struct BToKstar
        {
            constexpr static const char * label = "B->K^*";
        };
        constexpr const char * BToKstar::label;

        struct BsToPhi
        {
            constexpr static const char * label = "B_s->phi";
        };
        constexpr const char * BsToPhi::label;




        /*
         * Parametrize the entire formfactor, i.e., both leading and all sub-leading powers as described in [GvDV:2020].
         */
        template <typename Process_>
        class GvDV2020 :
            public NonlocalFormFactor<nc::PToV>
        {
            private:
                std::shared_ptr<FormFactors<PToV>> form_factors;

                // spectator quark option
                SwitchOption opt_q;

                //Polynomial expansion parameters
                UsedParameter re_alpha_0_perp;
                UsedParameter im_alpha_0_perp;
                UsedParameter re_alpha_1_perp;
                UsedParameter im_alpha_1_perp;
                UsedParameter re_alpha_2_perp;
                UsedParameter im_alpha_2_perp;

                UsedParameter re_alpha_0_para;
                UsedParameter im_alpha_0_para;
                UsedParameter re_alpha_1_para;
                UsedParameter im_alpha_1_para;
                UsedParameter re_alpha_2_para;
                UsedParameter im_alpha_2_para;

                UsedParameter re_alpha_0_long;
                UsedParameter im_alpha_0_long;
                UsedParameter re_alpha_1_long;
                UsedParameter im_alpha_1_long;
                UsedParameter re_alpha_2_long;
                UsedParameter im_alpha_2_long;

                //Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;

                // final state meson parameters
                UsedParameter m_V;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;

                std::string _final_state() const
                {
                    switch (opt_q.value()[0])
                    {
                        case 's':
                            return "phi";
                            break;

                        default:
                            return "K_d^*";
                    }
                }

            public:
                GvDV2020(const Parameters & p, const Options & o) :
                    form_factors(FormFactorFactory<PToV>::create(stringify(Process_::label) + "::" + o.get("form-factors", "BSZ2015"), p)),
                    opt_q(o, "q", { "d", "s" }, "d"),

                    re_alpha_0_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_0^perp}@GvDV2020"], *this),
                    im_alpha_0_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_0^perp}@GvDV2020"], *this),
                    re_alpha_1_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_1^perp}@GvDV2020"], *this),
                    im_alpha_1_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_1^perp}@GvDV2020"], *this),
                    re_alpha_2_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_2^perp}@GvDV2020"], *this),
                    im_alpha_2_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_2^perp}@GvDV2020"], *this),

                    re_alpha_0_para(p[stringify(Process_::label) + "ccbar::Re{alpha_0^para}@GvDV2020"], *this),
                    im_alpha_0_para(p[stringify(Process_::label) + "ccbar::Im{alpha_0^para}@GvDV2020"], *this),
                    re_alpha_1_para(p[stringify(Process_::label) + "ccbar::Re{alpha_1^para}@GvDV2020"], *this),
                    im_alpha_1_para(p[stringify(Process_::label) + "ccbar::Im{alpha_1^para}@GvDV2020"], *this),
                    re_alpha_2_para(p[stringify(Process_::label) + "ccbar::Re{alpha_2^para}@GvDV2020"], *this),
                    im_alpha_2_para(p[stringify(Process_::label) + "ccbar::Im{alpha_2^para}@GvDV2020"], *this),

                    re_alpha_0_long(p[stringify(Process_::label) + "ccbar::Re{alpha_0^long}@GvDV2020"], *this),
                    im_alpha_0_long(p[stringify(Process_::label) + "ccbar::Im{alpha_0^long}@GvDV2020"], *this),
                    re_alpha_1_long(p[stringify(Process_::label) + "ccbar::Re{alpha_1^long}@GvDV2020"], *this),
                    im_alpha_1_long(p[stringify(Process_::label) + "ccbar::Im{alpha_1^long}@GvDV2020"], *this),
                    re_alpha_2_long(p[stringify(Process_::label) + "ccbar::Re{alpha_2^long}@GvDV2020"], *this),
                    im_alpha_2_long(p[stringify(Process_::label) + "ccbar::Im{alpha_2^long}@GvDV2020"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_" + opt_q.value()], *this),

                    m_V(p["mass::" + _final_state()], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),

                    chiOPE(p["b->sccbar::chiOPE@GvDV2020"], *this)
                {
                    this->uses(*form_factors);
                }

                ~GvDV2020() = default;

                inline complex<double> phi(const double & q2, const unsigned phiParam[4]) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF                        a    b    c    d
                    // 0(P->P) aka plus          3    3    2    2
                    // perp(P->V) = par(P->V)    3    1    3    0
                    // 0(P->V) aka long          3    1    2    2

                    const double m_V2  = pow(m_V, 2);
                    const double m_B2  = pow(m_B, 2),  m_B4 =  pow(m_B, 4);
                    const double m_D02 = pow(m_D0, 2), m_D04 = pow(m_D0, 4);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nc_utils::z(q2, 4.0 * pow(m_D0, 2), s_0);
                    const double Q2   = this->t_s();
                    const double chi = this->chiOPE();

                    const double a = phiParam[0], b = phiParam[1], c = phiParam[2], d = phiParam[3];

                    const double Nlambda = 4*M_PI * pow(m_B2, 0.5*(a-b+c+d) - 1.) * pow( 2*(4*m_D02-s_0)/3/chi,0.5); //(C6)
                    const complex<double> phi1 = -pow(2*pow((4*m_D02-Q2)*(4*m_D02-s_0), 0.5) + 8*m_D02 - Q2 - s_0, 0.5) / 
                                                (2*pow((4*m_D02-Q2)*(4*m_D02-s_0), 0.5) + 8*m_D02 + Q2*(z-1.) - s_0*(z+1.)); //(C7)
                    const complex<double> phi2 = pow(m_B4*pow(z-1., 4) - 2*m_B2*pow(z-1., 2)*(-16*m_D02*z + m_V2*pow(z-1., 2) + s_0*pow(z+1., 2)) + 
                                                pow(16*m_D02*z + m_V2*pow(z-1., 2) - s_0*pow(z+1., 2), 2), 0.5);//(C8)
                    const complex<double> phi3 = pow(8*m_D02 + 4*pow(4*m_D04-s_0*m_D02, 0.5) - s_0, 0.5)/(-8*m_D02 - 4*pow(4*m_D04-s_0*m_D02, 0.5) + s_0*(z+1.)); //(C9)
                    const complex<double> phi4 = pow(s_0 * pow(z+1., 2.) - 16.*z*m_D02, -0.5); //(C10)

                    return Nlambda * pow(1.+z, 0.5) * pow(1.-z, a-b+c+d-1.5) * pow(phi1, a) * pow(phi2, 0.5*b) * pow(phi3, c) * pow(phi4, d); //(C5)
                }

                inline complex<double> H_residue_jpsi(const unsigned phiParam[4], const complex<double> & alpha_0, const complex<double> & alpha_1,
                                                      const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto zBV     = eos::nc_utils::z(pow(m_B+m_V, 2), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(m_Jpsi2,         s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(m_psi2S2,        s_p, s_0);

                    return eos::nc_utils::PGvDV2020(z_Jpsi, zBV, alpha_0, alpha_1, alpha_2) / phi(m_Jpsi2, phiParam) *
                            (1 - norm(z_Jpsi)) * (1. - z_Jpsi*std::conj(z_psi2S)) / (z_Jpsi - z_psi2S);
                };

                inline complex<double> H_residue_psi2s(const unsigned phiParam[4], const complex<double> & alpha_0, const complex<double> & alpha_1,
                                                       const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto zBV     = eos::nc_utils::z(pow(m_B+m_V, 2), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(m_Jpsi2,         s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(m_psi2S2,        s_p, s_0);

                    return eos::nc_utils::PGvDV2020(z_psi2S, zBV, alpha_0, alpha_1, alpha_2) / phi(m_psi2S2, phiParam) *
                            (1 - norm(z_psi2S)) * (1. - z_psi2S*std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi);
                };

            
                virtual complex<double> H_perp(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_perp, im_alpha_0_perp);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_perp, im_alpha_1_perp);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_perp, im_alpha_2_perp);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,              s_p, s_0);
                    const auto zBV     = eos::nc_utils::z(pow(m_B+m_V, 2), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(pow(m_Jpsi, 2),  s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(pow(m_psi2S, 2), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nc_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[4] = {3, 1, 3, 0};

                    return eos::nc_utils::PGvDV2020(z, zBV, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                }

                virtual complex<double> Hhat_perp(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_perp, im_alpha_0_perp);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_perp, im_alpha_1_perp);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_perp, im_alpha_2_perp);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,              s_p, s_0);
                    const auto zBV     = eos::nc_utils::z(pow(m_B+m_V, 2), s_p, s_0);

                    return eos::nc_utils::PGvDV2020(z, zBV, alpha_0, alpha_1, alpha_2);
                }

                virtual complex<double> H_para(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_para, im_alpha_0_para);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_para, im_alpha_1_para);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_para, im_alpha_2_para);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,              s_p, s_0);
                    const auto zBV     = eos::nc_utils::z(pow(m_B+m_V, 2), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(pow(m_Jpsi, 2),  s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(pow(m_psi2S, 2), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nc_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[4] = {3, 1, 3, 0};

                    return eos::nc_utils::PGvDV2020(z, zBV, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                }

                virtual complex<double> Hhat_para(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_para, im_alpha_0_para);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_para, im_alpha_1_para);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_para, im_alpha_2_para);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,              s_p, s_0);
                    const auto zBV     = eos::nc_utils::z(pow(m_B+m_V, 2), s_p, s_0);

                    return eos::nc_utils::PGvDV2020(z, zBV, alpha_0, alpha_1, alpha_2);
                }

                virtual complex<double> H_long(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_long, im_alpha_0_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_long, im_alpha_1_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_long, im_alpha_2_long);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,              s_p, s_0);
                    const auto zBV     = eos::nc_utils::z(pow(m_B+m_V, 2), s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(pow(m_Jpsi, 2),  s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(pow(m_psi2S, 2), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nc_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[4] = {3, 1, 2, 2};

                    return eos::nc_utils::PGvDV2020(z, zBV, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                }

                virtual complex<double> Hhat_long(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_long, im_alpha_0_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_long, im_alpha_1_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_long, im_alpha_2_long);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,              s_p, s_0);
                    const auto zBV     = eos::nc_utils::z(pow(m_B+m_V, 2), s_p, s_0);

                    return eos::nc_utils::PGvDV2020(z, zBV, alpha_0, alpha_1, alpha_2);
                }


                virtual complex<double> H_perp_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_perp, im_alpha_0_perp);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_perp, im_alpha_1_perp);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_perp, im_alpha_2_perp);

                    const unsigned phiParam[4] = {3, 1, 3, 0};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_perp_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_perp, im_alpha_0_perp);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_perp, im_alpha_1_perp);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_perp, im_alpha_2_perp);

                    const unsigned phiParam[4] = {3, 1, 3, 0};

                    return H_residue_psi2s(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_para_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_para, im_alpha_0_para);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_para, im_alpha_1_para);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_para, im_alpha_2_para);

                    const unsigned phiParam[4] = {3, 1, 3, 0};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_para_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_para, im_alpha_0_para);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_para, im_alpha_1_para);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_para, im_alpha_2_para);

                    const unsigned phiParam[4] = {3, 1, 3, 0};

                    return H_residue_psi2s(phiParam, alpha_0, alpha_1, alpha_2);
                };

                   virtual complex<double> H_long_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_long, im_alpha_0_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_long, im_alpha_1_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_long, im_alpha_2_long);

                    const unsigned phiParam[4] = {3, 1, 2, 2};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_long_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_long, im_alpha_0_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_long, im_alpha_1_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_long, im_alpha_2_long);

                    const unsigned phiParam[4] = {3, 1, 2, 2};

                    return H_residue_psi2s(phiParam, alpha_0, alpha_1, alpha_2);
                };


                virtual complex<double> normalized_moment_V1(const double & q2) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V2(const double & q2) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V23(const double & q2) const
                {
                    return 0.0;
                }


                static NonlocalFormFactorPtr<nc::PToV> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nc::PToV>(new GvDV2020<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    const unsigned phiParamlong[4] = {3, 1, 2, 2}; //long polarization
                    results.add({ real(1./this->phi(0.0, phiParamlong)), "Re{1/phi_long(q2 = 0.0)}" });
                    results.add({ imag(1./this->phi(0.0, phiParamlong)), "Im{1/phi_long(q2 = 0.0)}" });
                    results.add({ real(this->phi(16.0, phiParamlong)), "Re{phi_long(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phiParamlong)), "Im{phi_long(q2 = 16.0)}" });

                    const unsigned phiParamperp[4] = {3, 1, 3, 0}; //perp or para polarization
                    results.add({ real(this->phi(16.0, phiParamperp)), "Re{phi_perp(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phiParamperp)), "Im{phi_perp(q2 = 16.0)}" });

                    return results;
                }
        };



        /*
         * Parametrize the entire formfactor, i.e., both leading and all sub-leading powers as described in [GRvDV:2021].
         */
        template <typename Process_>
        class GRvDV2021 :
            public NonlocalFormFactor<nc::PToV>
        {
            private:
                std::shared_ptr<FormFactors<PToV>> form_factors;

                // spectator quark option
                SwitchOption opt_q;

                //Polynomial expansion parameters
                UsedParameter re_alpha_0_perp;
                UsedParameter im_alpha_0_perp;
                UsedParameter re_alpha_1_perp;
                UsedParameter im_alpha_1_perp;
                UsedParameter re_alpha_2_perp;
                UsedParameter im_alpha_2_perp;

                UsedParameter re_alpha_0_para;
                UsedParameter im_alpha_0_para;
                UsedParameter re_alpha_1_para;
                UsedParameter im_alpha_1_para;
                UsedParameter re_alpha_2_para;
                UsedParameter im_alpha_2_para;

                UsedParameter re_alpha_0_long;
                UsedParameter im_alpha_0_long;
                UsedParameter re_alpha_1_long;
                UsedParameter im_alpha_1_long;
                UsedParameter re_alpha_2_long;
                UsedParameter im_alpha_2_long;

                //Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;
                UsedParameter m_Bsst;

                // final state meson parameters
                UsedParameter m_V;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;

                std::string _final_state() const
                {
                    switch (opt_q.value()[0])
                    {
                        case 's':
                            return "phi";
                            break;

                        default:
                            return "K_d^*";
                    }
                }

            public:
                GRvDV2021(const Parameters & p, const Options & o) :
                    form_factors(FormFactorFactory<PToV>::create(stringify(Process_::label) + "::" + o.get("form-factors", "BSZ2015"), p)),
                    opt_q(o, "q", { "d", "s" }, "d"),

                    re_alpha_0_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_0^perp}@GRvDV2021"], *this),
                    im_alpha_0_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_0^perp}@GRvDV2021"], *this),
                    re_alpha_1_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_1^perp}@GRvDV2021"], *this),
                    im_alpha_1_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_1^perp}@GRvDV2021"], *this),
                    re_alpha_2_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_2^perp}@GRvDV2021"], *this),
                    im_alpha_2_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_2^perp}@GRvDV2021"], *this),

                    re_alpha_0_para(p[stringify(Process_::label) + "ccbar::Re{alpha_0^para}@GRvDV2021"], *this),
                    im_alpha_0_para(p[stringify(Process_::label) + "ccbar::Im{alpha_0^para}@GRvDV2021"], *this),
                    re_alpha_1_para(p[stringify(Process_::label) + "ccbar::Re{alpha_1^para}@GRvDV2021"], *this),
                    im_alpha_1_para(p[stringify(Process_::label) + "ccbar::Im{alpha_1^para}@GRvDV2021"], *this),
                    re_alpha_2_para(p[stringify(Process_::label) + "ccbar::Re{alpha_2^para}@GRvDV2021"], *this),
                    im_alpha_2_para(p[stringify(Process_::label) + "ccbar::Im{alpha_2^para}@GRvDV2021"], *this),

                    re_alpha_0_long(p[stringify(Process_::label) + "ccbar::Re{alpha_0^long}@GRvDV2021"], *this),
                    im_alpha_0_long(p[stringify(Process_::label) + "ccbar::Im{alpha_0^long}@GRvDV2021"], *this),
                    re_alpha_1_long(p[stringify(Process_::label) + "ccbar::Re{alpha_1^long}@GRvDV2021"], *this),
                    im_alpha_1_long(p[stringify(Process_::label) + "ccbar::Im{alpha_1^long}@GRvDV2021"], *this),
                    re_alpha_2_long(p[stringify(Process_::label) + "ccbar::Re{alpha_2^long}@GRvDV2021"], *this),
                    im_alpha_2_long(p[stringify(Process_::label) + "ccbar::Im{alpha_2^long}@GRvDV2021"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_" + opt_q.value()], *this),
                    m_Bsst(p["mass::B_s^*"], *this),

                    m_V(p["mass::" + _final_state()], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),

                    chiOPE(p["b->sccbar::chiOPE@GRvDV2021"], *this)
                {
                    this->uses(*form_factors);
                }

                ~GRvDV2021() = default;

                inline complex<double> phi(const double & q2, const unsigned phiParam[5]) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF                        a    b    c    d    e
                    // 0(P->P) aka plus          5    3    2    2    2
                    // perp(P->V) = par(P->V)    5    1    3    0    2
                    // 0(P->V) aka long          5    1    2    2    2

                    const double m_V2  = pow(m_V, 2);
                    const double m_Bsst2 = pow(m_Bsst, 2);
                    const double m_B2  = pow(m_B, 2),  m_B4 =  pow(m_B, 4);
                    const double m_D02 = pow(m_D0, 2), m_D04 = pow(m_D0, 4);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nc_utils::z(q2, 4.0 * pow(m_D0, 2), s_0);
                    const double Q2   = this->t_s();
                    const double chi = this->chiOPE();

                    const double a = phiParam[0], b = phiParam[1], c = phiParam[2], d = phiParam[3], e = phiParam[4];

                    const complex<double> Nlambda = 4. * M_PI * pow(m_B2, 0.5*(a-b+c+d-e) - 1.) * pow( 2.*(4.*m_D02-s_0)/3./chi, 0.5);
                    const complex<double> phi1 = -pow(2.*pow((4.*m_D02-Q2)*(4.*m_D02-s_0), 0.5) + 8.*m_D02 - Q2 - s_0, 0.5) /
                                                (2.*pow((4.*m_D02-Q2)*(4.*m_D02-s_0), 0.5) + 8.*m_D02 + Q2*(z-1.) - s_0*(z+1.));
                    const complex<double> phi2 = pow(m_B4*pow(z-1., 4.) - 2.*m_B2*pow(z-1., 2)*(-16*m_D02*z + m_V2*pow(z-1., 2) + s_0*pow(z+1., 2)) +
                                                pow(16*m_D02*z + m_V2*pow(z-1., 2) - s_0*pow(z+1., 2), 2), 0.5);
                    const complex<double> phi3 = pow(8*m_D02 + 4*pow(4*m_D04-s_0*m_D02, 0.5) - s_0, 0.5)/(-8*m_D02 - 4*pow(4*m_D04-s_0*m_D02, 0.5) + s_0*(z+1.));
                    const complex<double> phi4 = pow(s_0 * pow(z+1., 2.) - 16.*z*m_D02, -0.5);
                    const complex<double> phi5 = pow(s_0 * pow(z+1., 2.) - 16.*z*m_D02 - m_Bsst2*pow(-z+1., 2.), 0.5);

                    return Nlambda * pow(1.+z, 0.5) * pow(1.-z, a-b+c+d-e-1.5) * pow(phi1, a) * pow(phi2, 0.5*b) * pow(phi3, c) * pow(phi4, d) * pow(phi5, e);
                }

                inline complex<double> H_residue_jpsi(const unsigned phiParam[5], const complex<double> & alpha_0, const complex<double> & alpha_1,
                                                      const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z_Jpsi  = eos::nc_utils::z(m_Jpsi2,         s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(m_psi2S2,        s_p, s_0);

                    return eos::nc_utils::P(z_Jpsi, alpha_0, alpha_1, alpha_2) / phi(m_Jpsi2, phiParam) *
                            (1 - norm(z_Jpsi)) * (1. - z_Jpsi*std::conj(z_psi2S)) / (z_Jpsi - z_psi2S);
                };

                inline complex<double> H_residue_psi2s(const unsigned phiParam[5], const complex<double> & alpha_0, const complex<double> & alpha_1,
                                                       const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z_Jpsi  = eos::nc_utils::z(m_Jpsi2,         s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(m_psi2S2,        s_p, s_0);

                    return eos::nc_utils::P(z_psi2S, alpha_0, alpha_1, alpha_2) / phi(m_psi2S2, phiParam) *
                            (1 - norm(z_psi2S)) * (1. - z_psi2S*std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi);
                };

            
                virtual complex<double> H_perp(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_perp, im_alpha_0_perp);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_perp, im_alpha_1_perp);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_perp, im_alpha_2_perp);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,              s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(pow(m_Jpsi, 2),  s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(pow(m_psi2S, 2), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nc_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[5] = {5, 1, 3, 0, 2};

                    return eos::nc_utils::P(z, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                }

                virtual complex<double> Hhat_perp(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_perp, im_alpha_0_perp);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_perp, im_alpha_1_perp);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_perp, im_alpha_2_perp);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2, s_p, s_0);

                    return eos::nc_utils::P(z, alpha_0, alpha_1, alpha_2);
                }

                virtual complex<double> H_para(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_para, im_alpha_0_para);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_para, im_alpha_1_para);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_para, im_alpha_2_para);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,              s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(pow(m_Jpsi, 2),  s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(pow(m_psi2S, 2), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nc_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[5] = {5, 1, 3, 0, 2};

                    return eos::nc_utils::P(z, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                }

                virtual complex<double> Hhat_para(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_para, im_alpha_0_para);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_para, im_alpha_1_para);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_para, im_alpha_2_para);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,              s_p, s_0);

                    return eos::nc_utils::P(z, alpha_0, alpha_1, alpha_2);
                }

                virtual complex<double> H_long(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_long, im_alpha_0_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_long, im_alpha_1_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_long, im_alpha_2_long);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2,              s_p, s_0);
                    const auto z_Jpsi  = eos::nc_utils::z(pow(m_Jpsi, 2),  s_p, s_0);
                    const auto z_psi2S = eos::nc_utils::z(pow(m_psi2S, 2), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nc_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[5] = {5, 1, 2, 2, 2};

                    return eos::nc_utils::P(z, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                }

                virtual complex<double> Hhat_long(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_long, im_alpha_0_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_long, im_alpha_1_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_long, im_alpha_2_long);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nc_utils::z(q2, s_p, s_0);

                    return eos::nc_utils::P(z, alpha_0, alpha_1, alpha_2);
                }



                virtual complex<double> H_perp_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_perp, im_alpha_0_perp);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_perp, im_alpha_1_perp);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_perp, im_alpha_2_perp);

                    const unsigned phiParam[5] = {5, 1, 3, 0, 2};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_perp_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_perp, im_alpha_0_perp);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_perp, im_alpha_1_perp);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_perp, im_alpha_2_perp);

                    const unsigned phiParam[5] = {5, 1, 3, 0, 2};

                    return H_residue_psi2s(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_para_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_para, im_alpha_0_para);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_para, im_alpha_1_para);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_para, im_alpha_2_para);

                    const unsigned phiParam[5] = {5, 1, 3, 0, 2};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_para_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_para, im_alpha_0_para);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_para, im_alpha_1_para);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_para, im_alpha_2_para);

                    const unsigned phiParam[5] = {5, 1, 3, 0, 2};

                    return H_residue_psi2s(phiParam, alpha_0, alpha_1, alpha_2);
                };

                   virtual complex<double> H_long_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_long, im_alpha_0_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_long, im_alpha_1_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_long, im_alpha_2_long);

                    const unsigned phiParam[5] = {5, 1, 2, 2, 2};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> H_long_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_long, im_alpha_0_long);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_long, im_alpha_1_long);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_long, im_alpha_2_long);

                    const unsigned phiParam[5] = {5, 1, 2, 2, 2};

                    return H_residue_psi2s(phiParam, alpha_0, alpha_1, alpha_2);
                };


                virtual complex<double> normalized_moment_V1(const double & q2) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V2(const double & q2) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V23(const double & q2) const
                {
                    return 0.0;
                }


                static NonlocalFormFactorPtr<nc::PToV> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nc::PToV>(new GRvDV2021<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    const unsigned phiParamlong[5] = {5, 1, 2, 2, 2}; //long polarization
                    results.add({ real(this->phi(16.0, phiParamlong)), "Re{phi_long(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phiParamlong)), "Im{phi_long(q2 = 16.0)}" });

                    const unsigned phiParamperp[5] = {5, 1, 3, 0, 2}; //perp or para polarization
                    results.add({ real(this->phi(16.0, phiParamperp)), "Re{phi_perp(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phiParamperp)), "Im{phi_perp(q2 = 16.0)}" });

                    return results;
                }
        };
    }

    NonlocalFormFactorPtr<nc::PToV>
    NonlocalFormFactor<nc::PToV>::make(const QualifiedName & name, const Parameters & p, const Options & o)
    {
        typedef QualifiedName KeyType;
        typedef std::function<NonlocalFormFactorPtr<nc::PToV> (const Parameters &, const Options &)> ValueType;
        std::map<KeyType, ValueType> entries
        {
            // trivial
            std::make_pair("B->K^*::naive",         &nc_p_to_v::Naive::make),
            // analytic
            std::make_pair("B->K^*::QCDF",          &nc_p_to_v::QCDF::make),
            std::make_pair("B->K^*::LCSR",          &nc_p_to_v::LCSR::make),
            std::make_pair("B_s->phi::LCSR",        &nc_p_to_v::LCSR::make),
            // parametrizations
            std::make_pair("B->K^*::GvDV2020",      &nc_p_to_v::GvDV2020<nc_p_to_v::BToKstar>::make),
            std::make_pair("B->K^*::GRvDV2021",     &nc_p_to_v::GRvDV2021<nc_p_to_v::BToKstar>::make),
        };

        auto i = entries.find(name);
        if (entries.end() == i)
        {
            return NonlocalFormFactorPtr<nc::PToV>(nullptr);
        }

        return i->second(p, o);
    }

    template <typename Process_>
    struct Implementation<NonlocalFormFactorObservable<Process_, nc::PToV>>
    {
        NameOption opt_formfactor;
        NonlocalFormFactorPtr<nc::PToV> nc;

        Implementation(const Parameters & p, const Options & o) :
            opt_formfactor(o, "formfactor", qnp::Name("GvDV2020")),
            nc(NonlocalFormFactor<nc::PToV>::make(QualifiedName(qnp::Prefix(Process_::prefix), opt_formfactor.value()), p, o))
        {
        }

        ~Implementation() = default;
    };

    template <typename Process_>
    NonlocalFormFactorObservable<Process_, nc::PToV>::NonlocalFormFactorObservable(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<NonlocalFormFactorObservable<Process_, nc::PToV>>(new Implementation<NonlocalFormFactorObservable<Process_, nc::PToV>>(p, o))
    {
    }

    template <typename Process_>
    NonlocalFormFactorObservable<Process_, nc::PToV>::~NonlocalFormFactorObservable() = default;

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::re_H_perp(const double & q2) const
    {
        return real(this->_imp->nc->H_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::im_H_perp(const double & q2) const
    {
        return imag(this->_imp->nc->H_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::abs_H_perp(const double & q2) const
    {
        return abs(this->_imp->nc->H_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::re_Hhat_perp(const double & q2) const
    {
        return real(this->_imp->nc->Hhat_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::im_Hhat_perp(const double & q2) const
    {
        return imag(this->_imp->nc->Hhat_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::abs_Hhat_perp(const double & q2) const
    {
        return abs(this->_imp->nc->Hhat_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::re_H_para(const double & q2) const
    {
        return real(this->_imp->nc->H_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::im_H_para(const double & q2) const
    {
        return imag(this->_imp->nc->H_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::abs_H_para(const double & q2) const
    {
        return abs(this->_imp->nc->H_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::re_Hhat_para(const double & q2) const
    {
        return real(this->_imp->nc->Hhat_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::im_Hhat_para(const double & q2) const
    {
        return imag(this->_imp->nc->Hhat_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::abs_Hhat_para(const double & q2) const
    {
        return abs(this->_imp->nc->Hhat_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::re_H_long(const double & q2) const
    {
        return real(this->_imp->nc->H_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::im_H_long(const double & q2) const
    {
        return imag(this->_imp->nc->H_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::abs_H_long(const double & q2) const
    {
        return abs(this->_imp->nc->H_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::re_Hhat_long(const double & q2) const
    {
        return real(this->_imp->nc->Hhat_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::im_Hhat_long(const double & q2) const
    {
        return imag(this->_imp->nc->Hhat_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::abs_Hhat_long(const double & q2) const
    {
        return abs(this->_imp->nc->Hhat_long(q2));
    }


    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::re_normalized_moment_V1(const double & q2) const
    {
        return real(this->_imp->nc->normalized_moment_V1(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::re_normalized_moment_V2(const double & q2) const
    {
        return real(this->_imp->nc->normalized_moment_V2(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nc::PToV>::re_normalized_moment_V23(const double & q2) const
    {
        return real(this->_imp->nc->normalized_moment_V23(q2));
    }

    namespace nc
    {
        struct BToKstar
        {
            static constexpr const char * prefix = "B->K^*";
        };
    }

    template class NonlocalFormFactorObservable<nc::BToKstar, nc::PToV>;

    namespace nc
    {
        struct BsToPhi
        {
            static constexpr const char * prefix = "B_s->phi";
        };
    }

    template class NonlocalFormFactorObservable<nc::BsToPhi, nc::PToV>;
}
