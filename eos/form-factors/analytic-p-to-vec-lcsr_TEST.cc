/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#include <test/test.hh>
#include <eos/form-factors/analytic-p-to-vec-lcsr.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>

#include <vector>
#include <iostream>
#include <utility>

using namespace test;
using namespace eos;

class PToVLCSRFormFactorsTest :
    public TestCase
{
    public:
        PToVLCSRFormFactorsTest() :
            TestCase("p_to_v_lcsr_form_factors_test")
        {
        }

        virtual void run() const
        {
            /* B -> K^* diagnostic values */
            {
                Parameters p = Parameters::Defaults();
                // Checking with unphysical mu for ease of comparison with Stefan's Mathematica's implementation
                p["mass::B_d"]                     = 5.2795;
                p["mass::K_d^*"]                   = 0.896;
                p["mass::K_u^*"]                   = 0.896;
                p["decay-constant::B_d"]           = 0.180;
                p["B->K^*::mu@Kstar-LCSR"]         = 1.0;
                p["B->K^*::M^2@Kstar-LCSR"]        = 50.;
                p["B->K^*::s_0^A1,0@Kstar-LCSR"]   = 70.;
                p["B->K^*::s_0^A1,1@Kstar-LCSR"]   = 0.0;
                p["B->K^*::s_0^A2,0@Kstar-LCSR"]   = 70.;
                p["B->K^*::s_0^A2,1@Kstar-LCSR"]   = 0.0;
                p["B->K^*::s_0^A30,0@Kstar-LCSR"]  = 70.;
                p["B->K^*::s_0^A30,1@Kstar-LCSR"]  = 0.0;
                p["B->K^*::s_0^V,0@Kstar-LCSR"]    = 70.;
                p["B->K^*::s_0^V,1@Kstar-LCSR"]    = 0.0;
                p["B->K^*::s_0^T1,0@Kstar-LCSR"]   = 70.;
                p["B->K^*::s_0^T1,1@Kstar-LCSR"]   = 0.0;
                p["B->K^*::s_0^T23A,0@Kstar-LCSR"] = 70.;
                p["B->K^*::s_0^T23A,1@Kstar-LCSR"] = 0.0;
                p["B->K^*::s_0^T23B,0@Kstar-LCSR"] = 70.;
                p["B->K^*::s_0^T23B,1@Kstar-LCSR"] = 0.0;
                // p["B->K^*::mu@Kstar-LCSR"]         = 3.0;
                // p["B->K^*::M^2@Kstar-LCSR"]        = 35.;
                // p["B->K^*::s_0^A1,0@Kstar-LCSR"]   = 30.;
                // p["B->K^*::s_0^A2,0@Kstar-LCSR"]   = 30.;
                // p["B->K^*::s_0^A30,0@Kstar-LCSR"]  = 30.;
                // p["B->K^*::s_0^V,0@Kstar-LCSR"]    = 30.;
                // p["B->K^*::s_0^T1,0@Kstar-LCSR"]   = 30.;
                // p["B->K^*::s_0^T23A,0@Kstar-LCSR"] = 30.;
                // p["B->K^*::s_0^T23B,0@Kstar-LCSR"] = 30.;

                Options o
                {
                    {"2pt-twist"_ok,         "2|3|4|5"}
                };

                AnalyticFormFactorPToVLCSR<BToKstar> ff{ p, o };
                auto diagnostics = ff.diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                const double eps = 1e-5;

                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair( 6.64457,    eps), // m_b(mu)
                    std::make_pair( 0.126841,   eps), // m_s(mu)
                    std::make_pair( 0.0106813,  eps), // m_ud(mu)
                    std::make_pair(-0.00208224, eps), // eta(u = 0.3, q2 = 1.0)
                    std::make_pair(-0.0138584,  eps), // eta'(u = 0.3, q2 = 1.0)
                    std::make_pair(-0.0458859,  eps), // eta''(u = 0.3, q2 = 1.0)
                    std::make_pair( 0.05458868, eps), // exp(-s(0.3, 1.0) / M2())
                    std::make_pair( 0.04,       eps), // a1perp(mu)
                    std::make_pair( 1.16172,    eps), // phi2perp(u = 0.3)
                    std::make_pair( 1.1088,     eps), // phi2para(u = 0.3)
                    std::make_pair( 0.82814,    eps), // phi3perp(u = 0.3)
                    std::make_pair( 0.514421,   eps), // phi3para(u = 0.3)
                    std::make_pair( 0.,         eps), // phi4perp(u = 0.3)
                    std::make_pair( 1.27397,    eps), // phi4para(u = 0.3)
                    std::make_pair( 0.,         eps), // psi4perp(u = 0.3)
                    std::make_pair( 1.24213,    eps), // psi4para(u = 0.3)
                    std::make_pair( 0.231876,   eps), // Barphi2perp(u = 0.3)
                    std::make_pair( 0.242989,   eps), // Barphi2para(u = 0.3)
                    std::make_pair( 0.0265496,  eps), // BarBarphi2perp(u = 0.3)
                    std::make_pair( 0.0289051,  eps), // BarBarphi2para(u = 0.3)
                    std::make_pair( 0.294102,   eps), // Barphi3perp(u = 0.3)
                    std::make_pair( 0.334127,   eps), // Barphi3para(u = 0.3)
                    std::make_pair( 0.0476302,  eps), // BarBarphi3perp(u = 0.3)
                    std::make_pair( 0.0655413,  eps), // BarBarphi3para(u = 0.3)
                    std::make_pair( 0.213218,   eps), // Barphi4para(u = 0.3)
                    std::make_pair( 0.,         eps), // Barpsi4perp(u = 0.3)
                    std::make_pair( 0.151422,   eps), // Barpsi4para(u = 0.3)
                    std::make_pair( 0.,         eps), // BarBarpsi4perp(u = 0.3)
                    std::make_pair( 0.00724092, eps), // BarBarpsi4para(u = 0.3)
                    std::make_pair( 0.,         eps), // Barphi5perp(u = 0.3)
                    std::make_pair( 3.8724,     eps), // I_perp[0](u = 0.3, q2 = 1.0)
                    std::make_pair( 43.5454,    eps), // I_perp[1](u = 0.3, q2 = 1.0)
                    std::make_pair( 3.8724,     eps), // I_T5_q[0](u = 0.3, q2 = 1.0)
                    std::make_pair( 43.5454,    eps), // I_T5_q[1](u = 0.3, q2 = 1.0)
                    std::make_pair( 3.8724,     eps), // I_T5_q[2](u = 0.3, q2 = 1.0)
                    std::make_pair( 43.5454,    eps), // I_T5_q[3](u = 0.3, q2 = 1.0)
                    std::make_pair( 43.5454,    eps), // D_I_T5_q[2][1](u = 0.3, q2 = 1.0)
                    std::make_pair( 43.5454,    eps), // D_I_T5_q[3][1](u = 0.3, q2 = 1.0)
                    std::make_pair( 43.5454,    eps), // D_I_T5_q[3][2](u = 0.3, q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // a_1(q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // a_2(q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // a_30(q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // v(q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // t_1(q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // t_23A(q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // t_23B(q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // a_1_dsr(q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // a_2_dsr(q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // a_30_dsr(q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // v_dsr(q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // t_1_dsr(q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // t_23A_dsr(q2 = 1.0)
                    std::make_pair( 4.123105,   eps), // t_23B_dsr(q2 = 1.0)
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_RELATIVE_ERROR(ff.v(-5.0),     0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.a_0(-5.0),   0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.a_1(-5.0),   0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.a_2(-5.0),   0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.a_12(-5.0),  0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.t_1(-5.0),   0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.t_2(-5.0),   0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.t_23(-5.0),  0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.t_3(-5.0),   0.849651, eps);

                TEST_CHECK_RELATIVE_ERROR(ff.v_dsr(-5.0),     0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.a_1_dsr(-5.0),   0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.a_2_dsr(-5.0),   0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.a_30_dsr(-5.0),  0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.t_1_dsr(-5.0),   0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.t_23A_dsr(-5.0), 0.849651, eps);
                TEST_CHECK_RELATIVE_ERROR(ff.t_23B_dsr(-5.0), 0.849651, eps);
            }
        }
} p_to_v_form_factors_test;
