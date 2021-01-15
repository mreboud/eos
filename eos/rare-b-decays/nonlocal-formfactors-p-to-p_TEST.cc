/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#include <test/test.hh>
#include <eos/rare-b-decays/nonlocal-formfactors.hh>

using namespace test;
using namespace eos;

// class NonlocalFormFactorLCSRTest :
//     public TestCase
// {
//     public:
//         NonlocalFormFactorLCSRTest() :
//             TestCase("nonlocal_formfactor_lcsr_test")
//         {
//         }

//         virtual void run() const
//         {
//             static const double eps = 1e-5;

//             /* test cases*/
//             {
//                 Parameters p = Parameters::Defaults();
//                 p["B::1/lambda_B_p"]          = 2.173913;
//                 p["B::lambda_E^2"]            = 0.03;
//                 p["B::lambda_H^2"]            = 0.06;
//                 p["mass::B_d"]                = 5.27958;
//                 p["mass::K_d"]                = 0.497614;
//                 p["decay-constant::B_d"]      = 0.1905;
//                 p["decay-constant::K_d"]      = 0.1561;
//                 p["B->K::M^2@B-LCSR"]         = 1.0;
//                 p["B->K::s_0^+,0@B-LCSR"]     = 1.05;
//                 p["B->K::s_0^+,1@B-LCSR"]     = 0.0;
//                 p["B->K::mu@B-LCSR"]          = 1.0;
//                 p["b->sccbar::mu"]            = 1.0;
//                 p["b->sccbar::mu_c"]          = 1.0;
//                 // C_1_AK = C_2_EOS - ... C_1_EOS; setting c1 -> 0, c2 -> C_1_AK for this test-case only
//                 p["b->s::c1"]                 = 0.0;
//                 p["b->s::c2"]                 = 1.05873559;

//                 Options o = { { "model", "WilsonScan" } };

//                 auto nc = NonlocalFormFactor<nc::PToP>::make("B->K::LCSR", p, o);
//                 auto diagnostics = nc->diagnostics();

//                 std::cout << "Diagnostics:" << std::endl;
//                 for (auto & d : diagnostics)
//                 {
//                     std::cout << d.description << ": " << d.value << std::endl;
//                 }
//                 std::cout << "Diagnostics ended" << std::endl;

//                 static const std::vector<std::pair<double, double>> reference
//                 {
//                     /* quark masses in the MSbar scheme */
//                     std::make_pair(0.1268577,  1.0e-7),            //m_v(mu) in the MSbar scheme
//                     std::make_pair(1.4482600,  1.0e-6),            //m_c(mu) in the MSbar scheme
//                     std::make_pair(0.1561,     1.0e-6),            //final state decay constant

//                     /* I1_A_phi_3 */
//                     std::make_pair(-1.76201e-2,  1.0e-6 ),         // I1_A_phi_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-5.62386e-3,  1.0e-7 ),         // I1_A_phi_3(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair( 2.39547e-4,  1.0e-8 ),         // I1_A_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 6.18902e-5,  1.0e-9 ),         // I1_A_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-6.63902e-3,  1.0e-7 ),         // I2_A_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.45856e-3,  1.0e-7 ),         // I2_A_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair( 2.19348e-6,  1.0e-10),         // I1_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 9.56343e-8,  1.0e-12),         // I1_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-1.23581e-5,  1.0e-9 ),         // I2_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 2.09778e-6,  1.0e-10),         // I2_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-6.61329e-4,  1.0e-8 ),         // I3_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.01596e-4,  1.0e-4 ),         // I3_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-9.45221e-2,  1.0e-6 ),         // I3d1_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-2.98012e-2,  1.0e-6 ),         // I3d1_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* I1_A_phi_4 */
//                     std::make_pair(-4.04613e-2,  1.0e-6 ),         // I1_A_phi_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-3.02563e-2,  1.0e-6 ),         // I1_A_phi_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair( 5.80632e-3,  1.0e-7 ),         // I1_A_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 2.97265e-3,  1.0e-6 ),         // I1_A_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-7.30190e-2,  1.0e-6 ),         // I2_A_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-5.10359e-2,  1.0e-6 ),         // I2_A_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair( 9.06302e-5,  1.0e-9 ),         // I1_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 7.81572e-6,  1.0e-10),         // I1_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-5.35228e-4,  1.0e-8 ),         // I2_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 1.63031e-4,  1.0e-8 ),         // I2_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-2.69628e-2,  1.0e-6 ),         // I3_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-8.07360e-3,  1.0e-7 ),         // I3_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-2.55750e-0,  1.0e-4 ),         // I3d1_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.57621e-0,  1.0e-4 ),         // I3d1_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* I1_A_psi_4 */
//                     std::make_pair(-1.98959e-4,  1.0e-8 ),         // I1_A_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.42146e-5,  1.0e-9 ),         // I1_A_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-5.64130e-5,  1.0e-9 ),         // I2_A_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 8.44926e-5,  1.0e-9 ),         // I2_A_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair(-3.56595e-8,  1.0e-12),         // I1_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.40931e-9,  1.0e-13),         // I1_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 8.95686e-7,  1.0e-11),         // I2_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 3.82319e-8,  1.0e-12),         // I2_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 5.34015e-7,  1.0e-11),         // I3_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-3.88418e-7,  1.0e-11),         // I3_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 7.63683e-5,  1.0e-9 ),         // I3d1_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.15223e-4,  1.0e-8 ),         // I3d1_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* I1_A_chi_4 */
//                     std::make_pair(-4.25340e-4,  1.0e-8 ),         // I1_A_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-8.70429e-5,  1.0e-9 ),         // I1_A_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 8.49573e-3,  1.0e-7 ),         // I2_A_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 1.57675e-3,  1.0e-7 ),         // I2_A_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair(-1.01613e-5,  1.0e-9 ),         // I1_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-2.54200e-7,  1.0e-11),         // I1_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 5.86388e-5,  1.0e-9 ),         // I2_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-5.43771e-6,  1.0e-10),         // I2_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 3.04318e-3,  1.0e-7 ),         // I3_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 2.66276e-4,  1.0e-8 ),         // I3_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 4.34954e-1,  1.0e-5 ),         // I3d1_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 7.81042e-2,  1.0e-6 ),         // I3d1_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                 };
//                 TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus( 0.0)),  0.0,        1.0e-11);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus( 0.0)),  0.0,        1.0e-11);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus( 1.0)), -7.58964e-9, 1.0e-11);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus( 1.0)),  0.0,        1.0e-11);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus( 3.0)), -7.58213e-8, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus( 3.0)),  0.0,        1.0e-11);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus(-5.0)), -5.20521e-8, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus(-5.0)),  0.0,        1.0e-11);

//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_A( 0.0)),  0.734885, 1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_A( 0.0)),  0.0,      1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_A( 1.0)),  0.526179, 1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_A( 1.0)),  0.0,      1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_A( 3.0)),  0.468461, 1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_A( 3.0)),  0.0,      1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_A(-5.0)),  0.394032, 1.0e-3);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_A(-5.0)),  0.0,      1.0e-4);
//             }
//         }
// } nonlocal_formfactor_lcsr_test;

class NonlocalFormFactorGvDV2020Test :
    public TestCase
{
    public:
        NonlocalFormFactorGvDV2020Test() :
            TestCase("nonlocal_formfactor_GvDV2020_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"]                               = 5.27942;
                p["mass::K_d"]                               = 0.4936;
                p["mass::J/psi"]                             = 3.0969;
                p["mass::psi(2S)"]                           = 3.6860;
                p["mass::D^0"]                               = 1.86723;
                p["b->sccbar::t_0"]                          = 9.0;
                p["b->sccbar::t_s"]                          = -17.4724;
                p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;
                p["B->Kccbar::Re{alpha_0^plus}@GvDV2020"]    = 2.0;
                p["B->Kccbar::Im{alpha_0^plus}@GvDV2020"]    = 3.0;
                p["B->Kccbar::Re{alpha_1^plus}@GvDV2020"]    = 4.0;
                p["B->Kccbar::Im{alpha_1^plus}@GvDV2020"]    = 5.0;
                p["B->Kccbar::Re{alpha_2^plus}@GvDV2020"]    = 6.0;
                p["B->Kccbar::Im{alpha_2^plus}@GvDV2020"]    = 7.0;

                Options o = { { "model", "WilsonScan" } };

                auto nc = NonlocalFormFactor<nc::PToP>::make("B->K::GvDV2020", p, o);


                auto diagnostics = nc->diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {
                    /* outer functions */
                    std::make_pair( 0.0,      eps),            // Re{1/phi_+(q2 = 0.0)}
                    std::make_pair( 0.0,      eps),            // Im{1/phi_+(q2 = 0.0)}
                    std::make_pair( -16.3399, eps),            // Re{phi_+(q2 = 16.0)}
                    std::make_pair( 2.06922,  eps),            // Im{phi_+(q2 = 16.0)}

                    std::make_pair( 5.73999,  eps),            // Re{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, 2.0, 3.0, 4.0)}
                    std::make_pair( 0.0,      eps),            // Im{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, 2.0, 3.0, 4.0)}


                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus(-1.0)),  0.131447,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus(-1.0)),  0.156871,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus(0.0)),   0.,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus(0.0)),   0.,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus(4.0)),  -1.09013,    eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus(4.0)),  -1.2906,     eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus(12.0)), 14.5673,  10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus(12.0)), 16.9699,  10*eps);

                TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus_residue_jpsi()),  -1.09287,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus_residue_jpsi()),  -1.2799,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus_residue_psi2s()),  3.10823,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus_residue_psi2s()),  3.60742,  eps);
            }
        }
} nonlocal_formfactor_gvdv2020_test;


class NonlocalFormFactorGvDV2021Test :
    public TestCase
{
    public:
        NonlocalFormFactorGvDV2021Test() :
            TestCase("nonlocal_formfactor_GvDV2021_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"]                                = 5.27942;
                p["mass::K_d"]                                = 0.4936;
                p["mass::J/psi"]                              = 3.0969;
                p["mass::psi(2S)"]                            = 3.6860;
                p["mass::B_s^*"]                              = 5.4154;
                p["mass::D^0"]                                = 1.86723;
                p["b->sccbar::t_0"]                           = 9.0;
                p["b->sccbar::t_s"]                           = -17.4724;
                p["b->sccbar::chiOPE@GRvDV2021"]              = 1.81e-4;
                p["B->Kccbar::Re{alpha_0^plus}@GRvDV2021"]    = 2.0;
                p["B->Kccbar::Im{alpha_0^plus}@GRvDV2021"]    = 3.0;
                p["B->Kccbar::Re{alpha_1^plus}@GRvDV2021"]    = 4.0;
                p["B->Kccbar::Im{alpha_1^plus}@GRvDV2021"]    = 5.0;
                p["B->Kccbar::Re{alpha_2^plus}@GRvDV2021"]    = 6.0;
                p["B->Kccbar::Im{alpha_2^plus}@GRvDV2021"]    = 7.0;

                Options o = { { "model", "WilsonScan" } };

                auto nc = NonlocalFormFactor<nc::PToP>::make("B->K::GRvDV2021", p, o);


                auto diagnostics = nc->diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair( 5.31175,  eps),            // Re{phi_+(q2 = 16.0)}
                    std::make_pair( -3.8452,  eps),            // Im{phi_+(q2 = 16.0)}

                    std::make_pair( 1.16924,  eps),            // Re{P(q2 = 1.0, 2.0, 3.0, 4.0)}
                    std::make_pair( 0.0,      eps),            // Im{P(q2 = 1.0, 2.0, 3.0, 4.0)}
                    std::make_pair( 1.16924,  eps),            // Re{P(q2 = 1.0, (2.0,5.0), (3.0,6.0), (4.0,7.0))}
                    std::make_pair( 2.71519,  eps),            // Im{P(q2 = 1.0, (2.0,5.0), (3.0,6.0), (4.0,7.0))}

                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);


                TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus(-1.0)), -0.071174,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus(-1.0)), -0.0983567,  eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus(0.0)),   0.,         eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus(0.0)),   0.,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus(4.0)),   0.410774,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus(4.0)),   0.58287,    eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus(12.0)), -1.18903,    eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus(12.0)), -1.88908,    eps);

                TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus_residue_jpsi()),   0.180334,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus_residue_jpsi()),   0.273359,  eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_plus_residue_psi2s()), -0.155624,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_plus_residue_psi2s()), -0.225437,  eps);
            }
        }
} nonlocal_formfactor_grvdv2021_test;

