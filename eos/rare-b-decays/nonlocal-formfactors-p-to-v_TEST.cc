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

// class NonlocalFormFactorTest :
//     public TestCase
// {
//     public:
//         NonlocalFormFactorTest() :
//             TestCase("nonlocal_formfactor_test")
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
//                 p["mass::K_d^*"]              = 0.89594;
//                 p["decay-constant::B_d"]      = 0.1905;
//                 p["decay-constant::K_d^*"]    = 0.204;
//                 p["B->K^*::M^2@B-LCSR"]       = 1.0;
//                 p["B->K^*::s_0^V,0@B-LCSR"]   = 1.7;
//                 p["B->K^*::s_0^V,1@B-LCSR"]   = 0.0;
//                 p["B->K^*::s_0^A1,0@B-LCSR"]  = 1.7;
//                 p["B->K^*::s_0^A1,1@B-LCSR"]  = 0.0;
//                 p["B->K^*::s_0^A2,0@B-LCSR"]  = 1.7;
//                 p["B->K^*::s_0^A2,1@B-LCSR"]  = 0.0;
//                 p["B->K^*::mu@B-LCSR"]        = 1.0;
//                 p["b->sccbar::mu"]            = 1.0;
//                 p["b->sccbar::mu_c"]          = 1.0;
//                 // C_1_AK = C_2_EOS - ... C_1_EOS; setting c1 -> 0, c2 -> C_1_AK for this test-case only
//                 p["b->s::c1"]                 = 0.0;
//                 p["b->s::c2"]                 = 1.05873559;

//                 Options o = { { "model", "WilsonScan" } };

//                 auto nc = NonlocalFormFactor<nc::PToV>::make("B->K^*::LCSR", p, o);
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
//                     std::make_pair(0.204,      1.0e-6),            //final state decay constant

//                     /* V1 */

//                     /* I1_A_phi_3 */
//                     std::make_pair(-1.72092e-4,  1.0e-8 ),         // I1_A_phi_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 7.60581e-4,  1.0e-8 ),         // I1_A_phi_3(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair( 4.21716e-5,  1.0e-9 ),         // I1_A_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 2.39269e-6,  1.0e-10),         // I1_A_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 4.17710e-5,  1.0e-9 ),         // I2_A_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 3.29897e-5,  1.0e-9 ),         // I2_A_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair(-7.17168e-7,  1.0e-11),         // I1_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-4.77689e-9,  1.0e-13),         // I1_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 1.51033e-5,  1.0e-9 ),         // I2_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 1.21653e-7,  1.0e-11),         // I2_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 9.40093e-6,  1.0e-10),         // I3_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 7.21846e-8,  1.0e-12),         // I3_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 1.35323e-3,  1.0e-7 ),         // I3d1_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 2.122061e-5,  1.0e-9 ),         // I3d1_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* I1_A_phi_4 */
//                     std::make_pair( 4.83417e-2,  1.0e-6 ),         // I1_A_phi_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-3.22366e-3,  1.0e-7 ),         // I1_A_phi_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair(-1.43604e-3,  1.0e-7 ),         // I1_A_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.30395e-4,  1.0e-8 ),         // I1_A_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 2.88257e-2,  1.0e-6 ),         // I2_A_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.01379e-3,  1.0e-7 ),         // I2_A_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair( 4.77695e-6,  1.0e-10),         // I1_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 4.69068e-7,  1.0e-11),         // I1_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-5.05019e-5,  1.0e-9 ),         // I2_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 3.90288e-5,  1.0e-9 ),         // I2_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 3.92955e-4,  1.0e-8 ),         // I3_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.91119e-5,  1.0e-9 ),         // I3_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 3.76628e-2,  1.0e-6 ),         // I3d1_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-3.80580e-3,  1.0e-7 ),         // I3d1_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* I1_A_psi_4 */
//                     std::make_pair( 5.28700e-5,  1.0e-9 ),         // I1_A_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 2.73389e-6,  1.0e-10),         // I1_A_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-2.96432e-5,  1.0e-9 ),         // I2_A_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 4.44976e-5,  1.0e-9 ),         // I2_A_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair(-9.59579e-7,  1.0e-11),         // I1_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-6.99541e-9,  1.0e-13),         // I1_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 1.87943e-5,  1.0e-9 ),         // I2_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-2.40921e-7,  1.0e-11),         // I2_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-2.80607e-7,  1.0e-11),         // I3_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 2.04558e-7,  1.0e-11),         // I3_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-4.00993e-5,  1.0e-9 ),         // I3d1_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 6.06440e-5,  1.0e-9 ),         // I3d1_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* I1_A_chi_4 */
//                     std::make_pair(-9.41705e-5,  1.0e-9 ),         // I1_A_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-8.93572e-7,  1.0e-11),         // I1_A_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-1.25716e-4,  1.0e-8 ),         // I2_A_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 1.13976e-6,  1.0e-10),         // I2_A_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair( 1.37981e-6,  1.0e-10),         // I1_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.43442e-9,  1.0e-13),         // I1_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-3.18867e-5,  1.0e-9 ),         // I2_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-8.01616e-7,  1.0e-11),         // I2_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-4.38055e-5,  1.0e-9 ),         // I3_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 2.19374e-7,  1.0e-11),         // I3_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-6.30505e-3,  1.0e-7 ),         // I3d1_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 6.55081e-5,  1.0e-9 ),         // I3d1_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* V2 */

//                     /* I1_A_phi_3 */
//                     std::make_pair(-2.14517e-3,  1.0e-7 ),         // I1_A_phi_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 9.18189e-3,  1.0e-7 ),         // I1_A_phi_3(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair(-5.92753e-4,  1.0e-8 ),         // I1_A_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.04083e-5,  1.0e-9 ),         // I1_A_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 5.33252e-4,  1.0e-8 ),         // I2_A_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 3.86336e-4,  1.0e-8 ),         // I2_A_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair( 2.21919e-7,  1.0e-11),         // I1_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-3.30391e-10, 1.0e-14),         // I1_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-1.11523e-5,  1.0e-9 ),         // I2_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-2.38714e-7,  1.0e-11),         // I2_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 1.18335e-4,  1.0e-8 ),         // I3_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 1.07695e-6,  1.0e-10),         // I3_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 1.68933e-2,  1.0e-6 ),         // I3d1_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 3.15929e-4,  1.0e-8 ),         // I3d1_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* I1_A_phi_4 */
//                     std::make_pair( 6.08131e-1,  1.0e-5 ),         // I1_A_phi_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-2.64246e-2,  1.0e-6 ),         // I1_A_phi_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair(-3.46768e-2,  1.0e-6 ),         // I1_A_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-9.97687e-4,  1.0e-8 ),         // I1_A_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 3.62620e-1,  1.0e-5 ),         // I2_A_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-8.86237e-3,  1.0e-7 ),         // I2_A_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair(-1.18144e-6,  1.0e-10),         // I1_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-6.13453e-7,  1.0e-11),         // I1_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-3.56691e-4,  1.0e-8 ),         // I2_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.82647e-5,  1.0e-9 ),         // I2_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 4.94026e-3,  1.0e-7 ),         // I3_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-2.03064e-4,  1.0e-8 ),         // I3_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 4.67634e-1,  1.0e-5 ),         // I3d1_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-4.01671e-2,  1.0e-6 ),         // I3d1_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* I1_A_psi_4 */
//                     std::make_pair(-6.62324e-4,  1.0e-8 ),         // I1_A_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 4.34163e-5,  1.0e-9 ),         // I1_A_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-3.55015e-4,  1.0e-8 ),         // I2_A_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 5.17099e-4,  1.0e-8 ),         // I2_A_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair( 2.88554e-7,  1.0e-11),         // I1_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 4.82621e-9,  1.0e-13),         // I1_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-2.75707e-6,  1.0e-10),         // I2_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-6.71116e-9,  1.0e-13),         // I2_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-3.36063e-6,  1.0e-10),         // I3_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 2.37714e-6,  1.0e-10),         // I3_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-4.76041e-4,  1.0e-8 ),         // I3d1_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 7.01571e-4,  1.0e-8 ),         // I3d1_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* I1_A_chi_4 */
//                     std::make_pair( 1.25683e-3,  1.0e-7 ),         // I1_A_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 1.11500e-4,  1.0e-8 ),         // I1_A_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-1.58150e-3,  1.0e-7 ),         // I2_A_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 9.34275e-6,  1.0e-10),         // I2_A_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair(-4.43718e-7,  1.0e-11),         // I1_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 1.05209e-8,  1.0e-12),         // I1_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 4.57865e-5,  1.0e-9 ),         // I2_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 6.14055e-7,  1.0e-11),         // I2_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-5.51060e-4,  1.0e-8 ),         // I3_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 1.92345e-6,  1.0e-10),         // I3_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-7.86615e-2,  1.0e-6 ),         // I3d1_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 5.72701e-4,  1.0e-8 ),         // I3d1_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* V23 */

//                     /* I1_A_phi_3 */
//                     std::make_pair(-3.13913e-4,  1.0e-8 ),         // I1_A_phi_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-3.34136e-5,  1.0e-9 ),         // I1_A_phi_3(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair( 2.88167e-5,  1.0e-9 ),         // I1_A_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.90882e-6,  1.0e-10),         // I1_A_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-1.01339e-4,  1.0e-8 ),         // I2_A_phi_bar_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.12568e-5,  1.0e-9 ),         // I2_A_phi_bar_3(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair( 3.32078e-7,  1.0e-11),         // I1_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 4.38301e-8,  1.0e-12),         // I1_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-5.94698e-6,  1.0e-10),         // I2_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-9.76472e-7,  1.0e-11),         // I2_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-6.00794e-6,  1.0e-10),         // I3_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-5.14672e-7,  1.0e-11),         // I3_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-1.15534e-3,  1.0e-7 ),         // I3d1_A_phi_double_3(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.99873e-4,  1.0e-8 ),         // I3d1_A_phi_double_3(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* I1_A_phi_4 */
//                     std::make_pair(-2.13289e-3,  1.0e-7 ),         // I1_A_phi_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-2.42022e-3,  1.0e-7 ),         // I1_A_phi_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair( 7.29158e-4,  1.0e-8 ),         // I1_A_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.76269e-4,  1.0e-8 ),         // I1_A_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-1.68304e-3,  1.0e-7 ),         // I2_A_phi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-8.61187e-4,  1.0e-8 ),         // I2_A_phi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair( 1.44698e-5,  1.0e-9 ),         // I1_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 3.58227e-6,  1.0e-10),         // I1_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-2.34335e-4,  1.0e-8 ),         // I2_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-7.69839e-5,  1.0e-9 ),         // I2_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-2.44378e-4,  1.0e-8 ),         // I3_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-4.16472e-5,  1.0e-9 ),         // I3_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-3.52460e-2,  1.0e-6 ),         // I3d1_A_phi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.20905e-2,  1.0e-6 ),         // I3d1_A_phi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* I1_A_psi_4 */
//                     std::make_pair(-4.53872e-5,  1.0e-9 ),         // I1_A_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-7.26398e-6,  1.0e-10),         // I1_A_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-1.18673e-6,  1.0e-10),         // I2_A_psi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 9.08785e-7,  1.0e-11),         // I2_A_psi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair(-2.65409e-8,  1.0e-12),         // I1_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-6.47972e-10, 1.0e-14),         // I1_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-2.24612e-7,  1.0e-11),         // I2_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-8.77970e-9,  1.0e-13),         // I2_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-1.12338e-8,  1.0e-12),         // I3_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 4.17774e-9,  1.0e-13),         // I3_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair(-2.16183e-6,  1.0e-10),         // I3d1_A_psi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 1.63542e-6,  1.0e-10),         // I3d1_A_psi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     /* I1_A_chi_4 */
//                     std::make_pair(-2.05466e-4,  1.0e-8 ),         // I1_A_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-8.56729e-6,  1.0e-10),         // I1_A_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 2.45291e-4,  1.0e-8 ),         // I2_A_chi_bar_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 1.86762e-5,  1.0e-9 ),         // I2_A_chi_bar_4(0.01, 0.7, 0.9, 0.1, 3.0)

//                     std::make_pair(-1.58064e-6,  1.0e-10),         // I1_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair(-1.16506e-7,  1.0e-11),         // I1_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 2.69069e-5,  1.0e-9 ),         // I2_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 2.54917e-6,  1.0e-10),         // I2_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 2.76140e-5,  1.0e-9 ),         // I3_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 1.36121e-6,  1.0e-10),         // I3_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                     std::make_pair( 5.31024e-3,  1.0e-5 ),         // I3d1_A_chi_double_4(0.02, 0.4, 0.3, 0.7, 2.0)
//                     std::make_pair( 5.28651e-4,  1.0e-6 ),         // I3d1_A_chi_double_4(0.01, 0.7, 0.9, 0.1, 3.0)
//                 };
//                 TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp( 0.0)), -5.33392e-8, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp( 0.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp( 3.0)), -8.45168e-8, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp( 3.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp(-5.0)), -6.91408e-8, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp(-5.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp(-7.0)), -7.99300e-8, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp(-7.0)),  0.0,        1.0e-10);

//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_para( 0.0)), -1.38926e-7, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para( 0.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_para( 3.0)),  2.83843e-8, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para( 3.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_para(-5.0)), -2.45818e-7, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para(-5.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_para(-7.0)), -2.64257e-7, 1.0e-9 );
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para(-7.0)),  0.0,        1.0e-10);

//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_long( 0.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long( 0.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_long( 3.0)),  8.74998e-9, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long( 3.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_long(-5.0)), -1.28370e-7, 1.0e-9 );
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long(-5.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_long(-7.0)), -1.92647e-7, 1.0e-9 );
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long(-7.0)),  0.0,        1.0e-11);

//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V1( 0.0)),  0.590215, 1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V1( 0.0)),  0.0,      1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V1( 1.0)),  0.730058, 1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V1( 1.0)),  0.0,      1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V1( 3.0)),  1.022060, 1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V1( 3.0)),  0.0,      1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V1(-5.0)),  0.393984, 1.0e-3);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V1(-5.0)),  0.0,      1.0e-4);

//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V2( 0.0)),  0.800823, 1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V2( 0.0)),  0.0,      1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V2( 1.0)),  0.808517, 1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V2( 1.0)),  0.0,      1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V2( 3.0)),  0.747949, 1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V2( 3.0)),  0.0,      1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V2(-5.0)),  0.786164, 1.0e-3);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V2(-5.0)),  0.0,      1.0e-4);

//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V23( 0.0)),  0.892919, 1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V23( 0.0)),  0.0,      1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V23( 1.0)),  0.881786, 1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V23( 1.0)),  0.0,      1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V23( 3.0)),  0.671125, 1.0e-3);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V23( 3.0)),  0.0,      1.0e-4);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V23(-5.0)),  0.910152, 1.0e-3);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V23(-5.0)),  0.0,      1.0e-4);
//             }

//             {
//                 Parameters p = Parameters::Defaults();
//                 p["B_s::1/lambda_B_p"]          = 1.92308;
//                 p["B_s::lambda_E^2"]            = 0.03;
//                 p["B_s::lambda_H^2"]            = 0.06;
//                 p["mass::B_s"]                  = 5.36677;
//                 p["mass::phi"]                  = 1.019461;
//                 p["decay-constant::B_s"]        = 0.2307;
//                 p["decay-constant::phi"]        = 0.233;
//                 p["B_s->phi::M^2@B-LCSR"]       = 1.0;
//                 p["B_s->phi::s_0^V,0@B-LCSR"]   = 1.7;
//                 p["B_s->phi::s_0^V,1@B-LCSR"]   = 0.0;
//                 p["B_s->phi::s_0^A1,0@B-LCSR"]  = 1.7;
//                 p["B_s->phi::s_0^A1,1@B-LCSR"]  = 0.0;
//                 p["B_s->phi::s_0^A2,0@B-LCSR"]  = 1.7;
//                 p["B_s->phi::s_0^A2,1@B-LCSR"]  = 0.0;
//                 p["B_s->phi::mu@B-LCSR"]        = 1.0;
//                 p["b->sccbar::mu"]              = 1.0;
//                 p["b->sccbar::mu_c"]            = 1.0;
//                 // C_1_AK = C_2_EOS - ... C_1_EOS; setting c1 -> 0, c2 -> C_1_AK for this test-case only
//                 p["b->s::c1"]                   = 0.0;
//                 p["b->s::c2"]                   = 1.05873559;

//                 Options o = {
//                   { "model", "WilsonScan" },
//                   { "q",    "s"  }
//                 };

//                 auto nc = NonlocalFormFactor<nc::PToV>::make("B_s->phi::LCSR", p, o);

//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp( 0.0)), -4.90182e-8, 1.0e-11);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp( 0.0)),  0.0,        1.0e-11);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp( 3.0)), -8.55393e-8, 1.0e-11);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp( 3.0)),  0.0,        1.0e-11);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp(-5.0)), -5.90284e-8, 1.0e-11);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp(-5.0)),  0.0,        1.0e-11);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp(-7.0)), -6.84565e-8, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp(-7.0)),  0.0,        1.0e-11);

//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_para( 0.0)), -1.42157e-7, 1.0e-9 );
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para( 0.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_para( 3.0)),  5.86136e-9, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para( 3.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_para(-5.0)), -2.37766e-7, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para(-5.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_para(-7.0)), -2.54312e-7, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para(-7.0)),  0.0,        1.0e-10);

//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_long( 0.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long( 0.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_long( 3.0)),  9.2323e-10, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long( 3.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_long(-5.0)), -1.04636e-7, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long(-5.0)),  0.0,        1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(real(nc->H_long(-7.0)), -1.59028e-7, 1.0e-10);
//                 TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long(-7.0)),  0.0,        1.0e-11);

//                TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V1( 0.0)),  0.546004, 1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V1( 0.0)),  0.0,      1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V1( 3.0)),  1.004730, 1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V1( 3.0)),  0.0,      1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V1(-5.0)),  0.299404, 1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V1(-5.0)),  0.0,      1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V1(-7.0)),  0.333725, 1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V1(-7.0)),  0.0,      1.0e-3);

//                TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V2( 0.0)),  0.801087, 1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V2( 0.0)),  0.0,      1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V2( 3.0)),  0.325204, 1.0e-1);
//                TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V2( 3.0)),  0.0,      1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V2(-5.0)),  0.782896, 1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V2(-5.0)),  0.0,      1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V2(-7.0)),  0.779503, 1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V2(-7.0)),  0.0,      1.0e-3);

//                TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V23( 0.0)),  0.895153, 1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V23( 0.0)),  0.0,      1.0e-3);
// //               TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V23( 3.0)),  0.881786, 1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V23( 3.0)),  0.0,      1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V23(-5.0)),  0.915112, 1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V23(-5.0)),  0.0,      1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(real(nc->normalized_moment_V23(-7.0)),  0.918041, 1.0e-3);
//                TEST_CHECK_NEARLY_EQUAL(imag(nc->normalized_moment_V23(-7.0)),  0.0,      1.0e-3);
//             }
//         }
// } nonlocal_formfactor_test;


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
                p["mass::K_d^*"]                             = 0.89555;
                p["mass::J/psi"]                             = 3.0969;
                p["mass::psi(2S)"]                           = 3.6860;
                p["mass::D^0"]                               = 1.86723;
                p["b->sccbar::t_0"]                          = 9.0;
                p["b->sccbar::t_s"]                          = -17.4724;
                p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;
                p["B->K^*ccbar::Re{alpha_0^perp}@GvDV2020"]  = 2.0;
                p["B->K^*ccbar::Im{alpha_0^perp}@GvDV2020"]  = 3.0;
                p["B->K^*ccbar::Re{alpha_1^perp}@GvDV2020"]  = 4.0;
                p["B->K^*ccbar::Im{alpha_1^perp}@GvDV2020"]  = 5.0;
                p["B->K^*ccbar::Re{alpha_2^perp}@GvDV2020"]  = 6.0;
                p["B->K^*ccbar::Im{alpha_2^perp}@GvDV2020"]  = 7.0;
                p["B->K^*ccbar::Re{alpha_0^para}@GvDV2020"]  = 8.0;
                p["B->K^*ccbar::Im{alpha_0^para}@GvDV2020"]  = 9.0;
                p["B->K^*ccbar::Re{alpha_1^para}@GvDV2020"]  = 10.0;
                p["B->K^*ccbar::Im{alpha_1^para}@GvDV2020"]  = 11.0;
                p["B->K^*ccbar::Re{alpha_2^para}@GvDV2020"]  = 12.0;
                p["B->K^*ccbar::Im{alpha_2^para}@GvDV2020"]  = 13.0;
                p["B->K^*ccbar::Re{alpha_0^long}@GvDV2020"]  = 14.0;
                p["B->K^*ccbar::Im{alpha_0^long}@GvDV2020"]  = 15.0;
                p["B->K^*ccbar::Re{alpha_1^long}@GvDV2020"]  = 16.0;
                p["B->K^*ccbar::Im{alpha_1^long}@GvDV2020"]  = 17.0;
                p["B->K^*ccbar::Re{alpha_2^long}@GvDV2020"]  = 18.0;
                p["B->K^*ccbar::Im{alpha_2^long}@GvDV2020"]  = 19.0;

                Options o = { { "model", "WilsonScan" } };

                auto nc = NonlocalFormFactor<nc::PToV>::make("B->K^*::GvDV2020", p, o);


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
                    std::make_pair( 0.0,      eps),            // Re{1/phi_long(q2 = 0.0)}
                    std::make_pair( 0.0,      eps),            // Im{1/phi_long(q2 = 0.0)}

                    std::make_pair( -36.5755, 10*eps),         // Re{phi_long(q2 = 16.0)}
                    std::make_pair( 4.63177,  10*eps),         // Im{phi_long(q2 = 16.0)}

                    std::make_pair( 24.6148,  10*eps),         // Re{phi_perp(q2 = 16.0)}
                    std::make_pair( -13.2048, 10*eps)          // Im{phi_perp(q2 = 16.0)}
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp(16.0)),  -2.36353,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp(16.0)),  -1.27642,     eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_para(16.0)),  -4.48563,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para(16.0)),  -2.2198,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_long(16.0)),   5.53271,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long(16.0)),   0.443831,    eps);

                TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp_residue_jpsi()),   3.00101,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp_residue_jpsi()),   3.50296,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp_residue_psi2s()), -3.57057,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp_residue_psi2s()), -4.14177,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_para_residue_jpsi()),   6.01269,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para_residue_jpsi()),   6.51464,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_para_residue_psi2s()), -6.99773,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para_residue_psi2s()), -7.56893,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_long_residue_jpsi()),  -2.81615,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long_residue_jpsi()),  -2.97279,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_long_residue_psi2s()),  6.19006,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long_residue_psi2s()),  6.52922,   eps);
            }
        }
} nonlocal_formfactor_gvdv2020_test;


class NonlocalFormFactorGRvDV2021Test :
    public TestCase
{
    public:
        NonlocalFormFactorGRvDV2021Test() :
            TestCase("nonlocal_formfactor_GRvDV2021_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"]                                = 5.27942;
                p["mass::K_d^*"]                              = 0.89555;
                p["mass::J/psi"]                              = 3.0969;
                p["mass::psi(2S)"]                            = 3.6860;
                p["mass::B_s^*"]                              = 5.4154;
                p["mass::D^0"]                                = 1.86723;
                p["b->sccbar::t_0"]                           = 9.0;
                p["b->sccbar::t_s"]                           = -17.4724;
                p["b->sccbar::chiOPE@GRvDV2021"]              = 1.81e-4;
                p["B->K^*ccbar::Re{alpha_0^perp}@GRvDV2021"]  = 2.0;
                p["B->K^*ccbar::Im{alpha_0^perp}@GRvDV2021"]  = 3.0;
                p["B->K^*ccbar::Re{alpha_1^perp}@GRvDV2021"]  = 4.0;
                p["B->K^*ccbar::Im{alpha_1^perp}@GRvDV2021"]  = 5.0;
                p["B->K^*ccbar::Re{alpha_2^perp}@GRvDV2021"]  = 6.0;
                p["B->K^*ccbar::Im{alpha_2^perp}@GRvDV2021"]  = 7.0;
                p["B->K^*ccbar::Re{alpha_0^para}@GRvDV2021"]  = 8.0;
                p["B->K^*ccbar::Im{alpha_0^para}@GRvDV2021"]  = 9.0;
                p["B->K^*ccbar::Re{alpha_1^para}@GRvDV2021"]  = 10.0;
                p["B->K^*ccbar::Im{alpha_1^para}@GRvDV2021"]  = 11.0;
                p["B->K^*ccbar::Re{alpha_2^para}@GRvDV2021"]  = 12.0;
                p["B->K^*ccbar::Im{alpha_2^para}@GRvDV2021"]  = 13.0;
                p["B->K^*ccbar::Re{alpha_0^long}@GRvDV2021"]  = 14.0;
                p["B->K^*ccbar::Im{alpha_0^long}@GRvDV2021"]  = 15.0;
                p["B->K^*ccbar::Re{alpha_1^long}@GRvDV2021"]  = 16.0;
                p["B->K^*ccbar::Im{alpha_1^long}@GRvDV2021"]  = 17.0;
                p["B->K^*ccbar::Re{alpha_2^long}@GRvDV2021"]  = 18.0;
                p["B->K^*ccbar::Im{alpha_2^long}@GRvDV2021"]  = 19.0;

                Options o = { { "model", "WilsonScan" } };

                auto nc = NonlocalFormFactor<nc::PToV>::make("B->K^*::GRvDV2021", p, o);


                auto diagnostics = nc->diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair( 11.8899,  eps),         // Re{phi_long(q2 = 16.0)}
                    std::make_pair( -8.60714, eps),         // Im{phi_long(q2 = 16.0)}

                    std::make_pair( -6.07403, eps),         // Re{phi_perp(q2 = 16.0)}
                    std::make_pair( 9.3159,   eps)          // Im{phi_perp(q2 = 16.0)}
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp(-1.)), -2.9621,      eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp(-1.)), -4.09339,     eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_para(-1.)), -9.74979,  10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para(-1.)), -10.8811,  10*eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_long(-1.)), -0.412137,    eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long(-1.)), -0.44033,     eps);

                TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp(0.)), -2.96771,      eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp(0.)), -4.11808,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_para(0.)), -9.8699,    10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para(0.)), -11.0203,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_long(0.)),   0.,          eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long(0.)),   0.,          eps);

                TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp(4.)), -3.20673,      eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp(4.)), -4.55021,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_para(4.)), -11.2676,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para(4.)), -12.6111,   10*eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_long(4.)),   2.126,       eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long(4.)),   2.27378,     eps);

                TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp(12.)),   1.54707,    eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp(12.)),   2.4579,     eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_para(12.)),  7.01207,     eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para(12.)),  7.9229,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_long(12.)), -5.5288,      eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long(12.)), -5.93241,     eps);

                TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp_residue_jpsi()),  -0.381689,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp_residue_jpsi()),  -0.578581,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_perp_residue_psi2s()),  0.136974,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_perp_residue_psi2s()),  0.198421,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_para_residue_jpsi()),  -1.56304,    eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para_residue_jpsi()),  -1.75993,    eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_para_residue_psi2s()),  0.505653,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_para_residue_psi2s()),  0.5671,     eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_long_residue_jpsi()),   0.856417,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long_residue_jpsi()),   0.917859,   eps);
                TEST_CHECK_NEARLY_EQUAL(real(nc->H_long_residue_psi2s()), -0.519158,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(nc->H_long_residue_psi2s()), -0.555643,   eps);
            }
        }
} nonlocal_formfactor_grvdv2021_test;