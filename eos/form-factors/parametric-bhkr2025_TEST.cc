/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Méril Reboud
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
#include <eos/form-factors/parametric-bhkr2025.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class ParametricBHKR2025Test :
    public TestCase
{
    public:
        ParametricBHKR2025Test() :
            TestCase("parametric_BHKR2025_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {

                Parameters p = Parameters::Defaults();
                p["mass::pi^+"]                    =  0.13957;
                p["0->pipi::s_th@BHKR2025"]        =  0.976144;
                p["0->pipi::b_(+,1)^4@BHKR2025"]   =  0.1;
                p["0->pipi::b_(+,1)^5@BHKR2025"]   = -0.2;
                p["0->pipi::M_(+,0)@BHKR2025"]     =  0.77;
                p["0->pipi::Gamma_(+,0)@BHKR2025"] =  0.15;

                Options o
                {
                    {"n-resonances-1m"_ok, "1"}
                };

                /* 0->PP factory */
                {
                    std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->pipi::BHKR2025", p, o);

                    TEST_CHECK(nullptr != ff);
                }

                /* f_+ at timelike q2 > 0.0 */
                {
                    BHKR2025FormFactors<VacuumToPiPi> ff(p, o);

                    TEST_CHECK_RELATIVE_ERROR(real(ff.zeta( 0.0)),      0.144203,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.zeta( 0.0)),        0.0,        eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.zeta(+0.1)),        0.0,        eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.zeta(+0.1)),      0.0788823,  eps);
                    TEST_CHECK_RELATIVE_ERROR(real(ff.zeta(+1.1)),      0.34811 ,   eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.zeta(+1.1)),      0.937454,   eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.dzetadq2(+0.1)),    0.0,        eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.dzetadq2(+0.1)),  1.80858,    eps);
                    TEST_CHECK_RELATIVE_ERROR(real(ff.dzetadq2(+1.1)),  1.235,      eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.dzetadq2(+1.1)),  0.458601,   eps);

                    TEST_CHECK_RELATIVE_ERROR(real(ff.xi(complex<double>(0.5, 0.2))),  0.641326,  eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.xi(complex<double>(0.5, 0.2))),  0.169533,  eps);
                    TEST_CHECK_RELATIVE_ERROR(real(ff.xi(complex<double>(0.5, 0.5))),  0.696735,  eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.xi(complex<double>(0.5, 0.5))),  0.397384,  eps);
                    TEST_CHECK_RELATIVE_ERROR(real(ff.xi(-ff.zeta( 0.0))),            -1.0,       eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.xi(-ff.zeta( 0.0))),               0.0,       eps);

                    TEST_CHECK_RELATIVE_ERROR(real(ff.resonance_product_p(complex<double>(0.5, 0.2))),        1.65683,  eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.resonance_product_p(complex<double>(0.5, 0.2))),       -0.770885, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(ff.resonance_product_p_prime(complex<double>(0.5, 0.2))), -3.50436,  eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.resonance_product_p_prime(complex<double>(0.5, 0.2))),  2.08812,  eps);

                    TEST_CHECK_RELATIVE_ERROR(ff.b_0(),  0.274876,  eps);
                    TEST_CHECK_RELATIVE_ERROR(ff.b_1(),  0.0893618, eps);
                    TEST_CHECK_RELATIVE_ERROR(ff.b_2(), -0.518187,  eps);
                    TEST_CHECK_RELATIVE_ERROR(ff.b_3(), -0.0921298, eps);
                    TEST_CHECK_RELATIVE_ERROR(ff.b_4(), -0.0376633, eps);

                    TEST_CHECK_RELATIVE_ERROR(real(ff.f_p( 0.0)),  1.0,       eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p( 0.0)),    0.0,       eps);
                    TEST_CHECK_RELATIVE_ERROR(real(ff.f_p(+0.1)),  1.45629,   eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.f_p(+0.1)),  0.0380841, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(ff.f_p(+0.5)),  3.21735,   eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.f_p(+0.5)),  6.26233,   eps);

                    TEST_CHECK_RELATIVE_ERROR(ff.re_f_p_zeta(+0.1, -0.1), 1.20826,  eps);

                    // TEST_CHECK_NEARLY_EQUAL(ff.dFdq2_q2eq0(),                      1.93155736,    eps);
                    // TEST_CHECK_NEARLY_EQUAL(ff.r_pi_squared(),                     0.451265116,   eps);

                    // TEST_CHECK_NEARLY_EQUAL(ff.re_residue_rho(),                  -0.686247061,   eps);
                    // TEST_CHECK_NEARLY_EQUAL(ff.im_residue_rho(),                   0.272035152,   eps);

                    // TEST_CHECK_NEARLY_EQUAL(ff.re_residue_rho_q2(),               -0.708332619,   eps);
                    // TEST_CHECK_NEARLY_EQUAL(ff.im_residue_rho_q2(),                0.135266983,   eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.dispersive_integrand(complex<double>(0.0,0.5)),              5.27289,   eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.dispersive_integrand(complex<double>(sqrt(0.5),sqrt(0.5))),  0.0622617, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.saturation(),                                                1.65759,   eps);
                }
            }

            {

                Parameters p = Parameters::Defaults();
                p["mass::pi^+"]                    =  0.13957;
                p["0->pipi::s_th@BHKR2025"]        =  0.976144;
                p["0->pipi::b_(+,1)^4@BHKR2025"]   =  0.1;
                p["0->pipi::b_(+,1)^5@BHKR2025"]   = -0.2;
                p["0->pipi::M_(+,0)@BHKR2025"]     =  0.77;
                p["0->pipi::Gamma_(+,0)@BHKR2025"] =  0.15;
                p["0->pipi::M_(+,1)@BHKR2025"]     =  1.3;
                p["0->pipi::Gamma_(+,1)@BHKR2025"] =  0.5;

                Options o
                {
                    {"n-resonances-1m"_ok, "2"}
                };

                /* 0->PP factory */
                {
                    std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->pipi::BHKR2025", p, o);

                    TEST_CHECK(nullptr != ff);
                }

                /* f_+ at timelike q2 > 0.0 */
                {
                    BHKR2025FormFactors<VacuumToPiPi> ff(p, o);

                    TEST_CHECK_RELATIVE_ERROR(real(ff.resonance_product_p(complex<double>(0.5, 0.2))),        1.20704,  eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.resonance_product_p(complex<double>(0.5, 0.2))),       -0.736724, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(ff.resonance_product_p_prime(complex<double>(0.5, 0.2))), -3.22779,  eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.resonance_product_p_prime(complex<double>(0.5, 0.2))),  2.08132,  eps);
                }
            }

            {

                Parameters p = Parameters::Defaults();
                p["mass::pi^+"]                    =  0.13957;
                p["0->pipi::s_th@BHKR2025"]        =  0.976144;
                p["0->pipi::b_(+,1)^4@BHKR2025"]   =  0.1;
                p["0->pipi::b_(+,1)^5@BHKR2025"]   = -0.2;
                p["0->pipi::M_(+,0)@BHKR2025"]     =  0.77;
                p["0->pipi::Gamma_(+,0)@BHKR2025"] =  0.15;
                p["0->pipi::M_(+,1)@BHKR2025"]     =  1.3;
                p["0->pipi::Gamma_(+,1)@BHKR2025"] =  0.5;
                p["0->pipi::M_(+,2)@BHKR2025"]     =  1.7;
                p["0->pipi::Gamma_(+,2)@BHKR2025"] =  0.3;

                Options o
                {
                    {"n-resonances-1m"_ok, "3"}
                };

                /* 0->PP factory */
                {
                    std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->pipi::BHKR2025", p, o);

                    TEST_CHECK(nullptr != ff);
                }

                /* f_+ at timelike q2 > 0.0 */
                {
                    BHKR2025FormFactors<VacuumToPiPi> ff(p, o);

                    TEST_CHECK_RELATIVE_ERROR(real(ff.resonance_product_p(complex<double>(0.5, 0.2))),        1.32464,  eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.resonance_product_p(complex<double>(0.5, 0.2))),       -0.713213, eps);
                    TEST_CHECK_RELATIVE_ERROR(real(ff.resonance_product_p_prime(complex<double>(0.5, 0.2))), -3.07668,  eps);
                    TEST_CHECK_RELATIVE_ERROR(imag(ff.resonance_product_p_prime(complex<double>(0.5, 0.2))),  1.53294,  eps);
                }
            }
        }
} parametric_BHKR2025_test;
