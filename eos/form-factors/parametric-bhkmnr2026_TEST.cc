/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Fatemeh Nouri
 * Copyright (c) 2026 Méril Reboud
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
#include <eos/form-factors/parametric-bhkmnr2026.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class ParametricBHKMNR2026Test :
    public TestCase
{
    public:
        ParametricBHKMNR2026Test() :
            TestCase("parametric_BHKMNR2026_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-7;

            {

                Parameters p = Parameters::Defaults();
                p["mass::pi^+"]                        =  0.13957;
                p["0->pipi::s_0@BHKMNR2026"]           =  0.0;
                p["0->pipi::s_in@BHKMNR2026"]          =  1.0;
                p["0->pipi::a_(+,1)^4@BHKMNR2026"]     =  0.20;
                p["0->pipi::a_(+,1)^5@BHKMNR2026"]     =  0.12;
                p["0->pipi::a_(+,1)^6@BHKMNR2026"]     =  0.07;
                p["0->pipi::a_(+,1)^7@BHKMNR2026"]     =  0.02;
                p["0->pipi::a_(+,1)^8@BHKMNR2026"]     =  0;
                p["0->pipi::a_(+,1)^9@BHKMNR2026"]     =  0;
                p["0->pipi::a_(+,1)^10@BHKMNR2026"]    =  0;
                p["0->pipi::a_(+,1)^11@BHKMNR2026"]    =  0;
                p["0->pipi::a_(+,1)^12@BHKMNR2026"]    =  0;
                p["0->pipi::M_(+,1,0)@BHKMNR2026"]     =  0.760895;
                p["0->pipi::Gamma_(+,1,0)@BHKMNR2026"] =  0.146155;

                /* 0->PP factory */
                {
                    std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->pipi::BHKMNR2026", p, Options{ });

                    TEST_CHECK(nullptr != ff);
                }


                {
                    Options o{ { "n-resonances"_ok, "1" } };
                    BHKMNR2026FormFactors<VacuumToPiPi> ff(p, o);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(0.0)),                                                0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(0.0)),                                                0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(0.077919139600000)),                                 -0.24972020,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(0.077919139600000)),                                  0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(0.5)),                                               -0.02968043,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(0.5)),                                                0.55209340,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(1.0)),                                                0.12710102,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(1.0)),                                                0.99188978,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(complex<double>(0.5, 0.5))),                          0.22651223,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(complex<double>(0.5, 0.5))),                         -0.44265917,    eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(0.0)),                                                  2.53121472,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(0.0)),                                                  0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(0.5)),                                                  0.35146403,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(0.5)),                                                  0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(complex<double>(0.5, 0.5))),                           -0.53074826,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(complex<double>(0.5, 0.5))),                           -0.43239319,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(1.0)),                                                  0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(1.0)),                                                  0.00000000,    eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.dPdpsi(0.0)),                                            -5.91127510,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dPdpsi(0.0)),                                             0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dPdpsi(0.5)),                                            -1.96542644,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dPdpsi(0.5)),                                             0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dPdpsi(complex<double>(0.5, 0.5))),                       1.42575563,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dPdpsi(complex<double>(0.5, 0.5))),                       2.62969011,    eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms(0,0.0)),                                    -5.91127510,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms(0,0.0)),                                     0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms(1,0.0)),                                     2.53121472,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms(1,0.0)),                                     0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms(0,0.5)),                                    -1.96542644,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms(0,0.5)),                                     0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms(1,0.5)),                                    -0.63124919,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms(1,0.5)),                                     0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms(1,complex<double>(0.5, 0.5))),              -1.13271549,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms(1,complex<double>(0.5, 0.5))),               1.59532969,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms(2,complex<double>(0.5, 0.5))),              -1.41320012,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms(2,complex<double>(0.5, 0.5))),              -0.25026362,    eps);


                    const auto constrained_a = ff.constrained_a_fp_I1();
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[0],                                                 0.39506723,    eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[1],                                                 0.57968686,    eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[2],                                                 0.87223762,    eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[3],                                                 0.45874982,    eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(0.0)),                                                1.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(0.0)),                                                0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(0.5)),                                                3.44706708,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(0.5)),                                                1.48145360,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(complex<double>(0.5, 0.5))),                          0.67900787,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(complex<double>(0.5, 0.5))),                          0.88087501,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(1.0)),                                               -0.49337270,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(1.0)),                                                1.32130884,    eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(0.0)),                                        0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(0.0)),                                        0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(0.077919139600000)),                          0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(0.077919139600000)),                          0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(complex<double>(0.5, 0.5))),                 -0.87389460,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(complex<double>(0.5, 0.5))),                  0.22494469,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(1.0)),                                       -0.34127181,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(1.0)),                                        0.91396516,    eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.scattering_lenght_parameters()[0]),                       -1.84127173,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.scattering_lenght_parameters()[0]),                        0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.scattering_lenght_parameters()[1]),                      -14.56292832,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.scattering_lenght_parameters()[1]),                        0.00000000,    eps);

                    //needed for charged pion radius
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(0.0)),                                          -0.86803920,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_11(0.0)),                                           0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(0.077919139600000)),                             0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_11(0.077919139600000)),                             0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(complex<double>(0.5,0.5))),                     -3.00934401,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_11(complex<double>(0.5,0.5))),                     -1.07084185,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(1.0)),                                           0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(1.0)),                                           0.00000000,    eps);

                    // Lowering the precision to align with the integration precision
                    TEST_CHECK_RELATIVE_ERROR(ff.saturation(),                                                 0.26055906,    1e-5);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p_of_psi(complex<double>(0.5, 0.5))),                    0.08387160,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p_of_psi(complex<double>(0.5, 0.5))),                   -0.65190022,    eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p_of_psi(complex<double>(0.7, 0.3))),                    0.03490867,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p_of_psi(complex<double>(0.7, 0.3))),                   -0.25944037,    eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.abs2_f_p_of_psi( 0.5, 0.5),                                     0.43200834,    eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.abs2_f_p_of_psi( 0.7, 0.3),                                     0.06852792,    eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.abs2_f_p_of_psi( 0.0, 0.0),                                     1.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.abs2_f_p_of_psi(-1.0, 0.0),                                     1.29262147,    eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.arg_f_p_of_psi( 0.5, 0.5),                                     -1.44284212,    eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.arg_f_p_of_psi( 0.7, 0.3),                                     -1.43704590,    eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.arg_f_p_of_psi( 0.0, 0.0),                                      0.00000000,    eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.arg_f_p_of_psi(-1.0, 0.0),                                      0.00000000,    eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.re_residue_rho(),                                               0.04063694,    eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.im_residue_rho(),                                               0.25969762,    eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.root_penalty(),                                                 1.00000000,    eps);
                }

                p["mass::pi^+"]                        =  0.13957;
                p["0->pipi::s_0@BHKMNR2026"]           =  0.0;
                p["0->pipi::s_in@BHKMNR2026"]          =  0.85051;
                p["0->pipi::a_(+,1)^4@BHKMNR2026"]     = -0.54508094;
                p["0->pipi::a_(+,1)^5@BHKMNR2026"]     =  0.06982185;
                p["0->pipi::a_(+,1)^6@BHKMNR2026"]     = -0.3837015;
                p["0->pipi::a_(+,1)^7@BHKMNR2026"]     =  0.23947605;
                p["0->pipi::a_(+,1)^8@BHKMNR2026"]     = -0.23656627;
                p["0->pipi::a_(+,1)^9@BHKMNR2026"]     =  0.09873716;
                p["0->pipi::a_(+,1)^10@BHKMNR2026"]    =  0.;
                p["0->pipi::a_(+,1)^11@BHKMNR2026"]    =  0.001;
                p["0->pipi::a_(+,1)^12@BHKMNR2026"]    =  0;
                p["0->pipi::M_(+,1,0)@BHKMNR2026"]     =  0.76000036;
                p["0->pipi::Gamma_(+,1,0)@BHKMNR2026"] =  0.14483015;

                {
                    Options o{ { "n-resonances"_ok, "1" } };
                    BHKMNR2026FormFactors<VacuumToPiPi> ff(p, o);

                    TEST_CHECK_NEARLY_EQUAL(ff.root_penalty(),                                                 1.00000000,    eps);
                }
            }
        }
} parametric_BHKMNR2026_test;
