/*
 * Copyright (c) 2021 Méril Reboud
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
#include <eos/maths/complex.hh>
#include <eos/scattering/ee-to-ccbar.hh>
#include <eos/utils/kmatrix-impl.hh>

#include <array>
#include <memory>


using namespace test;
using namespace eos;

class eetoccbarTest :
    public TestCase
{
    public:
    eetoccbarTest() :
        TestCase("ee->ccbar scattering tests")
        {
        }

        virtual void run() const
        {

            constexpr double eps = 1e-5;

            Parameters p = Parameters::Defaults();
            p["ee->ccbar::g0(psi(2S),e^+e^-)"] = 0.1;
            p["ee->ccbar::g0(psi(3770),e^+e^-)"] = 0.2;
            p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"] = 0.3;
            p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"] = 0.4;
            p["ee->ccbar::g0(psi(2S),D^+D^-)"] = 0.5;
            p["ee->ccbar::g0(psi(3770),D^+D^-)"] = 0.6;


            // Set all cst to zero
            p["ee->ccbar::c(e^+e^-,e^+e^-)"] = 0.0;
            p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"] = 0.0;
            p["ee->ccbar::c(e^+e^-,D^+D^-)"] = 0.0;
            p["ee->ccbar::c(D^0Dbar^0,D^0Dbar^0)"] = 0.0;
            p["ee->ccbar::c(D^0Dbar^0,D^+D^-)"] = 0.0;
            p["ee->ccbar::c(D^+D^-,D^+D^-)"] = 0.0;

            // Build K Matrix
            auto psi2S = std::make_shared<CharmoniumResonance<3, 2>>("psi2S", p["mass::psi(2S)"], p["psi(2S)::q_R"]);
            auto psi3770 = std::make_shared<CharmoniumResonance<3, 2>>("psi3770", p["mass::psi(3770)", p["psi(3770)::q_R"]]);

            std::array<std::array<Parameter, 3>, 3> bkgcst {
                p["ee->ccbar::c(e^+e^-,e^+e^-)"],    p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"],    p["ee->ccbar::c(e^+e^-,D^+D^-)"],
                p["ee->ccbar::c(e^+e^-,D^0Dbar^0)"], p["ee->ccbar::c(D^0Dbar^0,D^0Dbar^0)"], p["ee->ccbar::c(D^0Dbar^0,D^+D^-)"],
                p["ee->ccbar::c(e^+e^-,D^+D^-)"],    p["ee->ccbar::c(D^0Dbar^0,D^+D^-)"],    p["ee->ccbar::c(D^+D^-,D^+D^-)"],
            };

            std::array<Parameter, 2> ee_g0s       {{p["ee->ccbar::g0(psi(2S),e^+e^-)"],       p["ee->ccbar::g0(psi(3770),e^+e^-)"]}};
            std::array<Parameter, 2> D0Dbar0_g0s  {{p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"],    p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"]}};
            std::array<Parameter, 2> DpDm_g0s     {{p["ee->ccbar::g0(psi(2S),D^+D^-)"],       p["ee->ccbar::g0(psi(3770),D^+D^-)"]}};

            auto ee        = std::make_shared<EffChannel<3, 2>>("ee", p["mass::e"], p["mass::e"], ee_g0s);
            auto D0Dbar0   = std::make_shared<PWavePPChannel<3, 2>>("D0Dbar0", p["mass::D^0"], p["mass::D^0"], D0Dbar0_g0s);
            auto DpDm      = std::make_shared<PWavePPChannel<3, 2>>("DpDm", p["mass::D^+"], p["mass::D^+"], DpDm_g0s);

            KMatrix<3, 2> KMatrix32({ee, D0Dbar0, DpDm}, {psi2S, psi3770}, bkgcst, "ee->ccbar");


            auto KMatrix32s0 = KMatrix32.tmatrix_row(0, 20.0);
            auto KMatrix32s1 = KMatrix32.tmatrix_row(0, 13.94);
            auto KMatrix32s2 = KMatrix32.tmatrix_row(0,  0.0001);

            TEST_CHECK_NEARLY_EQUAL(KMatrix32s0[0].real(), -0.0084963, eps);
            TEST_CHECK_NEARLY_EQUAL(KMatrix32s0[0].imag(),  0.0002653, eps);
            TEST_CHECK_NEARLY_EQUAL(KMatrix32s1[0].real(),  0.1034849, eps);
            TEST_CHECK_NEARLY_EQUAL(KMatrix32s1[0].imag(),  0.0108295, eps);
            TEST_CHECK_NEARLY_EQUAL(KMatrix32s2[0].real(),  0.0035448, eps);
            TEST_CHECK_NEARLY_EQUAL(KMatrix32s2[0].imag(),  0.0000125, eps);


            // Set ee -> ee cst to .5
            p["ee->ccbar::c(e^+e^-,e^+e^-)"] = 0.5;

            KMatrix32s0 = KMatrix32.tmatrix_row(0, 20.0);
            KMatrix32s1 = KMatrix32.tmatrix_row(0, 13.94);
            KMatrix32s2 = KMatrix32.tmatrix_row(0,  0.0001);

            TEST_CHECK_NEARLY_EQUAL(KMatrix32s0[0].real(),  0.395745, eps);
            TEST_CHECK_NEARLY_EQUAL(KMatrix32s0[0].imag(),  0.194664, eps);
            TEST_CHECK_NEARLY_EQUAL(KMatrix32s1[0].real(),  0.442759, eps);
            TEST_CHECK_NEARLY_EQUAL(KMatrix32s1[0].imag(),  0.267702, eps);
            TEST_CHECK_NEARLY_EQUAL(KMatrix32s2[0].real(),  0.404228, eps);
            TEST_CHECK_NEARLY_EQUAL(KMatrix32s2[0].imag(),  0.200366, eps);

            // Test the full K matrix

            Options oo;

            EEToCCBar c(p, oo);

            auto ir = c.prepare(4.5);
            TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoD0Dbar0(ir), 129.76254216, eps);

        }
} eetoccbar_test;