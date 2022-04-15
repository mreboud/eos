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
        p["ee->ccbar::g0(psi(2S),ee)"] = 0.1;
        p["ee->ccbar::g0(psi(3770),ee)"] = 0.2;
        p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"] = 0.3;
        p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"] = 0.4;
        p["ee->ccbar::g0(psi(2S),D^+D^-)"] = 0.5;
        p["ee->ccbar::g0(psi(3770),D^+D^-)"] = 0.6;


        //Set all cst to zero
        p["ee->ccbar::c_0(ee,ee)"] = 0.0;
        p["ee->ccbar::c_0(ee,D^0Dbar^0)"] = 0.0;
        p["ee->ccbar::c_0(ee,D^+D^-)"] = 0.0;
        p["ee->ccbar::c_0(D^0Dbar^0,D^0Dbar^0)"] = 0.0;
        p["ee->ccbar::c_0(D^0Dbar^0,D^+D^-)"] = 0.0;
        p["ee->ccbar::c_0(D^+D^-,D^+D^-)"] = 0.0;

        // Build K Matrix
        auto psi2S_res = std::make_shared<charmonium_resonance<3, 2, 0>>("psi2S_res", p["mass::psi(2S)"], p["size::psi(2S)"]);
        auto psi3770_res = std::make_shared<charmonium_resonance<3, 2, 0>>("psi3770_res", p["mass::psi(3770)"], p["size::psi(3770)"]);

        std::array<std::array<std::array<Parameter, 3>, 3>, 1> bkgcst {
            p["ee->ccbar::c_0(ee,ee)"],        p["ee->ccbar::c_0(ee,D^0Dbar^0)"],        p["ee->ccbar::c_0(ee,D^+D^-)"],
            p["ee->ccbar::c_0(ee,D^0Dbar^0)"], p["ee->ccbar::c_0(D^0Dbar^0,D^0Dbar^0)"], p["ee->ccbar::c_0(D^0Dbar^0,D^+D^-)"],
            p["ee->ccbar::c_0(ee,D^+D^-)"],    p["ee->ccbar::c_0(D^0Dbar^0,D^+D^-)"],    p["ee->ccbar::c_0(D^+D^-,D^+D^-)"],
        };

        std::array<Parameter, 2> ee_g0s       {{p["ee->ccbar::g0(psi(2S),ee)"],           p["ee->ccbar::g0(psi(3770),ee)"]}};
        std::array<Parameter, 2> D0Dbar0_g0s  {{p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"],    p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"]}};
        std::array<Parameter, 2> DpDm_g0s     {{p["ee->ccbar::g0(psi(2S),D^+D^-)"],       p["ee->ccbar::g0(psi(3770),D^+D^-)"]}};

        auto ee_chan        = std::make_shared<SPPchan<3, 2, 0>>("ee_chan", p["mass::e"], p["mass::e"], 3, ee_g0s);
        auto D0Dbar0_chan   = std::make_shared<SPPchan<3, 2, 0>>("D0Dbar0_chan", p["mass::D^0"], p["mass::D^0"], 3, D0Dbar0_g0s);
        auto DpDm_chan      = std::make_shared<SPPchan<3, 2, 0>>("DpDm_chan", p["mass::D^+"], p["mass::D^+"], 3, DpDm_g0s);

        KMatrix<3, 2, 0> KMatrix32({ee_chan, D0Dbar0_chan, DpDm_chan}, {psi2S_res, psi3770_res}, bkgcst, "KMatrix");


        auto KMatrix32s0 = KMatrix32.tmatrix_row(0, 20.0);
        auto KMatrix32s1 = KMatrix32.tmatrix_row(0, 13.94);
        auto KMatrix32s2 = KMatrix32.tmatrix_row(0,  0.0001);

        TEST_CHECK_NEARLY_EQUAL(KMatrix32s0[0].real(), -0.0084434, eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s0[0].imag(),  0.0007069, eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s1[0].real(),  0.1031590, eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s1[0].imag(),  0.0122801, eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s2[0].real(),  0.0035448, eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s2[0].imag(),  0.0000125, eps);


        //Set ee -> ee cst to .5
        p["ee->ccbar::c_0(ee,ee)"] = 0.5;

        KMatrix32s0 = KMatrix32.tmatrix_row(0, 20.0);
        KMatrix32s1 = KMatrix32.tmatrix_row(0, 13.94);
        KMatrix32s2 = KMatrix32.tmatrix_row(0,  0.0001);

        TEST_CHECK_NEARLY_EQUAL(KMatrix32s0[0].real(),  0.395486, eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s0[0].imag(),  0.194911, eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s1[0].real(),  0.441778, eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s1[0].imag(),  0.268202, eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s2[0].real(),  0.402543, eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s2[0].imag(),  0.201637, eps);

        // Test the full K matrix

        Options oo;

        EEToCCBar c(p, oo);

        auto ir = c.prepare(4.5);
        TEST_CHECK_RELATIVE_ERROR(c.sigma_eetoD0Dbar0(ir), 129.76254216, eps);

    }
} eetoccbar_test;