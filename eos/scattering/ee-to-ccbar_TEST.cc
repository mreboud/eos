#include <test/test.hh>
#include <eos/scattering/ee-to-ccbar.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/kmatrix-impl.hh>

#include <vector>
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


        // Build K Matrix
        auto psi2S_res = std::make_shared<charmonium_resonance<3, 2>>("psi2S_res", p["mass::psi(2S)"]);
        auto psi3770_res = std::make_shared<charmonium_resonance<3, 2>>("psi3770_res", p["mass::psi(3770)"]);

        std::vector<Parameter> ee_g0s       {{p["ee->ccbar::g0(psi(2S),ee)"],           p["ee->ccbar::g0(psi(3770),ee)"]}};
        std::vector<Parameter> D0Dbar0_g0s  {{p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"],    p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"]}};
        std::vector<Parameter> DpDm_g0s     {{p["ee->ccbar::g0(psi(2S),D^+D^-)"],       p["ee->ccbar::g0(psi(3770),D^+D^-)"]}};

        auto ee_chan        = std::make_shared<SPPchan<3, 2>>("ee_chan", p["mass::e"], p["mass::e"], 3, ee_g0s);
        auto D0Dbar0_chan   = std::make_shared<SPPchan<3, 2>>("D0Dbar0_chan", p["mass::D^0"], p["mass::D^0"], 3, D0Dbar0_g0s);
        auto DpDm_chan      = std::make_shared<SPPchan<3, 2>>("DpDm_chan", p["mass::D^+"], p["mass::D^+"], 3, DpDm_g0s);

        KMatrix<3, 2> KMatrix32({ee_chan, D0Dbar0_chan, DpDm_chan}, {psi2S_res, psi3770_res}, "KMatrix");



        auto KMatrix32s0 = KMatrix32.tmatrix_row(0, 9.0);
        auto KMatrix32s1 = KMatrix32.tmatrix_row(0, 1.5);
        auto KMatrix32s2 = KMatrix32.tmatrix_row(0, 0.001);


        TEST_CHECK_NEARLY_EQUAL(KMatrix32s0[0].real(), 0.00870921,  eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s0[0].imag(), 0.000075856, eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s1[0].real(), 0.0037429,   eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s1[0].imag(), 0.0000140095,eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s2[0].real(), 0.00338731,  eps);
        TEST_CHECK_NEARLY_EQUAL(KMatrix32s2[0].imag(), 0.000011468, eps);



    }
} eetoccbar_test;