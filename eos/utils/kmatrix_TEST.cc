#include <test/test.hh>
#include <eos/utils/kmatrix-impl.hh>
#include <eos/utils/complex.hh>
#include <iostream>
#include <memory>

#include <cmath>


using namespace test;
using namespace eos;


///////////////////////
/// Simplest K matrix, 1 channel, 1 resonance
///////////////////////
struct PPchan11 : public KMatrix<1,1>::Channel
{
    PPchan11(std::string name, double m1, double m2) : Channel(name, m1, m2)
    {
    };

    //Usefull definitions for beta and rho
    double mm = this->_m1 - this->_m2;
    double mp = this->_m1 + this->_m2;
    //sqrt of the Källen factor, defined with an absolute value
    double sqlk(const double & s) {
	return std::sqrt(std::abs((s-mp*mp)*(s-mm*mm)));
    }
    const double pi = M_PI;
    const complex<double> i = complex<double>(0.0, 1.0);

    double beta(const double & s) {
	if (s < mp*mp) {//Kinematic threshold
	    return 0.;
	}
	else {
	    return sqlk(s)/s;
	}
    }

    complex<double> rho(const double & s) {
	complex<double> result = 0.0;
	if (s < mm*mm){
	    result += 2*mm*(mp*mp-s)*std::log((mm+mp)/(mp-mm));
	    result += mp*sqlk(s)*std::log((mm*mm+mp*mp-2*s-2*sqlk(s))/(mm*mm+mp*mp-2*s+2*sqlk(s)));
	    result *= i/(2*mp*pi*s);
	    return result;
	}
	else if (s < mp*mp){
	    result += 2*mm*(mp*mp-s)*std::log((mm+mp)/(mp-mm));
	    //the pi is due to the simplification of the logs
	    result += 2*mp*sqlk(s)*(pi+std::atan(2*sqlk(s)/(mm*mm+mp*mp-2*s)));
	    result *= i/(2*mp*pi*s);
	    return result;
	}
	else {
	    result += 2*mm*(mp*mp-s)*std::log((mm+mp)/(mp-mm));
	    result += mp*sqlk(s)*std::log((mm*mm+mp*mp-2*s-2*sqlk(s))/(mm*mm+mp*mp-2*s+2*sqlk(s)));
	    result *= i/(2*mp*pi*s);
	    result += sqlk(s)/s;
	    return result;
	}
    }
};

struct res11 : public KMatrix<1,1>::Resonance
{
    res11(std::string name, double m) : Resonance(name, m)
    {
    };
};

double BreitWigner(double const s, double const M, double const Ga)
{
    return M*M*Ga*Ga / (pow(s - M*M, 2) + M*M*Ga*Ga);
}








///////////////////////
/// 1 channel, 2 resonances
///////////////////////
struct PPchan12 : public KMatrix<1,2>::Channel
{

    PPchan12(std::string name, double m1, double m2) : Channel(name, m1, m2)
    {
    };

    //Usefull definitions for beta and rho
    double mm = this->_m1 - this->_m2;
    double mp = this->_m1 + this->_m2;
    //sqrt of the Källen factor, defined with an absolute value
    double sqlk(const double & s) {
	return std::sqrt(std::abs((s-mp*mp)*(s-mm*mm)));
    }
    const double pi = M_PI;
    const complex<double> i = complex<double>(0.0, 1.0);

    double beta(const double & s) {
	if (s < mp*mp) {//Kinematic threshold
	    return 0.;
	}
	else {
	    return sqlk(s)/s;
	}
    }

    complex<double> rho(const double & s) {
	complex<double> result = 0.0;
	if (s < mm*mm){
	    result += 2*mm*(mp*mp-s)*std::log((mm+mp)/(mp-mm));
	    result += mp*sqlk(s)*std::log((mm*mm+mp*mp-2*s-2*sqlk(s))/(mm*mm+mp*mp-2*s+2*sqlk(s)));
	    result *= i/(2*mp*pi*s);
	    return result;
	}
	else if (s < mp*mp){
	    result += 2*mm*(mp*mp-s)*std::log((mm+mp)/(mp-mm));
	    //the pi is due to the simplification of the logs
	    result += 2*mp*sqlk(s)*(pi+std::atan(2*sqlk(s)/(mm*mm+mp*mp-2*s)));
	    result *= i/(2*mp*pi*s);
	    return result;
	}
	else {
	    result += 2*mm*(mp*mp-s)*std::log((mm+mp)/(mp-mm));
	    result += mp*sqlk(s)*std::log((mm*mm+mp*mp-2*s-2*sqlk(s))/(mm*mm+mp*mp-2*s+2*sqlk(s)));
	    result *= i/(2*mp*pi*s);
	    result += sqlk(s)/s;
	    return result;
	}
    }
};

struct res12 : public KMatrix<1,2>::Resonance
{
    res12(std::string name, double m) : Resonance(name, m)
    {
    };
};















///////////////////////
/// TESTS
///////////////////////

class KMatrixTest :
    public TestCase
{
public:
    KMatrixTest() :
	TestCase("KMatrix tests")
    {
    }

    virtual void run() const
    {

	constexpr double eps = 1e-4;



	// One channel, one resonnance:
	// In the limit where the resonance mass is much larger
	// than the channel masses, one recover a simple Breit-Wigner distribution
	{
	    auto Res = std::make_shared<res11>("Res", 15.0);
	    auto Chan = std::make_shared<PPchan11>("Chan", 0.7, 0.8);

	    Chan->g0s[0] = 1.0;

	    KMatrix<1,1> simplestKM({Chan}, {Res}, "simplestKM");

	    //Check |T-Matrix|^2 against Breit-Wigner
	    // The mass of the BW gets correction from the channels loop

	    double m = Res->_m;
	    double M = std::sqrt(m*m + Chan->rho(m*m).imag());

	    //s = 100. (everybody is ~zero there)
	    double Trowsq = std::norm(simplestKM.tmatrix_row(0, 100.)[0]);
	    TEST_CHECK_NEARLY_EQUAL(BreitWigner(100., M, 1./M),	Trowsq,	eps);

	    //s = 222. (resonance region)
	    Trowsq = std::norm(simplestKM.tmatrix_row(0, 222.)[0]);
	    TEST_CHECK_NEARLY_EQUAL(BreitWigner(222., M, 1./M),	Trowsq,	eps);

	    //s = 235. (resonance region)
	    Trowsq = std::norm(simplestKM.tmatrix_row(0, 235.)[0]);
	    TEST_CHECK_NEARLY_EQUAL(BreitWigner(235., M, 1./M),	Trowsq,	eps);
      
      
	    //s = 300. (everybody is ~zero there)
	    Trowsq = std::norm(simplestKM.tmatrix_row(0, 300.)[0]);
	    TEST_CHECK_NEARLY_EQUAL(BreitWigner(300., M, 1./M),	Trowsq,	eps);
      
      

	}





    
	// One channel, two resonances
	{
	    auto myres1 = std::make_shared<res12>("myres1", 1.0);
	    auto myres2 = std::make_shared<res12>("myres2", 2.0);


	    auto myPPchan = std::make_shared<PPchan12>("myPPchan", 0.7, 0.8);

	    myPPchan->g0s[0] = 1.1;
	    myPPchan->g0s[1] = 2.2;


	    //Test beta and rho for a PP channel
	    TEST_CHECK_NEARLY_EQUAL(myPPchan->beta(9.0),	0.865544,	eps);
	    TEST_CHECK_NEARLY_EQUAL(myPPchan->beta(1.5),	0.,	eps);

	    TEST_CHECK_NEARLY_EQUAL(myPPchan->rho(9.0).real(),	0.865544,	eps);
	    TEST_CHECK_NEARLY_EQUAL(myPPchan->rho(9.0).imag(),	0.724611,	eps);
	    TEST_CHECK_NEARLY_EQUAL(myPPchan->rho(1.5).real(),	0.,	eps);
	    TEST_CHECK_NEARLY_EQUAL(myPPchan->rho(1.5).imag(),	0.429317,	eps);
	    TEST_CHECK_NEARLY_EQUAL(myPPchan->rho(0.001).real(),	0.,	eps);
	    TEST_CHECK_NEARLY_EQUAL(myPPchan->rho(0.001).imag(),	0.635582,	eps);


	    //Test KMatrix inversion into TMatrix
	    KMatrix<1,2> myKM({myPPchan}, {myres1, myres2}, "myKM");

	    auto myK0ats0 = myKM.tmatrix_row(0, 9.0);
	    auto myK0ats1 = myKM.tmatrix_row(0, 1.5);
	    auto myK0ats2 = myKM.tmatrix_row(0, 0.001);

	    TEST_CHECK_NEARLY_EQUAL(myK0ats0[0].real(),	-0.217114,	eps);
	    TEST_CHECK_NEARLY_EQUAL(myK0ats0[0].imag(),	1.11299,	eps);
	    TEST_CHECK_NEARLY_EQUAL(myK0ats1[0].real(),	-0.610949,	eps);
	    TEST_CHECK_NEARLY_EQUAL(myK0ats1[0].imag(),	0.,	eps);
	    TEST_CHECK_NEARLY_EQUAL(myK0ats2[0].real(),	0.953701,	eps);
	    TEST_CHECK_NEARLY_EQUAL(myK0ats2[0].imag(),	0.,	eps);

	}





    
	TEST_CHECK_NEARLY_EQUAL(3.,	9.,	eps);


    }
} kmatrix_test;









// namespace ee_to_ccbar
// {
//   class PPChannel :
//     public KMatrix::Channel
//   {
//     const double m1;
//     const double m2;

//     virtual double beta(const double & s) const [[overwrite]]
//     {
//       if (s < pow(m1 + m2, 2))
//         return 0;
//       else
//         return std::sqrt((1.0 - (m1+m2)*(m1+m2)/s)*(1.0 - (m1-m2)*(m1-m2)/s))
//     }

      
//     virtual complex<double> rho(const double & s) const [[overwrite]]
//     {
//       if (s < pow(m1 - m2, 2))
//         //CHANGE THIS
//         complex<double> result = 1. + 2i
//           return result;
//       else if (s < pow(m1 + m2, 2)) 
//         //CHANGE THIS
//         complex<double> result = 1. + 2i
//           return result;
//       else 
//         //CHANGE THIS
//         complex<double> result = 1. + 2i
//           return result;
//     }
//   };
//
//
// // D^0 Dbar^0 channel
// template <>
// class DzDbarzChannel<EEToCCbarScattering::nresonances> :
//     public PPChannel<EEToCCbarScattering::nresonances>
// {
//     std::array<Parameter, nresonances> _g0;

//     std::array<Parameter, nresonances> _r;

//     std::string _par_name(i)

//     DzDbarzChannel(const Parameters & p) :
//         PPChannel(p["mass::D^0"], p["mass::D^0"]),
//         _g0
//         {{
//             p["ee->ccbar::g0(D^0Dbar^0,psi(2S)"],
//             p["ee->ccbar::g0(D^0Dbar^0,psi(2S)"],
//             p["ee->ccbar::g0(D^0Dbar^0,psi(2S)"],
//             p["ee->ccbar::g0(D^0Dbar^0,psi(2S)"],
//             p["ee->ccbar::g0(D^0Dbar^0,psi(2S)"],
//             ...
//         }}
//         _r
//         {{
//             p["ee->ccbar::r0(D^+D^-/D^0Dbar^0,psi(2S)"],
//             ...
//         }}
//     {
//     }

//     std:array<complex<double>, nresonances> g0() const
//     {
//         auto result = std::array<complex<double>, nresonances>;

//         for (unsigned i = 0 ; i < _g0.size() ; ++i)
//         {
//             result[i] = _g0[i]; // here

//             result[i] = _g0[i] * _r[i]; // for D+ D-
//         }
//     }
// };

//}
//class ee_to_ccbar
