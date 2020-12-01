#include <test/test.hh>
#include <eos/utils/kmatrix-impl.hh>
#include <eos/utils/complex.hh>
#include <iostream>
#include <memory>

#include <cmath>


using namespace test;
using namespace eos;


struct PPchan : public KMatrix<1,2>::Channel
{

  PPchan(std::string name, double m1, double m2) : Channel(name, m1, m2)
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

struct myres : public KMatrix<1,2>::Resonance
{

  myres(std::string name, double m) : Resonance(name, m)
  {
  };

};

class MyTest :
  public TestCase
{
public:
  MyTest() :
    TestCase("my_test")
  {
  }

  virtual void run() const
  {

    constexpr double eps = 1e-4;


    auto myres1 = std::make_shared<myres>("myres1", 1.0);
    auto myres2 = std::make_shared<myres>("myres2", 2.0);


    auto myPPchan = std::make_shared<PPchan>("myPPchan", 0.7, 0.8);

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


    KMatrix<1,2> myKM({myPPchan}, {myres1, myres2}, "myKM");

    auto myK0ats0 = myKM.tmatrix_row(0, 9.0);
    auto myK0ats1 = myKM.tmatrix_row(0, 1.5);
    auto myK0ats2 = myKM.tmatrix_row(0, 0.001);

    TEST_CHECK_NEARLY_EQUAL(myK0ats0[0].real(),	0.635582,	eps);
    TEST_CHECK_NEARLY_EQUAL(myK0ats0[0].imag(),	0.635582,	eps);


    TEST_CHECK_NEARLY_EQUAL(3.,	9.,	eps);


  }
} my_test;









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
