#ifndef EOS_GUARD_EOS_SCATTERING_EETOCCBAR_HH
#define EOS_GUARD_EOS_SCATTERING_EETOCCBAR_HH 1


namespace ee_to_ccbar
{
    template <unsigned nchannels_, unsigned nresonances_>
    struct PPchan : public KMatrix<nchannels_, nresonances_>::Channel
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
		return pow(sqlk(s),3)/s/s;
	    }
        }
	
	complex<double> rho(const double & s) {
	    complex<double> result = 0.0;
	    if (s < mm*mm){
		result += mm*sqlk(s)*(2*mm*mp*s-(-3*mp*mp*s+mm*mm*(2*mp*mp+s))*std::log(_m1/_m2));
		result += 2*pow(mp,3)*pow(mm*mm-s,2)*std::log((-2*(s+sqlk(s))+mm*mm+mp*mp)/(mp*mp-mm*mm));
		result *= i*pow(mp*mp-s,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(mm*mm-s));
		return result;
	    }
	    else if (s < mp*mp){
		result += mm*sqlk(s)*(2*mm*mp*s-(-3*mp*mp*s+mm*mm*(2*mp*mp+s))*std::log(_m1/_m2));
		result += 4*pow(mp,3)*pow(s-mm*mm,2)*std::arctan(sqrt((s-mm*mm)/(mp*mp-s)));
		result *= i*pow(mp*mp-s,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(s-mm*mm));
		return result;
	    }
	    else {
		result += mm*sqlk(s)*(-2*mm*mp*s+(-3*mp*mp*s+mm*mm*(2*mp*mp+s))*std::log(_m1/_m2));
		result += 2*pow(mp,3)*pow(s-mm*mm,2)*std::log((2*(s+sqlk(s))-mm*mm-mp*mp)/(mp*mp-mm*mm));
		result *= i*pow(s-mp*mp,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(s-mm*mm));
		result += pow(sqlk(s),3)/s/s;
		return result;
	    }
	}
    };




    // Hardcode 3 channels, 2 resonances
    struct eeChannel32 : public PPchan<3, 2>::Channel
    {
	std::array<double, nresonances_> g0s;

	eeChannel32(const Parameters & p) :
	    PPchan(p["mass::e"], p["mass::e"]),
	    _g0
	    {{
		p["ee->ccbar::g0(psi(2S),ee)"],
		p["ee->ccbar::g0(psi(3770),ee)"]
	    }}
	{
	}
    };

    struct D0Dbar0Channel32 : public PPchan<3, 2>::Channel
    {
	std::array<double, nresonances_> g0s;

	D0Dbar0Channel32(const Parameters & p) :
	    PPchan(p["mass::D^0"], p["mass::D^0"]),
	    _g0
	    {{
		p["ee->ccbar::g0(psi(2S),D^0Dbar^0)"],
		p["ee->ccbar::g0(psi(3770),D^0Dbar^0)"]
	    }}
	{
	}
    };


    struct DpDmChannel32 : public PPchan<3, 2>::Channel
    {
	std::array<double, nresonances_> g0s;

	DpDmChannel32(const Parameters & p) :
	    PPchan(p["mass::D^+"], p["mass::D^+"]),
	    _g0
	    {{
		p["ee->ccbar::g0(psi(2S),D^pD^m)"],
		p["ee->ccbar::g0(psi(3770),D^pD^m)"]
	    }}
	{
	}	
    };
}
