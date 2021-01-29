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

#ifndef EOS_GUARD_EOS_SCATTERING_EETOCCBAR_HH
#define EOS_GUARD_EOS_SCATTERING_EETOCCBAR_HH 1

#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/kmatrix-impl.hh>
#include <vector>

namespace eos

{
    // S -> PP channel
    template <unsigned nchannels_, unsigned nresonances_>
    struct SPPchan :
    public KMatrix<nchannels_, nresonances_>::Channel
    {

        SPPchan(std::string name, double m1, double m2, unsigned N_orbital, std::vector<Parameter> g0s) :
        KMatrix<nchannels_, nresonances_>::Channel(name, m1, m2, N_orbital, g0s)
        {
        };

        //Usefull definitions for beta and rho
        const double mm = this->_m1 - this->_m2;
        const double mp = this->_m1 + this->_m2;
        //sqrt of the Källen factor, defined with an absolute value
        double sqlk(const double & s) {
            return std::sqrt(std::abs((s-mp*mp)*(s-mm*mm)));
        }
        const double pi = M_PI;
        const complex<double> i =  complex<double>(0.0, 1.0);

        const unsigned N_orbital = 3;

        double beta(const double & s){
            if (s < mp*mp) { return 0.; } //Kinematic threshold
            else { return sqlk(s)/s; }
        }

        complex<double> rho(const double & s) {
            complex<double> result = 0.0;
            if (s < mm*mm){
                // result += mm*(mp*mp-s)*std::log((mm+mp)/(mp-mm));
                // result += mp*sqlk(s)*std::log((mm*mm+mp*mp-2*s-2*sqlk(s))/(mp*mp-mm*mm));
                // result *= i/(mp*pi*s);
                return result;
            }
            else if (s < mp*mp){
                // result += mm*(mp*mp-s)*std::log((mm+mp)/(mp-mm));
                // result += 2*mp*sqlk(s)*std::atan(std::sqrt((s-mm*mm)/(mp*mp-s)));
                // result *= i/(mp*pi*s);
                return result;
            }
            else {
                // result += mm*(mp*mp-s)*std::log((mm+mp)/(mp-mm));
                // result += mp*sqlk(s)*std::log((2*s+2*sqlk(s)-mm*mm-mp*mp)/(mp*mp-mm*mm));
                // result *= i/(mp*pi*s);
                result += sqlk(s)/s;
                return result;
            }
        }
    };


    // V -> PP channel
    template <unsigned nchannels_, unsigned nresonances_>
    struct PPchan :
    public KMatrix<nchannels_, nresonances_>::Channel
    {

        PPchan(std::string name, double m1, double m2, unsigned N_orbital, std::vector<Parameter> g0s) :
        KMatrix<nchannels_, nresonances_>::Channel(name, m1, m2, N_orbital, g0s)
        {
        };

        //Usefull definitions for beta and rho
        const double mm = this->_m1 - this->_m2;
        const double mp = this->_m1 + this->_m2;
        //sqrt of the Källen factor, defined with an absolute value
        double sqlk(const double & s) {
            return std::sqrt(std::abs((s-mp*mp)*(s-mm*mm)));
        }
        const double pi = M_PI;
        const complex<double> i =  complex<double>(0.0, 1.0);

        const unsigned N_orbital = 3;

        double beta(const double & s){
            if (s < mp*mp) { return 0.; } //Kinematic threshold
            else { return pow(sqlk(s)/s, 3.); }
        }

        complex<double> rho(const double & s) {
            complex<double> result = 0.0;
            if (s < mm*mm){
                // result += mm*sqlk(s)*(2*mm*mp*s-(-3*mp*mp*s+mm*mm*(2*mp*mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 2*pow(mp,3)*pow(mm*mm-s,2)*std::log((-2*(s+sqlk(s))+mm*mm+mp*mp)/(mp*mp-mm*mm));
                // result *= i*pow(mp*mp-s,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(mm*mm-s));
                return result;
            }
            else if (s < mp*mp){
                // result += mm*sqlk(s)*(2*mm*mp*s-(-3*mp*mp*s+mm*mm*(2*mp*mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 4*pow(mp,3)*pow(s-mm*mm,2)*std::atan(sqrt((s-mm*mm)/(mp*mp-s)));
                // result *= i*pow(mp*mp-s,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(s-mm*mm));
                return result;
            }
            else {
                // result += mm*sqlk(s)*(-2*mm*mp*s+(-3*mp*mp*s+mm*mm*(2*mp*mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 2*pow(mp,3)*pow(s-mm*mm,2)*std::log((2*(s+sqlk(s))-mm*mm-mp*mp)/(mp*mp-mm*mm));
                // result *= i*pow(s-mp*mp,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(s-mm*mm));
                result += pow(sqlk(s)/s, 3.);
                return result;
           }
        }
    };


    // V -> VP channel
    template <unsigned nchannels_, unsigned nresonances_>
    struct VPchan :
    public KMatrix<nchannels_, nresonances_>::Channel
    {

        VPchan(std::string name, double m1, double m2, unsigned N_orbital, std::vector<Parameter> g0s) :
        KMatrix<nchannels_, nresonances_>::Channel(name, m1, m2, N_orbital, g0s)
        {
        };

        //Usefull definitions for beta and rho
        const double mm = this->_m1 - this->_m2;
        const double mp = this->_m1 + this->_m2;
        //sqrt of the Källen factor, defined with an absolute value
        double sqlk(const double & s) {
            return std::sqrt(std::abs((s-mp*mp)*(s-mm*mm)));
        }
        const double pi = M_PI;
        const complex<double> i =  complex<double>(0.0, 1.0);

        const unsigned N_orbital = 3;

        double beta(const double & s){
            if (s < mp*mp) { return 0.; } //Kinematic threshold
            else { return sqlk(s)/s; }
        }

        complex<double> rho(const double & s) {
            complex<double> result = 0.0;
            if (s < mm*mm){
                // result += mm*(mp*mp-s)*std::log((mm+mp)/(mp-mm));
                // result += mp*sqlk(s)*std::log((mm*mm+mp*mp-2*s-2*sqlk(s))/(mp*mp-mm*mm));
                // result *= i/(mp*pi*s);
                return result;
            }
            else if (s < mp*mp){
                // result += mm*(mp*mp-s)*std::log((mm+mp)/(mp-mm));
                // result += 2*mp*sqlk(s)*std::atan(std::sqrt((s-mm*mm)/(mp*mp-s)));
                // result *= i/(mp*pi*s);
                return result;
            }
            else {
                // result += mm*(mp*mp-s)*std::log((mm+mp)/(mp-mm));
                // result += mp*sqlk(s)*std::log((2*s+2*sqlk(s)-mm*mm-mp*mp)/(mp*mp-mm*mm));
                // result *= i/(mp*pi*s);
                result += sqlk(s)/s;
                return result;
            }
        }
    };


    template <unsigned nchannels_, unsigned nresonances_>
    struct charmonium_resonance :
    public KMatrix<nchannels_, nresonances_>::Resonance
    {
        charmonium_resonance(std::string name, double m) :
        KMatrix<nchannels_, nresonances_>::Resonance(name, m)
        {
        };
    };

    class EEToCCBar :
        public ParameterUser,
        public PrivateImplementationPattern<EEToCCBar>
    {
        public:
            EEToCCBar(const Parameters & parameters, const Options & options);
            ~EEToCCBar();

            // psi2S widths
            double psi2S_ee_width() const;
            double psi2S_eff_width() const;
            double psi2S_total_width() const;

            // sigma(ee -> channel)
            double sigma_eetoee(const double & E) const;
            double sigma_eetoeff(const double & E) const;
            double sigma_eetoD0Dbar0(const double & E) const;
            double sigma_eetoDpDm(const double & E) const;

            // Rc ratio
            double Rc(const double & E) const;

    };

}
#endif
