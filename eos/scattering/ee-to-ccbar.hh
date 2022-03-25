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
#include <eos/utils/kmatrix.hh>
#include <eos/utils/concrete-cacheable-observable.hh>

#include <vector>

namespace eos

{
    /*
    Channels follow the following convention
    #   name          type         Nf      copy
    0   ee            PP (P)       3       -
    1   eff(2S)       PP (P)       3       -
    2   eff(3770)     PP (P)       3       -
    3   eff(4040)     PP (P)       3       -
    4   eff(4160)     PP (P)       3       -
    5   eff(4415)     PP (P)       3       -
    6   D0   D0bar    PP (P)       3       -
    7   D+   D-       PP (P)       3       6 (isospin)
    8   D0   D*0bar   VP           3       -
    9   D*0  D0bar    VP           3       8 (c.c.)
    10  D+   D*-      VP           3       8 (isospin)
    11  D*+  D-       VP           3       8 (c.c.)
    12  Ds   Ds       PP (P)       3       -
    13  D*0  D*0bar   VV (P, S=0)  3       -
    14  D*0  D*0bar   VV (P, S=2)  3       -
    15  D*0  D*0bar   VV (F, S=2)  7       -
    16  D*+  D*-      VV (P, S=0)  3       13 (isospin)
    17  D*+  D*-      VV (P, S=2)  3       14 (isospin)
    18  D*+  D*-      VV (F, S=2)  7       15 (isospin)
    19  Ds+  Ds*-     VP           3       -
    20  Ds*+ Ds-      VP           3       19 (c.c.)
    21  Ds*+ Ds*-     VV (P, S=0)  3       -
    22  Ds*+ Ds*-     VV (P, S=2)  3       -
    23  Ds*+ Ds*-     VV (F, S=2)  7       -
    */

    // Effective channel
    template <unsigned nchannels_, unsigned nresonances_>
    struct EffChannel :
    public KMatrix<nchannels_, nresonances_>::Channel
    {

        EffChannel(std::string name, double m1, double m2, std::array<Parameter, nresonances_> g0s) :
            KMatrix<nchannels_, nresonances_>::Channel(name, m1, m2, 0, g0s)
        {
        };

        // Usefull definitions for beta and rho
        const double mm = this->_m1 - this->_m2;
        const double mp = this->_m1 + this->_m2;
        // Square root of the Källen factor, defined with an absolute value
        double sqlk(const double & s) {
            return std::sqrt(std::abs((s - mp * mp) * (s - mm * mm)));
        }
        const double pi = M_PI;
        const complex<double> i =  complex<double>(0.0, 1.0);

        double beta(const double & s){
            if (s < mp * mp) { return 0.; } // Kinematic threshold
            else { return pow(sqlk(s)/s, 3.); }
        }

        complex<double> rho(const double & s) {
            complex<double> result = 0.0;
            if (s < mm * mm){
                // result += mm*sqlk(s)*(2*mm*mp*s-(-3*mp * mp*s+mm * mm*(2*mp * mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 2*pow(mp,3)*pow(mm * mm-s,2)*std::log((-2*(s+sqlk(s))+mm * mm+mp * mp)/(mp * mp-mm * mm));
                // result *= i*pow(mp * mp-s,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(mm * mm-s));
                return result;
            }
            else if (s < mp * mp){
                // result += mm*sqlk(s)*(2*mm*mp*s-(-3*mp * mp*s+mm * mm*(2*mp * mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 4*pow(mp,3)*pow(s-mm * mm,2)*std::atan(sqrt((s-mm * mm)/(mp * mp-s)));
                // result *= i*pow(mp * mp-s,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(s-mm * mm));
                return result;
            }
            else {
                // result += mm*sqlk(s)*(-2*mm*mp*s+(-3*mp * mp*s+mm * mm*(2*mp * mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 2*pow(mp,3)*pow(s-mm * mm,2)*std::log((2*(s+sqlk(s))-mm * mm-mp * mp)/(mp * mp-mm * mm));
                // result *= i*pow(s-mp * mp,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(s-mm * mm));
                result += pow(sqlk(s)/s, 3.);
                return result;
           }
        }
    };


    // V -> PP channel
    template <unsigned nchannels_, unsigned nresonances_>
    struct PWavePPChannel :
    public KMatrix<nchannels_, nresonances_>::Channel
    {

        PWavePPChannel(std::string name, double m1, double m2, std::array<Parameter, nresonances_> g0s) :
            KMatrix<nchannels_, nresonances_>::Channel(name, m1, m2, 1, g0s)
        {
        };

        // Usefull definitions for beta and rho
        const double mm = this->_m1 - this->_m2;
        const double mp = this->_m1 + this->_m2;
        // Square root of the Källen factor, defined with an absolute value
        double sqlk(const double & s) {
            return std::sqrt(std::abs((s - mp * mp) * (s - mm * mm)));
        }
        const double pi = M_PI;
        const complex<double> i =  complex<double>(0.0, 1.0);

        double beta(const double & s){
            if (s < mp * mp) { return 0.; } // Kinematic threshold
            else { return pow(sqlk(s)/s, 3.); }
        }

        complex<double> rho(const double & s) {
            complex<double> result = 0.0;
            if (s < mm * mm){
                // result += mm*sqlk(s)*(2*mm*mp*s-(-3*mp * mp*s+mm * mm*(2*mp * mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 2*pow(mp,3)*pow(mm * mm-s,2)*std::log((-2*(s+sqlk(s))+mm * mm+mp * mp)/(mp * mp-mm * mm));
                // result *= i*pow(mp * mp-s,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(mm * mm-s));
                return result;
            }
            else if (s < mp * mp){
                // result += mm*sqlk(s)*(2*mm*mp*s-(-3*mp * mp*s+mm * mm*(2*mp * mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 4*pow(mp,3)*pow(s-mm * mm,2)*std::atan(sqrt((s-mm * mm)/(mp * mp-s)));
                // result *= i*pow(mp * mp-s,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(s-mm * mm));
                return result;
            }
            else {
                // result += mm*sqlk(s)*(-2*mm*mp*s+(-3*mp * mp*s+mm * mm*(2*mp * mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 2*pow(mp,3)*pow(s-mm * mm,2)*std::log((2*(s+sqlk(s))-mm * mm-mp * mp)/(mp * mp-mm * mm));
                // result *= i*pow(s-mp * mp,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(s-mm * mm));
                result += pow(sqlk(s)/s, 3.);
                return result;
           }
        }
    };


    // V -> VP channel
    template <unsigned nchannels_, unsigned nresonances_>
    struct PWaveVPChannel :
    public KMatrix<nchannels_, nresonances_>::Channel
    {

        PWaveVPChannel(std::string name, double m1, double m2, std::array<Parameter, nresonances_> g0s) :
            KMatrix<nchannels_, nresonances_>::Channel(name, m1, m2, 1, g0s)
        {
        };

        // Usefull definitions for beta and rho
        const double mm = this->_m1 - this->_m2;
        const double mp = this->_m1 + this->_m2;
        // Square root of the Källen factor, defined with an absolute value
        double sqlk(const double & s) {
            return std::sqrt(std::abs((s - mp * mp) * (s - mm * mm)));
        }
        const double pi = M_PI;
        const complex<double> i =  complex<double>(0.0, 1.0);

        double beta(const double & s){
            if (s < mp * mp) { return 0.; } // Kinematic threshold
            else { return pow(sqlk(s)/s, 3.); }
        }

        complex<double> rho(const double & s) {
            complex<double> result = 0.0;
            if (s < mm * mm){
                // result += mm*(mp * mp-s)*std::log((mm+mp)/(mp-mm));
                // result += mp*sqlk(s)*std::log((mm * mm+mp * mp-2*s-2*sqlk(s))/(mp * mp-mm * mm));
                // result *= i/(mp*pi*s);
                return result;
            }
            else if (s < mp * mp){
                // result += mm*(mp * mp-s)*std::log((mm+mp)/(mp-mm));
                // result += 2*mp*sqlk(s)*std::atan(std::sqrt((s-mm * mm)/(mp * mp-s)));
                // result *= i/(mp*pi*s);
                return result;
            }
            else {
                // result += mm*(mp * mp-s)*std::log((mm+mp)/(mp-mm));
                // result += mp*sqlk(s)*std::log((2*s+2*sqlk(s)-mm * mm-mp * mp)/(mp * mp-mm * mm));
                // result *= i/(mp*pi*s);
                result += pow(sqlk(s)/s, 3.);
                return result;
            }
        }
    };


    // V -> VV channel
    template <unsigned nchannels_, unsigned nresonances_>
    struct PWaveVVChannel :
    public KMatrix<nchannels_, nresonances_>::Channel
    {

        PWaveVVChannel(std::string name, double m1, double m2, std::array<Parameter, nresonances_> g0s) :
            KMatrix<nchannels_, nresonances_>::Channel(name, m1, m2, 1, g0s)
        {
        };

        // Usefull definitions for beta and rho
        const double mm = this->_m1 - this->_m2;
        const double mp = this->_m1 + this->_m2;
        // Square root of the Källen factor, defined with an absolute value
        double sqlk(const double & s) {
            return std::sqrt(std::abs((s - mp * mp) * (s - mm * mm)));
        }
        const double pi = M_PI;
        const complex<double> i =  complex<double>(0.0, 1.0);

        double beta(const double & s){
            if (s < mp * mp) { return 0.; } // Kinematic threshold
            else { return pow(sqlk(s)/s, 3.); }
        }

        complex<double> rho(const double & s) {
            complex<double> result = 0.0;
            if (s < mm * mm){
                // result += mm*sqlk(s)*(2*mm*mp*s-(-3*mp * mp*s+mm * mm*(2*mp * mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 2*pow(mp,3)*pow(mm * mm-s,2)*std::log((-2*(s+sqlk(s))+mm * mm+mp * mp)/(mp * mp-mm * mm));
                // result *= i*pow(mp * mp-s,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(mm * mm-s));
                return result;
            }
            else if (s < mp * mp){
                // result += mm*sqlk(s)*(2*mm*mp*s-(-3*mp * mp*s+mm * mm*(2*mp * mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 4*pow(mp,3)*pow(s-mm * mm,2)*std::atan(sqrt((s-mm * mm)/(mp * mp-s)));
                // result *= i*pow(mp * mp-s,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(s-mm * mm));
                return result;
            }
            else {
                // result += mm*sqlk(s)*(-2*mm*mp*s+(-3*mp * mp*s+mm * mm*(2*mp * mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 2*pow(mp,3)*pow(s-mm * mm,2)*std::log((2*(s+sqlk(s))-mm * mm-mp * mp)/(mp * mp-mm * mm));
                // result *= i*pow(s-mp * mp,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(s-mm * mm));
                result += pow(sqlk(s)/s, 3.);
                return result;
            }
        }
    };

    // V -> VV channel
    template <unsigned nchannels_, unsigned nresonances_>
    struct FWaveVVChannel :
    public KMatrix<nchannels_, nresonances_>::Channel
    {

        FWaveVVChannel(std::string name, double m1, double m2, std::array<Parameter, nresonances_> g0s) :
            KMatrix<nchannels_, nresonances_>::Channel(name, m1, m2, 2, g0s)
        {
        };

        // Usefull definitions for beta and rho
        const double mm = this->_m1 - this->_m2;
        const double mp = this->_m1 + this->_m2;
        // Square root of the Källen factor, defined with an absolute value
        double sqlk(const double & s) {
            return std::sqrt(std::abs((s - mp * mp) * (s - mm * mm)));
        }
        const double pi = M_PI;
        const complex<double> i =  complex<double>(0.0, 1.0);

        double beta(const double & s){
            if (s < mp * mp) { return 0.; } // Kinematic threshold
            else { return pow(sqlk(s)/s, 7.); }
        }

        complex<double> rho(const double & s) {
            complex<double> result = 0.0;
            if (s < mm * mm){
                // result += mm*sqlk(s)*(2*mm*mp*s-(-3*mp * mp*s+mm * mm*(2*mp * mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 2*pow(mp,3)*pow(mm * mm-s,2)*std::log((-2*(s+sqlk(s))+mm * mm+mp * mp)/(mp * mp-mm * mm));
                // result *= i*pow(mp * mp-s,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(mm * mm-s));
                return result;
            }
            else if (s < mp * mp){
                // result += mm*sqlk(s)*(2*mm*mp*s-(-3*mp * mp*s+mm * mm*(2*mp * mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 4*pow(mp,3)*pow(s-mm * mm,2)*std::atan(sqrt((s-mm * mm)/(mp * mp-s)));
                // result *= i*pow(mp * mp-s,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(s-mm * mm));
                return result;
            }
            else {
                // result += mm*sqlk(s)*(-2*mm*mp*s+(-3*mp * mp*s+mm * mm*(2*mp * mp+s))*std::log((mp+mm)/(mp-mm)));
                // result += 2*pow(mp,3)*pow(s-mm * mm,2)*std::log((2*(s+sqlk(s))-mm * mm-mp * mp)/(mp * mp-mm * mm));
                // result *= i*pow(s-mp * mp,1.5)/(2*pow(mp,3)*pi*s*s*sqrt(s-mm * mm));
                result += pow(sqlk(s)/s, 7.);
                return result;
            }
        }
    };


    template <unsigned nchannels_, unsigned nresonances_>
    struct CharmoniumResonance :
    public KMatrix<nchannels_, nresonances_>::Resonance
    {
        CharmoniumResonance(std::string name, Parameter m, Parameter q) :
        KMatrix<nchannels_, nresonances_>::Resonance(name, m, q)
        {
        };
    };

    class EEToCCBar :
        public ParameterUser,
        public PrivateImplementationPattern<EEToCCBar>
    {
        public:

            const static long unsigned nchannels = 24;
            const static long unsigned nresonances = 5;

            struct IntermediateResult :
                public CacheableObservable::IntermediateResult
            {
                std::shared_ptr<KMatrix<nchannels, nresonances>> K;
                std::array<complex<double>, nchannels> tmatrix_row_0;

                double E;
                double s;
            };

            EEToCCBar(const Parameters & parameters, const Options & options);
            ~EEToCCBar();

            // Observables
            const IntermediateResult * prepare(const double & E) const;

            // double evaluate(const IntermediateResult *) const;

            // psi2S widths
            double psi2S_ee_width() const;
            double psi2S_eff_width() const;
            double psi2S_total_width() const;
            double psi3770_D0Dbar0_width() const;
            double psi3770_DpDm_width() const;
            double psi3770_eff_width() const;
            double psi3770_total_width() const;
            double psi4040_total_width() const;
            double psi4160_total_width() const;
            double psi4415_total_width() const;
            double psi4040_DD_width() const;
            double psi4040_DDst_width() const;
            double psi4040_DstDst_width() const;
            double psi4160_DD_width() const;
            double psi4160_DDst_width() const;
            double psi4160_DstDst_width() const;
            double psi4415_DD_width() const;
            double psi4415_DDst_width() const;
            double psi4415_DstDst_width() const;

            // sigma(ee -> channel)
            double sigma_eetoee(const IntermediateResult *) const;
            double sigma_eetoeff(const IntermediateResult *) const;
            double sigma_eetoD0Dbar0(const IntermediateResult *) const;
            double sigma_eetoDpDm(const IntermediateResult *) const;
            double sigma_eetoD0Dbarst0(const IntermediateResult *) const;
            double sigma_eetoDpDstm(const IntermediateResult *) const;
            double sigma_eetoDspDsm(const IntermediateResult *) const;
            double sigma_eetoDst0Dbarst0(const IntermediateResult *) const;
            double sigma_eetoDstpDstm(const IntermediateResult *) const;
            double sigma_eetoDstpTDstmT(const IntermediateResult *) const;
            double sigma_eetoDstpTDstmL(const IntermediateResult *) const;
            double sigma_eetoDstpLDstmL(const IntermediateResult *) const;
            double sigma_eetoDspDsstm(const IntermediateResult *) const;
            double sigma_eetoDsstpDsstm(const IntermediateResult *) const;

            // R ratios
            double R(const IntermediateResult *) const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();

    };
}
#endif
