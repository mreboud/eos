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
    1   effJpsi       PP (P)       3       -
    2   eff(2S)       PP (P)       3       -
    3   eff(3770)     PP (P)       3       -
    4   D0   D0bar    PP (P)       3       -
    5   D+   D-       PP (P)       3       4 (isospin)
    */

    // Effective channel
    template <unsigned nchannels_, unsigned nresonances_, unsigned order_>
    struct EffChannel :
    public KMatrix<nchannels_, nresonances_, order_>::Channel
    {

        EffChannel(std::string name, double m1, double m2, std::array<Parameter, nresonances_> g0s) :
            KMatrix<nchannels_, nresonances_, order_>::Channel(name, m1, m2, 0, g0s)
        {
        };

        // Usefull definitions for beta and rho
        const double mm = this->_m1 - this->_m2;
        const double mp = this->_m1 + this->_m2;
        // Square root of the Källen factor, defined with an absolute value
        double sqlk(const double & s)
        {
            return std::sqrt(std::abs((s - mp * mp) * (s - mm * mm)));
        }
        const double pi = M_PI;
        const complex<double> i =  complex<double>(0.0, 1.0);

        double beta(const double & s)
        {
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
    template <unsigned nchannels_, unsigned nresonances_, unsigned order_>
    struct PWavePPChannel :
    public KMatrix<nchannels_, nresonances_, order_>::Channel
    {

        PWavePPChannel(std::string name, double m1, double m2, std::array<Parameter, nresonances_> g0s) :
            KMatrix<nchannels_, nresonances_, order_>::Channel(name, m1, m2, 1, g0s)
        {
        };

        // Usefull definitions for beta and rho
        const double mm = this->_m1 - this->_m2;
        const double mp = this->_m1 + this->_m2;
        // Square root of the Källen factor, defined with an absolute value
        double sqlk(const double & s)
        {
            return std::sqrt(std::abs((s - mp * mp) * (s - mm * mm)));
        }
        const double pi = M_PI;
        const complex<double> i =  complex<double>(0.0, 1.0);

        double beta(const double & s)
        {
            if (s < mp * mp) { return 0.; } // Kinematic threshold
            else { return pow(sqlk(s)/s, 3.); }
        }

        // complex<double> rho(const double & s) {
        //     complex<double> result = 0.0;
        //     if (s < mm * mm)
        //     {
        //         return result;
        //     }
        //     else if (s < mp * mp)
        //     {
        //         return result;
        //     }
        //     else
        //     {
        //         result += pow(sqlk(s)/s, 3.);
        //         return result;
        //     }

        complex<double> rho(const double & s) {
            if (mm != 0)
            {
                eos::InternalError("Not implemented");
            }
            complex<double> result = 0.0;
            if (s < 0)
            {
                result = std::atanh(sqlk(s) / s) + std::atanh((mp * mp / 2 - s) / sqlk(s));
                return - 0.5 * i * power_of<3>(sqlk(s)) / pi / s / s * result;
            }
            if (s < mp * mp)
            {
                result += std::atan(sqlk(s) / s) + std::atan((s - mp * mp / 2) / sqlk(s));
                return - 0.5 * i * power_of<3>(sqlk(s)) / pi / s / s * result;
            }
            else
            {
                result += pi / 2 - i * std::atanh(sqlk(s) / s) + i * std::atanh(sqlk(s) / (s - mp * mp / 2));
                return 0.5 * power_of<3>(sqlk(s)) / pi / s / s * result;
           }
        }
    };


    template <unsigned nchannels_, unsigned nresonances_, unsigned order_>
    struct CharmoniumResonance :
    public KMatrix<nchannels_, nresonances_, order_>::Resonance
    {
        CharmoniumResonance(std::string name, Parameter m, Parameter q) :
        KMatrix<nchannels_, nresonances_, order_>::Resonance(name, m, q)
        {
        };
    };

    class EEToCCBar :
        public ParameterUser,
        public PrivateImplementationPattern<EEToCCBar>
    {
        public:

            const static long unsigned nchannels = 7;
            const static long unsigned nresonances = 3;

            const static long unsigned order = 1;

            struct IntermediateResult :
                public CacheableObservable::IntermediateResult
            {
                std::shared_ptr<KMatrix<nchannels, nresonances, order>> K;
                std::array<complex<double>, nchannels> tmatrix_row_0;

                double E;
                double s;
            };

            EEToCCBar(const Parameters & parameters, const Options & options);
            ~EEToCCBar();

            // Observables
            const IntermediateResult * prepare(const double & E) const;

            // double evaluate(const IntermediateResult *) const;

            // resonances widths
            double Jpsi_ee_width() const;
            double Jpsi_eff_width() const;
            double Jpsi_total_width() const;
            double psi2S_ee_width() const;
            double psi2S_eff_width() const;
            double psi2S_Jpsipipi_width() const;
            double psi2S_total_width() const;
            double psi3770_D0Dbar0_width() const;
            double psi3770_DpDm_width() const;
            double psi3770_eff_width() const;
            double psi3770_Jpsipipi_width() const;
            double psi3770_total_width() const;

            // sigma(ee -> channel)
            double sigma_eetoee(const IntermediateResult *) const;
            double sigma_eetoeff(const IntermediateResult *) const;
            double sigma_eetoJpsipipi(const IntermediateResult *) const;
            double sigma_eetoD0Dbar0(const IntermediateResult *) const;
            double sigma_eetoDpDm(const IntermediateResult *) const;

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
