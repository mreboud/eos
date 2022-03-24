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

#include <eos/scattering/ee-to-ccbar.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/kmatrix-impl.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <vector>
#include <array>
#include <memory>

namespace eos
{

    template <>
    struct Implementation<EEToCCBar>
    {
        /* Channels and resonances
        Up to DD*       8  channels 2 resonances
        Up to D*D*      14 channels 3 resonances
        Up to 4.8 GeV   24 channels 5 resonances
         */
        const static long unsigned nchannels = 24;
        const static long unsigned nresonances = 5;

        UsedParameter hbar;
        UsedParameter alpha_em;
        UsedParameter m_e;
        UsedParameter m_D0;
        UsedParameter m_D;
        UsedParameter m_Dst0;
        UsedParameter m_Dst;
        UsedParameter m_Ds;
        UsedParameter m_Dsst;

        static const std::vector<std::string> resonance_names;
        static const std::vector<std::string> channel_names;

        //Charmonium masses and sizes
        std::array<UsedParameter, nresonances> m;
        std::array<UsedParameter, nresonances> q;

        // Channel-Resonance couplings
        std::array<std::array<UsedParameter, nchannels>, nresonances> g0;

        // Non-cc contribution to the Rc ratio
        std::array<std::array<UsedParameter, nchannels>, nchannels> bkgcst;

        // R_uds
        UsedParameter Rconstant;

        // Normalization of the exclusive channels
        UsedParameter exclusive_norm;

        std::shared_ptr<KMatrix<nchannels, nresonances>> K;

        using IntermediateResult = EEToCCBar::IntermediateResult;
        IntermediateResult _intermediate_result;

        template <typename T, T... indices>
        auto _resonance_masses(const Parameters & p, ParameterUser & u, std::integer_sequence<T, indices...>)
            -> std::array<UsedParameter, sizeof...(indices)>
        {
            if (sizeof...(indices) > resonance_names.size())
                throw InternalError("The number of requested resonances is larger than the number of resonance names.");

            return std::array<UsedParameter, sizeof...(indices)>
            {{
                UsedParameter(p["mass::" + resonance_names[indices]], u)...
            }};
        }

        template <typename T, T... indices>
        auto _resonance_sizes(const Parameters & p, ParameterUser & u, std::integer_sequence<T, indices...>)
            -> std::array<UsedParameter, sizeof...(indices)>
        {
            if (sizeof...(indices) > resonance_names.size())
                throw InternalError("The number of requested resonances is larger than the number of resonance names.");

            return std::array<UsedParameter, sizeof...(indices)>
            {{
                UsedParameter(p["size::" + resonance_names[indices]], u)...
            }};
        }

        long unsigned _filter_channel_index(unsigned channel_index)
        {
            switch (channel_index)
            {
                case 7: // D+D-
                    return 6; // D0Dbar0

                case 9: // D*0D0bar
                case 10: // D+D*-
                case 11: // D*+D-
                    return 8; // D0D*0bar

                case 15: // D*0D*0bar S=2 F wave
                    return 14; // D*0D*0bar S=2 P wave

                case 16: // D*+D*- S=0
                    return 13; // D*0D*0bar S=0

                case 17: // D*+D*- S=2 P wave
                case 18: // D*+D*- S=2 F wave
                    return 14; // D*0D*0bar S=2

                case 20: // Ds*+Ds-
                    return 19; // Ds+Ds*-

                case 23: // Ds*+Ds*- S=2 F wave
                    return 22; // Ds+Ds*- S=2 P wave

                default:
                    return channel_index;
            }
        }

        template <typename T, T... column_indices>
        auto _g0_row(const Parameters & p, ParameterUser & u, T row_index, std::integer_sequence<T, column_indices...>)
            -> std::array<UsedParameter, sizeof...(column_indices)>
        {
            if (sizeof...(column_indices) > channel_names.size())
                throw InternalError("The number of requested channels is larger than the number of channel names.");

            return std::array<UsedParameter, sizeof...(column_indices)>
            {{
                UsedParameter(p["ee->ccbar::g0(" + resonance_names[row_index] + "," + channel_names[_filter_channel_index(column_indices)] + ")"], u)...
            }};
        }

        template <typename T, T... row_indices, T... column_indices>
        auto _g0_matrix(const Parameters & p, ParameterUser & u, std::integer_sequence<T, row_indices...>, std::integer_sequence<T, column_indices...> column_seq)
            -> std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>
        {
            if (sizeof...(row_indices) > resonance_names.size())
                throw InternalError("The number of requested resonances is larger than the number of resonance names.");

            return std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>
            {{
                _g0_row(p, u, row_indices, column_seq)...
            }};
        }

        template <typename T, T... column_indices>
        auto _c_row(const Parameters & p, ParameterUser & u, T row_index, std::integer_sequence<T, column_indices...>)
            -> std::array<UsedParameter, sizeof...(column_indices)>
        {
            if (sizeof...(column_indices) > channel_names.size())
                throw InternalError("The number of requested channels (column) is larger than the number of channel names.");

            return std::array<UsedParameter, sizeof...(column_indices)>
            {{
                UsedParameter(p["ee->ccbar::c(" + channel_names[_filter_channel_index(row_index)] + ","
                        + channel_names[_filter_channel_index(column_indices)] + ")"], u)...
            }};
        }

        template <typename T, T... row_indices, T... column_indices>
        auto _c_matrix(const Parameters & p, ParameterUser & u, std::integer_sequence<T, row_indices...>, std::integer_sequence<T, column_indices...> column_seq)
            -> std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>
        {
            if (sizeof...(row_indices) > channel_names.size())
                throw InternalError("The number of requested channels (row) is larger than the number of channel names.");

            return std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>
            {{
                _c_row(p, u, row_indices, column_seq)...
            }};
        }

        template <typename T, T... row_indices>
        auto _get_g0_column(T column_index, std::integer_sequence<T, row_indices...>)
            -> std::array<Parameter, sizeof...(row_indices)>
        {
            return std::array<Parameter, sizeof...(row_indices)>
            {{
                g0[row_indices][column_index]...
            }};
        }

        template <typename T, T... column_indices>
        auto _get_bkgcst_row(T row_index, std::integer_sequence<T, column_indices...>)
            -> std::array<Parameter, sizeof...(column_indices)>
        {
            return std::array<Parameter, sizeof...(column_indices)>
            {{
                bkgcst[row_index][column_indices]...
            }};
        }

        template <typename T, T... row_indices, T... column_indices>
        auto _get_bkgcst(std::integer_sequence<T, row_indices...>, std::integer_sequence<T, column_indices...> column_seq)
            -> std::array<std::array<Parameter, sizeof...(column_indices)>, sizeof...(row_indices)>
        {
            return std::array<std::array<Parameter, sizeof...(column_indices)>, sizeof...(row_indices)>
            {{
                _get_bkgcst_row(row_indices, column_seq)...
            }};
        }

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            hbar(p["QM::hbar"], u),
            alpha_em(p["QED::alpha_e(0)"], u),
            m_e(p["mass::e"], u),
            m_D0(p["mass::D^0"], u),
            m_D(p["mass::D^+"], u),
            m_Dst0(p["mass::D_u^*"], u),
            m_Dst(p["mass::D_d^*"], u),
            m_Ds(p["mass::D_s"], u),
            m_Dsst(p["mass::D_s^*"], u),

            m(_resonance_masses(p, u, std::make_index_sequence<nresonances>())),
            q(_resonance_sizes(p, u, std::make_index_sequence<nresonances>())),
            g0(_g0_matrix(p, u, std::make_index_sequence<nresonances>(), std::make_index_sequence<nchannels>())),
            bkgcst(_c_matrix(p, u, std::make_index_sequence<nchannels>(), std::make_index_sequence<nchannels>())),

            Rconstant(p["ee->ccbar::Rconstant"], u),
            exclusive_norm(p["ee->ccbar::exclusive_norm"], u)
        {
            std::array<std::shared_ptr<KMatrix<nchannels, nresonances>::Resonance>, nresonances> resonance_array;
            for (unsigned i = 0; i<nresonances; i++)
            {
                resonance_array[i] = std::make_shared<charmonium_resonance<nchannels, nresonances>>(resonance_names[i] + "_res", m[i], q[i]);
            }

            std::array<std::shared_ptr<KMatrix<nchannels, nresonances>::Channel>, nchannels> channel_array;
            for (unsigned i = 0; i<nchannels; i++)
            {
                switch (i)
                {
                    case 0: // ee channel
                    case 1: // eff channels
                    case 2:
                    case 3:
                    case 4:
                    case 5:
                        channel_array[i] = std::make_shared<PP_Pwavechan<nchannels, nresonances>>(channel_names[i] + "_chan", m_e, m_e, 1, _get_g0_column(_filter_channel_index(i), std::make_index_sequence<nresonances>()));
                        break;
                    case 6:
                        channel_array[i] = std::make_shared<PP_Pwavechan<nchannels, nresonances>>(channel_names[i] + "_chan", m_D0, m_D0, 1, _get_g0_column(_filter_channel_index(i), std::make_index_sequence<nresonances>()));
                        break;
                    case 7:
                        channel_array[i] = std::make_shared<PP_Pwavechan<nchannels, nresonances>>(channel_names[i] + "_chan", m_D, m_D, 1, _get_g0_column(_filter_channel_index(i), std::make_index_sequence<nresonances>()));
                        break;
                    case 8:
                    case 9: //D0Dbarst0
                        channel_array[i] = std::make_shared<VP_Pwavechan<nchannels, nresonances>>(channel_names[i] + "_chan", m_D0, m_Dst0, 1, _get_g0_column(_filter_channel_index(i), std::make_index_sequence<nresonances>()));
                        break;
                    case 10:
                    case 11: //DpDstm
                        channel_array[i] = std::make_shared<VP_Pwavechan<nchannels, nresonances>>(channel_names[i] + "_chan", m_D, m_Dst, 1, _get_g0_column(_filter_channel_index(i), std::make_index_sequence<nresonances>()));
                        break;
                    case 12:
                        channel_array[i] = std::make_shared<PP_Pwavechan<nchannels, nresonances>>(channel_names[i] + "_chan", m_Ds, m_Ds, 1, _get_g0_column(_filter_channel_index(i), std::make_index_sequence<nresonances>()));
                        break;
                    case 13:
                    case 14:
                        channel_array[i] = std::make_shared<VV_Pwavechan<nchannels, nresonances>>(channel_names[i] + "_chan", m_Dst0, m_Dst0, 1, _get_g0_column(_filter_channel_index(i), std::make_index_sequence<nresonances>()));
                        break;
                    case 15:
                        channel_array[i] = std::make_shared<VV_Fwavechan<nchannels, nresonances>>(channel_names[i] + "_chan", m_Dst0, m_Dst0, 3, _get_g0_column(_filter_channel_index(i), std::make_index_sequence<nresonances>()));
                        break;
                    case 16:
                    case 17:
                        channel_array[i] = std::make_shared<VV_Pwavechan<nchannels, nresonances>>(channel_names[i] + "_chan", m_Dst, m_Dst, 1, _get_g0_column(_filter_channel_index(i), std::make_index_sequence<nresonances>()));
                        break;
                    case 18:
                        channel_array[i] = std::make_shared<VV_Fwavechan<nchannels, nresonances>>(channel_names[i] + "_chan", m_Dst, m_Dst, 3, _get_g0_column(_filter_channel_index(i), std::make_index_sequence<nresonances>()));
                        break;
                    case 19:
                    case 20:
                        channel_array[i] = std::make_shared<VP_Pwavechan<nchannels, nresonances>>(channel_names[i] + "_chan", m_Ds, m_Dsst, 1, _get_g0_column(_filter_channel_index(i), std::make_index_sequence<nresonances>()));
                        break;
                    case 21:
                    case 22:
                        channel_array[i] = std::make_shared<VV_Pwavechan<nchannels, nresonances>>(channel_names[i] + "_chan", m_Dsst, m_Dsst, 1, _get_g0_column(_filter_channel_index(i), std::make_index_sequence<nresonances>()));
                        break;
                    case 23:
                        channel_array[i] = std::make_shared<VV_Fwavechan<nchannels, nresonances>>(channel_names[i] + "_chan", m_Dsst, m_Dsst, 3, _get_g0_column(_filter_channel_index(i), std::make_index_sequence<nresonances>()));
                        break;

                    default:
                        throw InternalError("The number of requested channels (array) is larger than the number of channel names.");
                }
            }

            K = std::shared_ptr<KMatrix<nchannels, nresonances>> (
                new KMatrix<nchannels, nresonances>(
                    channel_array,
                    resonance_array,
                    _get_bkgcst(std::make_index_sequence<nchannels>(),std::make_index_sequence<nchannels>()),
                    "KMatrix"
                    )
                );
        }

        const IntermediateResult * prepare(const double & E)
        {
            _intermediate_result.tmatrix_row_0 = K->tmatrix_row(0, E*E);

            _intermediate_result.E = E;
            _intermediate_result.s = E*E;

            return &_intermediate_result;
        }


        inline double sigma_eetomumu(const double & E)
        {
            // Conversion factor between GeV^2 and nb
            const double speedoflight = 299792458.; //Exact value
            const double GeVtonb = 10 * pow( 1.e18 * hbar * speedoflight, 2.0);

            return GeVtonb * 4.0 * M_PI * alpha_em*alpha_em / (3.0 * E*E);
        }

        // amplitude of ee -> channel
        complex<double> amplitude_eetochannel(const IntermediateResult * intermediate_result, const unsigned & channel)
        {
            // Conversion factor between GeV^2 and nb
            const double speedoflight = 299792458.; //Exact value
            const double GeVtonb = 10 * pow( 1.e18 * hbar * speedoflight, 2.0);

            // Channel properties
            const double Nf = 2 * K->_channels[channel]->_l_orbital + 1;
            const double rhof = real(K->_channels[channel]->rho(intermediate_result->s));

            // Get T-matrix[ee, channel]
            const complex<double> T1f = intermediate_result->tmatrix_row_0[channel];

            return sqrt(GeVtonb * 16. * M_PI / intermediate_result->s * Nf * rhof) * T1f;
        }

        double sigma_eetochannel(const IntermediateResult * intermediate_result, const unsigned & channel)
        {
            return norm(amplitude_eetochannel(intermediate_result, channel));
        }

        // Widths
        double res_partial_width(unsigned res, unsigned channel)
        {
            return K->partial_width(res, channel);
        }

        double res_total_width(unsigned res)
        {
            return K->width(res);
        }

        // Ratios
        double R(const IntermediateResult * intermediate_result)
        {
            double total_xsec = 0.0;

            for (unsigned i = 1; i < nchannels; i++)
            {
                total_xsec += sigma_eetochannel(intermediate_result, i);
            }

            return total_xsec / sigma_eetomumu(intermediate_result->E) + Rconstant; //Add constant term
        }

        double Rc(const IntermediateResult * intermediate_result)
        {
            double total_xsec = 0.0;

            for (unsigned i = 6; i < nchannels; i++)
            {
                total_xsec += sigma_eetochannel(intermediate_result, i);
            }

            return total_xsec / sigma_eetomumu(intermediate_result->E);
        }

        // Rudsc constraint
        double Rudsc_prior(const IntermediateResult * intermediate_result)
        {
            const double value = R(intermediate_result)/3.6; //th value taken from Topsy-Turvy paper

            if (value < 0.0)
            {
                throw InternalError("R ratio was found to be negative!");
            }
            else if ((0.0 <= value) && (value < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.036; // 1% uncertainty, to be discussed
                return -pow((value - 1.0) / sigma, 2) / 2.0;
            }
        }

        // psi(4040) width constraint
        double psi4040_total_width_prior()
        {
            const double value = res_total_width(2)/80; //th value taken from Topsy-Turvy paper

            if (value < 0.0)
            {
                throw InternalError("psi(4040) total width was found to be negative!");
            }
            else if ((0.0 <= value) && (value < 1.0))
            {
                return 0.0;
            }
            else
            {
                // add an r-fit like penalty
                static const double sigma = 0.01; //small value
                return -pow((value - 1.0) / sigma, 2) / 2.0;
            }
        }
    };

    const std::vector<std::string>
    Implementation<EEToCCBar>::resonance_names = { "psi(2S)", "psi(3770)", "psi(4040)", "psi(4160)", "psi(4415)" };

    const std::vector<std::string>
    Implementation<EEToCCBar>::channel_names = { "ee", "eff(2S)", "eff(3770)", "eff(4040)", "eff(4160)", "eff(4415)",
            "D^0Dbar^0", "D^+D^-", "D^0Dbar^*0", "D^*0Dbar^0", "D^+D^*-", "D^*+D^-", "D_s^+D_s^-", "D^*0Dbar^*0P0",
            "D^*0Dbar^*0P2", "D^*0Dbar^*0F2", "D^*+D^*-P0", "D^*+D^*-P2", "D^*+D^*-F2", "D_s^+D_s^*-", "D_s^*+D_s^-",
            "D_s^*+D_s^*-P0", "D_s^*+D_s^*-P2", "D_s^*+D_s^*-F2" };


    EEToCCBar::EEToCCBar(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<EEToCCBar>(new Implementation<EEToCCBar>(parameters, options, *this))
    {
    }

    EEToCCBar::~EEToCCBar()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<EEToCCBar>::options
    {
    };

    const EEToCCBar::IntermediateResult *
    EEToCCBar::prepare(const double & E) const
    {
        return _imp->prepare(E);
    }

    double
    EEToCCBar::psi2S_ee_width() const
    {
        return _imp->res_partial_width(0, 0);
    }

    double
    EEToCCBar::psi2S_eff_width() const
    {
        return _imp->res_partial_width(0, 1);
    }

    double
    EEToCCBar::psi2S_total_width() const
    {
        return _imp->res_total_width(0);
    }

    double
    EEToCCBar::psi3770_total_width() const
    {
        return _imp->res_total_width(1);
    }

    double
    EEToCCBar::psi4040_total_width() const
    {
        return _imp->res_total_width(2);
    }

    double
    EEToCCBar::psi4160_total_width() const
    {
        return _imp->res_total_width(3);
    }

    double
    EEToCCBar::psi4415_total_width() const
    {
        return _imp->res_total_width(4);
    }


    double
    EEToCCBar::sigma_eetoee(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, 0);
    }

    double
    EEToCCBar::sigma_eetoeff(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (_imp->sigma_eetochannel(ir, 1) + _imp->sigma_eetochannel(ir, 2)
        + _imp->sigma_eetochannel(ir, 3) + _imp->sigma_eetochannel(ir, 4) + _imp->sigma_eetochannel(ir, 5));
    }

    double
    EEToCCBar::sigma_eetoD0Dbar0(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, 6);
    }

    double
    EEToCCBar::sigma_eetoDpDm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, 7);
    }

    double
    EEToCCBar::sigma_eetoD0Dbarst0(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (_imp->sigma_eetochannel(ir, 8) + _imp->sigma_eetochannel(ir, 9));
    }

    double
    EEToCCBar::sigma_eetoDpDstm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (_imp->sigma_eetochannel(ir, 10) + _imp->sigma_eetochannel(ir, 11));
    }

    double
    EEToCCBar::sigma_eetoDspDsm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, 12);
    }

    double
    EEToCCBar::sigma_eetoDst0Dbarst0(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (_imp->sigma_eetochannel(ir, 13) + _imp->sigma_eetochannel(ir, 14) + _imp->sigma_eetochannel(ir, 15));
    }

    double
    EEToCCBar::sigma_eetoDstpDstm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (_imp->sigma_eetochannel(ir, 16) + _imp->sigma_eetochannel(ir, 17) + _imp->sigma_eetochannel(ir, 18));
    }

    double
    EEToCCBar::sigma_eetoDstpLDstmL(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (1.0 / 75.0) * norm(
                5.0 * _imp->amplitude_eetochannel(ir, 16) - pow(30.0, 0.5) * _imp->amplitude_eetochannel(ir, 18)
                + 2. * pow(5.0, 0.5) * _imp->amplitude_eetochannel(ir, 17)
            );
    }

    double
    EEToCCBar::sigma_eetoDstpTDstmL(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (1.0 / 5.0) * norm(
                pow(2.0, 0.5) * _imp->amplitude_eetochannel(ir, 18)
                + pow(3.0, 0.5) * _imp->amplitude_eetochannel(ir, 17)
            );
    }

    double
    EEToCCBar::sigma_eetoDstpTDstmT(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (1.0 / 150.0) * norm(
                10.0 * _imp->amplitude_eetochannel(ir, 16) + pow(30.0, 0.5) * _imp->amplitude_eetochannel(ir, 18)
                - 2. * pow(5.0, 0.5) * _imp->amplitude_eetochannel(ir, 17)
            );
    }

    double
    EEToCCBar::sigma_eetoDspDsstm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (_imp->sigma_eetochannel(ir, 19) + _imp->sigma_eetochannel(ir, 20));
    }

    double
    EEToCCBar::sigma_eetoDsstpDsstm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (_imp->sigma_eetochannel(ir, 21) + _imp->sigma_eetochannel(ir, 22) + _imp->sigma_eetochannel(ir, 23));
    }

    double
    EEToCCBar::R(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->R(ir);
    }

    double
    EEToCCBar::Rc(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->Rc(ir);
    }

    double
    EEToCCBar::Rudsc_prior(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->Rudsc_prior(ir);
    }

    double
    EEToCCBar::psi4040_total_width_prior() const
    {
        return _imp->psi4040_total_width_prior();
    }

    const std::set<ReferenceName>
    EEToCCBar::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    EEToCCBar::begin_options()
    {
        return Implementation<EEToCCBar>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    EEToCCBar::end_options()
    {
        return Implementation<EEToCCBar>::options.cend();
    }
}
