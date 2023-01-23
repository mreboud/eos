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

#include <eos/scattering/ee-to-ccbar-full.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/destringify.hh>
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
    struct Implementation<EEToCCBarFull>
    {
        UsedParameter hbar;
        UsedParameter alpha_em;
        UsedParameter m_e;
        UsedParameter m_D0;
        UsedParameter m_D;
        UsedParameter m_Dst0;
        UsedParameter m_Dst;
        UsedParameter m_Ds;
        UsedParameter m_Dsst;

        bool assume_isospin;

        static const inline std::vector<std::string> resonance_names =
        {
            "J/psi", "psi(2S)", "psi(3770)", "psi(4040)", "psi(4160)", "psi(4415)"
        };

        enum FullResonances
        {
            Jpsi = 0, psi2S, psi3770, psi4040, psi4160, psi4415
        };

        static const inline std::vector<std::string> channel_names =
        {
            "e^+e^-", "eff(Jpsi)", "eff(2S)", "D^0Dbar^0", "D^+D^-",
            "eff(3770)", "D^0Dbar^*0", "D^*0Dbar^0", "D^+D^*-", "D^*+D^-", "D_s^+D_s^-",
            "D^*0Dbar^*0P0", "D^*0Dbar^*0P2", "D^*0Dbar^*0F2", "D^*+D^*-P0", "D^*+D^*-P2", "D^*+D^*-F2",
            "eff(4040)", "D_s^+D_s^*-", "D_s^*+D_s^-",
            "eff(4160)", "D_s^*+D_s^*-P0", "D_s^*+D_s^*-P2", "D_s^*+D_s^*-F2", "eff(4415)"
        };

        enum FullChannels
        {
            ee = 0, effJpsi, eff2S, D0Dbar0, DpDm,
            eff3770, D0Dbarst0, Dst0Dbar0, DpDstm, DstpDm, DspDsm,
            Dst0Dbarst0P0, Dst0Dbarst0P2, Dst0Dbarst0F2, DstpDstmP0, DstpDstmP2, DstpDstmF2,
            eff4040, DspDsstm, DsstpDsm,
            eff4160, DsstpDsstmP0, DsstpDsstmP2, DsstpDsstmF2, eff4415
        };

        // Charmonium masses and effective sizes
        std::array<UsedParameter, EEToCCBarFull::nresonances> m;
        std::array<UsedParameter, EEToCCBarFull::nresonances> q;

        // Channel-Resonance couplings
        std::array<std::array<UsedParameter, EEToCCBarFull::nchannels>, EEToCCBarFull::nresonances> g0;

        // Non-cc contribution to the Rc ratio
        std::array<std::array<std::array<UsedParameter, EEToCCBarFull::nchannels>, EEToCCBarFull::nchannels>, EEToCCBarFull::order + 1> bkgpol;

        // R_uds
        UsedParameter Rconstant;

        // Normalization of the exclusive channels
        UsedParameter exclusive_norm;

        std::shared_ptr<KMatrix<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>> K;

        using IntermediateResult = EEToCCBarFull::IntermediateResult;
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
                UsedParameter(p["q_R::" + resonance_names[indices]], u)...
            }};
        }

        long unsigned _filter_channel_index(FullChannels channel)
        {
            if (assume_isospin)
            {
                switch (channel)
                {
                    case DpDm:
                        return D0Dbar0;

                    case Dst0Dbar0:
                    case DpDstm:
                    case DstpDm:
                        return D0Dbarst0;

                    case DstpDstmP0:
                        return Dst0Dbarst0P0;

                    case DstpDstmP2:
                        return Dst0Dbarst0P2;

                    case DstpDstmF2:
                        return Dst0Dbarst0F2;

                    case DsstpDsm:
                        return DspDsstm;

                    default:
                        return channel;
                }
            }
            else
            {
                switch (channel)
                {
                    case Dst0Dbar0:
                        return D0Dbarst0;

                    case DstpDm:
                        return DpDstm;

                    case DsstpDsm:
                        return DspDsstm;

                    default:
                        return channel;
                }
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
                UsedParameter(p["ee->ccbar::g0(" + resonance_names[row_index] + "," + channel_names[_filter_channel_index(FullChannels(column_indices))] + ")"], u)...
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
        auto _c_row(const Parameters & p, ParameterUser & u, T power, T row_index, std::integer_sequence<T, column_indices...>)
            -> std::array<UsedParameter, sizeof...(column_indices)>
        {
            if (sizeof...(column_indices) > channel_names.size())
                throw InternalError("The number of requested channels (column) is larger than the number of channel names.");

            return std::array<UsedParameter, sizeof...(column_indices)>
            {{
                UsedParameter(p["ee->ccbar::c_" + std::to_string(power) + "(" + channel_names[_filter_channel_index(FullChannels(row_index))] + ","
                        + channel_names[_filter_channel_index(FullChannels(column_indices))] + ")"], u)...
            }};
        }

        template <typename T, T... row_indices, T... column_indices>
        auto _c_matrix(const Parameters & p, ParameterUser & u, T power, std::integer_sequence<T, row_indices...>, std::integer_sequence<T, column_indices...> column_seq)
            -> std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>
        {
            if (sizeof...(row_indices) > channel_names.size())
                throw InternalError("The number of requested channels (row) is larger than the number of channel names.");

            return std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>
            {{
                _c_row(p, u, power, row_indices, column_seq)...
            }};
        }

        template <typename T, T... row_indices, T... column_indices, T... power_indices>
        auto _c_array(const Parameters & p, ParameterUser & u, std::integer_sequence<T, power_indices...>, std::integer_sequence<T, row_indices...> row_seq, std::integer_sequence<T, column_indices...> column_seq)
            -> std::array<std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>, sizeof...(power_indices)>
        {
            if (sizeof...(power_indices) >EEToCCBarFull::order + 1)
                throw InternalError("The number of background constants matrices is larger than the requested order.");

            return std::array<std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>, sizeof...(power_indices)>
            {{
                _c_matrix(p, u, power_indices, row_seq, column_seq)...
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
        auto _get_bkgpol_row(T power, T row_index, std::integer_sequence<T, column_indices...>)
            -> std::array<Parameter, sizeof...(column_indices)>
        {
            return std::array<Parameter, sizeof...(column_indices)>
            {{
                bkgpol[power][row_index][column_indices]...
            }};
        }

        template <typename T, T... row_indices, T... column_indices>
        auto _get_bkgpol_matrix(T power, std::integer_sequence<T, row_indices...>, std::integer_sequence<T, column_indices...> column_seq)
            -> std::array<std::array<Parameter, sizeof...(column_indices)>, sizeof...(row_indices)>
        {
            return std::array<std::array<Parameter, sizeof...(column_indices)>, sizeof...(row_indices)>
            {{
                _get_bkgpol_row(power, row_indices, column_seq)...
            }};
        }

        template <typename T, T... row_indices, T... column_indices, T... power_indices>
        auto _get_bkgpol(std::integer_sequence<T, power_indices...>, std::integer_sequence<T, row_indices...> row_seq, std::integer_sequence<T, column_indices...> column_seq)
            -> std::array<std::array<std::array<Parameter, sizeof...(column_indices)>, sizeof...(row_indices)>, sizeof...(power_indices)>
        {
            return std::array<std::array<std::array<Parameter, sizeof...(column_indices)>, sizeof...(row_indices)>, sizeof...(power_indices)>
            {{
                _get_bkgpol_matrix(power_indices, row_seq, column_seq)...
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

            assume_isospin(destringify<bool>(o.get("assume_isospin", "false"))),

            m(_resonance_masses(p, u, std::make_index_sequence<EEToCCBarFull::nresonances>())),
            q(_resonance_sizes(p, u, std::make_index_sequence<EEToCCBarFull::nresonances>())),
            g0(_g0_matrix(p, u, std::make_index_sequence<EEToCCBarFull::nresonances>(), std::make_index_sequence<EEToCCBarFull::nchannels>())),
            bkgpol(_c_array(p, u, std::make_index_sequence<EEToCCBarFull::order + 1>(), std::make_index_sequence<EEToCCBarFull::nchannels>(), std::make_index_sequence<EEToCCBarFull::nchannels>())),

            Rconstant(p["ee->ccbar::Rconstant"], u),
            exclusive_norm(p["ee->ccbar::exclusive_norm"], u)
        {
            std::array<std::shared_ptr<KMatrix<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>::Resonance>, EEToCCBarFull::nresonances> resonance_array;
            for (unsigned i = 0; i<EEToCCBarFull::nresonances; i++)
            {
                resonance_array[i] = std::make_shared<CharmoniumResonanceFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(resonance_names[i], m[i], q[i]);
            }

            std::array<std::shared_ptr<KMatrix<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>::Channel>, EEToCCBarFull::nchannels> channel_array;
            for (unsigned i = 0; i<EEToCCBarFull::nchannels; i++)
            {
                switch (FullChannels(i))
                {
                    case ee:
                    case effJpsi:
                    case eff2S:
                    case eff3770:
                    case eff4040:
                    case eff4160:
                    case eff4415:
                        channel_array[i] = std::make_shared<EffChannelFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(channel_names[i], m_e, m_e, _get_g0_column(_filter_channel_index(FullChannels(i)), std::make_index_sequence<EEToCCBarFull::nresonances>()));
                        break;
                    case D0Dbar0:
                        channel_array[i] = std::make_shared<PWavePPChannelFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(channel_names[i], m_D0, m_D0, _get_g0_column(_filter_channel_index(FullChannels(i)), std::make_index_sequence<EEToCCBarFull::nresonances>()));
                        break;
                    case DpDm:
                        channel_array[i] = std::make_shared<PWavePPChannelFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(channel_names[i], m_D, m_D, _get_g0_column(_filter_channel_index(FullChannels(i)), std::make_index_sequence<EEToCCBarFull::nresonances>()));
                        break;
                    case D0Dbarst0:
                    case Dst0Dbar0:
                        channel_array[i] = std::make_shared<PWaveVPChannelFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(channel_names[i], m_D0, m_Dst0, _get_g0_column(_filter_channel_index(FullChannels(i)), std::make_index_sequence<EEToCCBarFull::nresonances>()));
                        break;
                    case DpDstm:
                    case DstpDm:
                        channel_array[i] = std::make_shared<PWaveVPChannelFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(channel_names[i], m_D, m_Dst, _get_g0_column(_filter_channel_index(FullChannels(i)), std::make_index_sequence<EEToCCBarFull::nresonances>()));
                        break;
                    case DspDsm:
                        channel_array[i] = std::make_shared<PWavePPChannelFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(channel_names[i], m_Ds, m_Ds, _get_g0_column(_filter_channel_index(FullChannels(i)), std::make_index_sequence<EEToCCBarFull::nresonances>()));
                        break;
                    case Dst0Dbarst0P0:
                    case Dst0Dbarst0P2:
                        channel_array[i] = std::make_shared<PWaveVVChannelFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(channel_names[i], m_Dst0, m_Dst0, _get_g0_column(_filter_channel_index(FullChannels(i)), std::make_index_sequence<EEToCCBarFull::nresonances>()));
                        break;
                    case Dst0Dbarst0F2:
                        channel_array[i] = std::make_shared<FWaveVVChannelFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(channel_names[i], m_Dst0, m_Dst0, _get_g0_column(_filter_channel_index(FullChannels(i)), std::make_index_sequence<EEToCCBarFull::nresonances>()));
                        break;
                    case DstpDstmP0:
                    case DstpDstmP2:
                        channel_array[i] = std::make_shared<PWaveVVChannelFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(channel_names[i], m_Dst, m_Dst, _get_g0_column(_filter_channel_index(FullChannels(i)), std::make_index_sequence<EEToCCBarFull::nresonances>()));
                        break;
                    case DstpDstmF2:
                        channel_array[i] = std::make_shared<FWaveVVChannelFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(channel_names[i], m_Dst, m_Dst, _get_g0_column(_filter_channel_index(FullChannels(i)), std::make_index_sequence<EEToCCBarFull::nresonances>()));
                        break;
                    case DspDsstm:
                    case DsstpDsm:
                        channel_array[i] = std::make_shared<PWaveVPChannelFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(channel_names[i], m_Ds, m_Dsst, _get_g0_column(_filter_channel_index(FullChannels(i)), std::make_index_sequence<EEToCCBarFull::nresonances>()));
                        break;
                    case DsstpDsstmP0:
                    case DsstpDsstmP2:
                        channel_array[i] = std::make_shared<PWaveVVChannelFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(channel_names[i], m_Dsst, m_Dsst, _get_g0_column(_filter_channel_index(FullChannels(i)), std::make_index_sequence<EEToCCBarFull::nresonances>()));
                        break;
                    case DsstpDsstmF2:
                        channel_array[i] = std::make_shared<FWaveVVChannelFull<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>>(channel_names[i], m_Dsst, m_Dsst, _get_g0_column(_filter_channel_index(FullChannels(i)), std::make_index_sequence<EEToCCBarFull::nresonances>()));
                        break;

                    default:
                        throw InternalError("The number of requested channels (array) is larger than the number of known channel names.");
                }
            }

            K = std::shared_ptr<KMatrix<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>> (
                new KMatrix<EEToCCBarFull::nchannels, EEToCCBarFull::nresonances, EEToCCBarFull::order>(
                    channel_array,
                    resonance_array,
                    _get_bkgpol(std::make_index_sequence<EEToCCBarFull::order + 1>(), std::make_index_sequence<EEToCCBarFull::nchannels>(), std::make_index_sequence<EEToCCBarFull::nchannels>()),
                    "e^+e^-->ccbar"
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
            const double speedoflight = 299792458.; // Exact value
            const double GeVtonb = 10 * pow( 1.e18 * hbar * speedoflight, 2.0);

            return GeVtonb * 4.0 * M_PI * alpha_em * alpha_em / (3.0 * E*E);
        }

        // amplitude of ee -> channel
        complex<double> amplitude_eetochannel(const IntermediateResult * intermediate_result, const FullChannels & channel)
        {
            // Conversion factor between GeV^2 and nb
            const double speedoflight = 299792458.; //Exact value
            const double GeVtonb = 10 * pow( 1.e18 * hbar * speedoflight, 2.0);

            // Channel properties
            const double Nf = 2 * K->_channels[FullChannels(channel)]->_l_orbital + 1;
            const double rhof = real(K->_channels[FullChannels(channel)]->rho(intermediate_result->s));

            // Get T-matrix[ee, channel]
            const complex<double> T1f = intermediate_result->tmatrix_row_0[FullChannels(channel)];

            return sqrt(GeVtonb * 16. * M_PI / intermediate_result->s * Nf * rhof) * T1f;
        }

        double sigma_eetochannel(const IntermediateResult * intermediate_result, const FullChannels & channel)
        {
            return norm(amplitude_eetochannel(intermediate_result, channel));
        }

        // Widths
        double res_partial_width(const FullResonances & resonance, const FullChannels & channel)
        {
            return K->partial_width(resonance, channel);
        }

        double res_total_width(const FullResonances & resonance)
        {
            return K->width(resonance);
        }

        // Ratios
        double R(const IntermediateResult * intermediate_result)
        {
            double total_xsec = 0.0;

            for (unsigned i = 1; i < EEToCCBarFull::nchannels; i++)
            {
                total_xsec += sigma_eetochannel(intermediate_result, FullChannels(i));
            }

            return total_xsec / sigma_eetomumu(intermediate_result->E) + Rconstant; //Add constant term
        }
    };

    EEToCCBarFull::EEToCCBarFull(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<EEToCCBarFull>(new Implementation<EEToCCBarFull>(parameters, options, *this))
    {
    }

    EEToCCBarFull::~EEToCCBarFull()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<EEToCCBarFull>::options
    {
        {"assume_isospin", { "true", "false" }, "false"},
    };

    const EEToCCBarFull::IntermediateResult *
    EEToCCBarFull::prepare(const double & E) const
    {
        return _imp->prepare(E);
    }

    using FullResonances = Implementation<eos::EEToCCBarFull>::FullResonances;
    using FullChannels = Implementation<eos::EEToCCBarFull>::FullChannels;

    double
    EEToCCBarFull::Jpsi_ee_width() const
    {
        return _imp->res_partial_width(FullResonances::Jpsi, FullChannels::ee);
    }

    double
    EEToCCBarFull::Jpsi_eff_width() const
    {
        return _imp->res_partial_width(FullResonances::Jpsi, FullChannels::eff2S);
    }

    double
    EEToCCBarFull::Jpsi_total_width() const
    {
        return _imp->res_total_width(FullResonances::Jpsi);
    }

    double
    EEToCCBarFull::psi2S_ee_width() const
    {
        return _imp->res_partial_width(FullResonances::psi2S, FullChannels::ee);
    }

    double
    EEToCCBarFull::psi2S_eff_width() const
    {
        return _imp->res_partial_width(FullResonances::psi2S, FullChannels::eff2S);
    }

    double
    EEToCCBarFull::psi2S_total_width() const
    {
        return _imp->res_total_width(FullResonances::psi2S);
    }

    double
    EEToCCBarFull::psi3770_total_width() const
    {
        return _imp->res_total_width(FullResonances::psi3770);
    }

    double
    EEToCCBarFull::psi4040_total_width() const
    {
        return _imp->res_total_width(FullResonances::psi4040);
    }

    double
    EEToCCBarFull::psi4160_total_width() const
    {
        return _imp->res_total_width(FullResonances::psi4160);
    }

    double
    EEToCCBarFull::psi4415_total_width() const
    {
        return _imp->res_total_width(FullResonances::psi4415);
    }

    double
    EEToCCBarFull::psi3770_D0Dbar0_width() const
    {
        return _imp->res_partial_width(FullResonances::psi3770, FullChannels::D0Dbar0);
    }

    double
    EEToCCBarFull::psi3770_DpDm_width() const
    {
        return _imp->res_partial_width(FullResonances::psi3770, FullChannels::DpDm);
    }

    double
    EEToCCBarFull::psi3770_eff_width() const
    {
        return _imp->res_partial_width(FullResonances::psi3770, FullChannels::eff3770);
    }

    double
    EEToCCBarFull::psi4040_DD_width() const
    {
        return _imp->res_partial_width(FullResonances::psi4040, FullChannels::D0Dbar0)
             + _imp->res_partial_width(FullResonances::psi4040, FullChannels::DpDm);
    }

    double
    EEToCCBarFull::psi4040_DDst_width() const
    {
        return _imp->res_partial_width(FullResonances::psi4040, FullChannels::D0Dbarst0)
             + _imp->res_partial_width(FullResonances::psi4040, FullChannels::Dst0Dbar0)
             + _imp->res_partial_width(FullResonances::psi4040, FullChannels::DpDstm)
             + _imp->res_partial_width(FullResonances::psi4040, FullChannels::DstpDm);
    }

    double
    EEToCCBarFull::psi4040_DstDst_width() const
    {
        return _imp->res_partial_width(FullResonances::psi4040, FullChannels::Dst0Dbarst0P0)
             + _imp->res_partial_width(FullResonances::psi4040, FullChannels::Dst0Dbarst0P2)
             + _imp->res_partial_width(FullResonances::psi4040, FullChannels::Dst0Dbarst0F2)
             + _imp->res_partial_width(FullResonances::psi4040, FullChannels::DstpDstmP0)
             + _imp->res_partial_width(FullResonances::psi4040, FullChannels::DstpDstmP2)
             + _imp->res_partial_width(FullResonances::psi4040, FullChannels::DstpDstmF2);
    }

    double
    EEToCCBarFull::psi4160_DD_width() const
    {
        return _imp->res_partial_width(FullResonances::psi4160, FullChannels::D0Dbar0)
             + _imp->res_partial_width(FullResonances::psi4160, FullChannels::DpDm);
    }

    double
    EEToCCBarFull::psi4160_DDst_width() const
    {
        return _imp->res_partial_width(FullResonances::psi4160, FullChannels::D0Dbarst0)
             + _imp->res_partial_width(FullResonances::psi4160, FullChannels::Dst0Dbar0)
             + _imp->res_partial_width(FullResonances::psi4160, FullChannels::DpDstm)
             + _imp->res_partial_width(FullResonances::psi4160, FullChannels::DstpDm);
    }

    double
    EEToCCBarFull::psi4160_DstDst_width() const
    {
        return _imp->res_partial_width(FullResonances::psi4160, FullChannels::Dst0Dbarst0P0)
             + _imp->res_partial_width(FullResonances::psi4160, FullChannels::Dst0Dbarst0P2)
             + _imp->res_partial_width(FullResonances::psi4160, FullChannels::Dst0Dbarst0F2)
             + _imp->res_partial_width(FullResonances::psi4160, FullChannels::DstpDstmP0)
             + _imp->res_partial_width(FullResonances::psi4160, FullChannels::DstpDstmP2)
             + _imp->res_partial_width(FullResonances::psi4160, FullChannels::DstpDstmF2);
    }

    double
    EEToCCBarFull::psi4415_DD_width() const
    {
        return _imp->res_partial_width(FullResonances::psi4415, FullChannels::D0Dbar0)
             + _imp->res_partial_width(FullResonances::psi4415, FullChannels::DpDm);
    }

    double
    EEToCCBarFull::psi4415_DDst_width() const
    {
        return _imp->res_partial_width(FullResonances::psi4415, FullChannels::D0Dbarst0)
             + _imp->res_partial_width(FullResonances::psi4415, FullChannels::Dst0Dbar0)
             + _imp->res_partial_width(FullResonances::psi4415, FullChannels::DpDstm)
             + _imp->res_partial_width(FullResonances::psi4415, FullChannels::DstpDm);
    }

    double
    EEToCCBarFull::psi4415_DstDst_width() const
    {
        return _imp->res_partial_width(FullResonances::psi4415, FullChannels::Dst0Dbarst0P0)
             + _imp->res_partial_width(FullResonances::psi4415, FullChannels::Dst0Dbarst0P2)
             + _imp->res_partial_width(FullResonances::psi4415, FullChannels::Dst0Dbarst0F2)
             + _imp->res_partial_width(FullResonances::psi4415, FullChannels::DstpDstmP0)
             + _imp->res_partial_width(FullResonances::psi4415, FullChannels::DstpDstmP2)
             + _imp->res_partial_width(FullResonances::psi4415, FullChannels::DstpDstmF2);
    }

    double
    EEToCCBarFull::sigma_eetoee(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, FullChannels::ee);
    }

    double
    EEToCCBarFull::sigma_eetoeff(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, FullChannels::effJpsi)
            + _imp->sigma_eetochannel(ir, FullChannels::eff2S)
            + _imp->sigma_eetochannel(ir, FullChannels::eff3770)
            + _imp->sigma_eetochannel(ir, FullChannels::eff4040)
            + _imp->sigma_eetochannel(ir, FullChannels::eff4160)
            + _imp->sigma_eetochannel(ir, FullChannels::eff4415)
            );
    }

    double
    EEToCCBarFull::sigma_eetoD0Dbar0(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, FullChannels::D0Dbar0);
    }

    double
    EEToCCBarFull::sigma_eetoDpDm(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, FullChannels::DpDm);
    }

    double
    EEToCCBarFull::sigma_eetoD0Dbarst0(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, FullChannels::D0Dbarst0)
            + _imp->sigma_eetochannel(ir, FullChannels::Dst0Dbar0)
            );
    }

    double
    EEToCCBarFull::sigma_eetoDpDstm(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, FullChannels::DpDstm)
            + _imp->sigma_eetochannel(ir, FullChannels::DstpDm)
            );
    }

    double
    EEToCCBarFull::sigma_eetoDspDsm(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, FullChannels::DspDsm);
    }

    double
    EEToCCBarFull::sigma_eetoDst0Dbarst0(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, FullChannels::Dst0Dbarst0P0)
            + _imp->sigma_eetochannel(ir, FullChannels::Dst0Dbarst0P2)
            + _imp->sigma_eetochannel(ir, FullChannels::Dst0Dbarst0F2)
            );
    }

    double
    EEToCCBarFull::sigma_eetoDstpDstm(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, FullChannels::DstpDstmP0)
            + _imp->sigma_eetochannel(ir, FullChannels::DstpDstmP2)
            + _imp->sigma_eetochannel(ir, FullChannels::DstpDstmF2)
            );
    }

    double
    EEToCCBarFull::sigma_eetoDstpLDstmL(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (1.0 / 75.0) * norm(
                                  5.0 * _imp->amplitude_eetochannel(ir, FullChannels::DstpDstmP0)
                     - pow(30.0, 0.5) * _imp->amplitude_eetochannel(ir, FullChannels::DstpDstmF2)
                + 2. * pow( 5.0, 0.5) * _imp->amplitude_eetochannel(ir, FullChannels::DstpDstmP2)
            );
    }

    double
    EEToCCBarFull::sigma_eetoDstpTDstmL(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (1.0 / 5.0) * norm(
                  pow(2.0, 0.5) * _imp->amplitude_eetochannel(ir, FullChannels::DstpDstmF2)
                + pow(3.0, 0.5) * _imp->amplitude_eetochannel(ir, FullChannels::DstpDstmP2)
            );
    }

    double
    EEToCCBarFull::sigma_eetoDstpTDstmT(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (1.0 / 150.0) * norm(
                                 10.0 * _imp->amplitude_eetochannel(ir, FullChannels::DstpDstmP0)
                     + pow(30.0, 0.5) * _imp->amplitude_eetochannel(ir, FullChannels::DstpDstmF2)
                - 2. * pow( 5.0, 0.5) * _imp->amplitude_eetochannel(ir, FullChannels::DstpDstmP2)
            );
    }

    double
    EEToCCBarFull::sigma_eetoDspDsstm(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, FullChannels::DspDsstm)
            + _imp->sigma_eetochannel(ir, FullChannels::DsstpDsm));
    }

    double
    EEToCCBarFull::sigma_eetoDsstpDsstm(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, FullChannels::DsstpDsstmP0)
            + _imp->sigma_eetochannel(ir, FullChannels::DsstpDsstmP2)
            + _imp->sigma_eetochannel(ir, FullChannels::DsstpDsstmF2)
            );
    }

    double
    EEToCCBarFull::R(const EEToCCBarFull::IntermediateResult * ir) const
    {
        return _imp->R(ir);
    }

    const std::set<ReferenceName>
    EEToCCBarFull::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    EEToCCBarFull::begin_options()
    {
        return Implementation<EEToCCBarFull>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    EEToCCBarFull::end_options()
    {
        return Implementation<EEToCCBarFull>::options.cend();
    }
}
