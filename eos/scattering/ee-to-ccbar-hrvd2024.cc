/*
 * Copyright (c) 2024 Méril Reboud
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

#include <eos/scattering/ee-to-ccbar-hrvd2024.hh>
#include <eos/scattering/ee-to-ccbar.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/polylog.hh>
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
    struct Implementation<EEToCCBarHRvD>
    {
        UsedParameter hbar;
        UsedParameter alpha_em;
        UsedParameter m_e;
        UsedParameter m_eff;
        UsedParameter m_D0;
        UsedParameter m_Dp;
        UsedParameter m_Dst0;
        UsedParameter m_Dstp;
        UsedParameter m_Dsp;
        UsedParameter m_Dsstp;

        BooleanOption assume_isospin_opt;
        BooleanOption assume_uspin_opt;
        bool assume_isospin;
        bool assume_uspin;

        static const inline std::array<std::string, EEToCCBarHRvD::nresonances> resonance_names =
        {
            "psi(2S)", "psi(3770)", "psi(4040)", "psi(4160)"
        };

        enum Resonances
        {
            psi2S = 0, psi3770, psi4040, psi4160
        };

        static const inline std::array<std::string, EEToCCBarHRvD::nchannels> channel_names =
        {
            "e^+e^-", "eff(2S)", "D^0Dbar^0", "D^+D^-", "eff(3770)",
            "D^0Dbar^*0", "D^*0Dbar^0", "D^+D^*-", "D^*+D^-", "D_s^+D_s^-",
            "D^*0Dbar^*0P0", "D^*0Dbar^*0P2", "D^*0Dbar^*0F2",
            "D^*+D^*-P0", "D^*+D^*-P2", "D^*+D^*-F2", "eff(4040)",
            "D_s^+D_s^*-", "D_s^*+D_s^-", "eff(4160)",
            // "D_s^*+D_s^*-P0", "D_s^*+D_s^*-P2", "D_s^*+D_s^*-F2"
        };

        enum Channels
        {
            ee = 0, eff2S, D0Dbar0, DpDm, eff3770,
            D0Dbarst0, Dst0Dbar0, DpDstm, DstpDm, DspDsm,
            Dst0Dbarst0P0, Dst0Dbarst0P2, Dst0Dbarst0F2,
            DstpDstmP0, DstpDstmP2, DstpDstmF2, eff4040,
            DspDsstm, DsstpDsm, eff4160,
            // DsstpDsstmP0, DsstpDsstmP2, DsstpDsstmF2
        };

        // Resonance masses
        std::array<UsedParameter, EEToCCBarHRvD::nresonances> m;

        // Channel-Resonance couplings
        std::array<std::array<UsedParameter, EEToCCBarHRvD::nchannels>, EEToCCBarHRvD::nresonances> g0;

        // Channels barrier factors scales
        std::array<UsedParameter, EEToCCBarHRvD::nchannels> q;

        // Non resonant contributions to the K matrix
        std::array<std::array<UsedParameter, EEToCCBarHRvD::nchannels>, EEToCCBarHRvD::nchannels> bkgcst;

        // R_uds
        UsedParameter Rconstant;

        std::shared_ptr<KMatrix<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>> K;

        // QED corrections
        BooleanOption opt_qed;
        std::function<double(const double &)> eta;

        using IntermediateResult = EEToCCBarHRvD::IntermediateResult;
        IntermediateResult _intermediate_result;

        template <typename T, T... indices>
        auto _resonance_masses(const Parameters & p, ParameterUser & u, std::integer_sequence<T, indices...>)
            -> std::array<UsedParameter, sizeof...(indices)>
        {
            static_assert(sizeof...(indices) <= resonance_names.size(), "The number of requested resonances is larger than the number of resonance names.");

            return std::array<UsedParameter, sizeof...(indices)>
            {{
                UsedParameter(p["mass::" + resonance_names[indices]], u)...
            }};
        }

        template <typename T, T... indices>
        auto _channel_effective_momentum(const Parameters & p, ParameterUser & u, std::integer_sequence<T, indices...>)
            -> std::array<UsedParameter, sizeof...(indices)>
        {
            static_assert(sizeof...(indices) <= channel_names.size(), "The number of requested channels is larger than the number of channel names.");

            return std::array<UsedParameter, sizeof...(indices)>
            {{
                UsedParameter(p["ee->ccbar::q_0(" + channel_names[indices] + ")"], u)...
            }};
        }

        long unsigned _filter_channel_index(Channels channel)
        {
            switch (channel)
            {
                case DpDm:
                    return D0Dbar0;

                case Dst0Dbar0:
                    return D0Dbarst0;

                case DpDstm:
                    if (assume_isospin) return D0Dbarst0;
                    break;

                case DstpDm:
                    if (assume_isospin) return D0Dbarst0;
                    return DpDstm;

                case DspDsm:
                    if (assume_uspin) return D0Dbar0;
                    break;

                case DstpDstmP0:
                    if (assume_isospin) return Dst0Dbarst0P0;
                    break;

                case DstpDstmP2:
                    if (assume_isospin) return Dst0Dbarst0P2;
                    break;

                case DstpDstmF2:
                    if (assume_isospin) return Dst0Dbarst0F2;
                    break;

                case DspDsstm:
                    if (assume_uspin) return D0Dbarst0;
                    break;

                case DsstpDsm:
                    if (assume_uspin) return D0Dbarst0;
                    return DspDsstm;

                // case DsstpDsstmP0:
                //     if (assume_uspin) return Dst0Dbarst0P0;
                //     break;

                // case DsstpDsstmP2:
                //     if (assume_uspin) return Dst0Dbarst0P2;
                //     break;

                // case DsstpDsstmF2:
                //     if (assume_uspin) return Dst0Dbarst0F2;
                //     break;

                default:
                    return channel;
            }

            return channel;
        }

        template <typename T, T... column_indices>
        auto _g0_row(const Parameters & p, ParameterUser & u, T row_index, std::integer_sequence<T, column_indices...>)
            -> std::array<UsedParameter, sizeof...(column_indices)>
        {
            static_assert(sizeof...(column_indices) <= channel_names.size(), "The number of requested channels is larger than the number of channel names.");

            return std::array<UsedParameter, sizeof...(column_indices)>
            {{
                UsedParameter(p["ee->ccbar::g0(" + resonance_names[row_index] + "," + channel_names[_filter_channel_index(Channels(column_indices))] + ")"], u)...
            }};
        }

        template <typename T, T... row_indices, T... column_indices>
        auto _g0_matrix(const Parameters & p, ParameterUser & u, std::integer_sequence<T, row_indices...>, std::integer_sequence<T, column_indices...> column_seq)
            -> std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>
        {
            static_assert(sizeof...(row_indices) <= resonance_names.size(), "The number of requested resonances is larger than the number of resonance names.");

            return std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>
            {{
                _g0_row(p, u, row_indices, column_seq)...
            }};
        }

        template <typename T>
        std::string _channel_name_tuple(const T & a, const T & b)
        {
            if (b > a)
                return "(" + channel_names[_filter_channel_index(Channels(a))] + "," + channel_names[_filter_channel_index(Channels(b))] + ")";

            return "(" + channel_names[_filter_channel_index(Channels(b))] + "," + channel_names[_filter_channel_index(Channels(a))] + ")";
        }

        template <typename T, T... column_indices>
        auto _c_row(const Parameters & p, ParameterUser & u, T row_index, std::integer_sequence<T, column_indices...>)
            -> std::array<UsedParameter, sizeof...(column_indices)>
        {
            static_assert(sizeof...(column_indices) <= channel_names.size(), "The number of requested channels (column) is larger than the number of channel names.");

            return std::array<UsedParameter, sizeof...(column_indices)>
            {{
                UsedParameter(p["ee->ccbar::c" + _channel_name_tuple(row_index, column_indices)], u)...
            }};
        }

        template <typename T, T... row_indices, T... column_indices>
        auto _c_matrix(const Parameters & p, ParameterUser & u, std::integer_sequence<T, row_indices...>, std::integer_sequence<T, column_indices...> column_seq)
            -> std::array<std::array<UsedParameter, sizeof...(column_indices)>, sizeof...(row_indices)>
        {
            static_assert(sizeof...(row_indices) <= channel_names.size(), "The number of requested channels (row) is larger than the number of channel names.");

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
        auto _get_bkgcst_matrix(std::integer_sequence<T, row_indices...>, std::integer_sequence<T, column_indices...> column_seq)
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
            m_eff(p["ee->ccbar::effective_mass"], u),
            m_D0(p["mass::D^0"], u),
            m_Dp(p["mass::D^+"], u),
            m_Dst0(p["mass::D_u^*"], u),
            m_Dstp(p["mass::D_d^*"], u),
            m_Dsp(p["mass::D_s"], u),
            m_Dsstp(p["mass::D_s^*"], u),
            assume_isospin_opt(o ,options, "assume_isospin"),
            assume_uspin_opt(o ,options, "assume_uspin"),
            assume_isospin(assume_isospin_opt.value()),
            assume_uspin(assume_uspin_opt.value()),
            m(_resonance_masses(p, u, std::make_index_sequence<EEToCCBarHRvD::nresonances>())),
            g0(_g0_matrix(p, u, std::make_index_sequence<EEToCCBarHRvD::nresonances>(), std::make_index_sequence<EEToCCBarHRvD::nchannels>())),
            q(_channel_effective_momentum(p, u, std::make_index_sequence<EEToCCBarHRvD::nchannels>())),
            bkgcst(_c_matrix(p, u, std::make_index_sequence<EEToCCBarHRvD::nchannels>(), std::make_index_sequence<EEToCCBarHRvD::nchannels>())),
            Rconstant(p["ee->ccbar::Rconstant"], u),
            opt_qed(o, options, "QED"),
            eta([this](const double &) { return 0.0; })
        {
            std::array<std::shared_ptr<KMatrix<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>::Resonance>, EEToCCBarHRvD::nresonances> resonance_array;
            for (unsigned i = 0 ; i < EEToCCBarHRvD::nresonances ; i++)
            {
                resonance_array[i] = std::make_shared<CharmoniumResonance<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(resonance_names[i], m[i]);
            }

            std::array<std::shared_ptr<KMatrix<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>::Channel>, EEToCCBarHRvD::nchannels> channel_array;
            for (unsigned i = 0 ; i < EEToCCBarHRvD::nchannels ; i++)
            {
                switch (Channels(i))
                {
                    case ee:
                        channel_array[i] = std::make_shared<EEChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_e, m_e, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                        break;
                    case eff2S:
                    case eff3770:
                    case eff4040:
                    case eff4160:
                        channel_array[i] = std::make_shared<EffChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_eff, m_eff, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                        break;
                    case D0Dbar0:
                        channel_array[i] = std::make_shared<PWavePPChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_D0, m_D0, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                        break;
                    case DpDm:
                        channel_array[i] = std::make_shared<PWavePPChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_Dp, m_Dp, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                        break;
                    case D0Dbarst0:
                    case Dst0Dbar0:
                        channel_array[i] = std::make_shared<PWaveVPChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_Dst0, m_D0, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                        break;
                    case DpDstm:
                    case DstpDm:
                        channel_array[i] = std::make_shared<PWaveVPChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_Dstp, m_Dp, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                        break;
                    case DspDsm:
                        channel_array[i] = std::make_shared<PWavePPChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_Dsp, m_Dsp, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                        break;
                    case Dst0Dbarst0P0:
                    case Dst0Dbarst0P2:
                        channel_array[i] = std::make_shared<PWaveVVChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_Dst0, m_Dst0, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                        break;
                    case Dst0Dbarst0F2:
                        channel_array[i] = std::make_shared<FWaveVVChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_Dst0, m_Dst0, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                        break;
                    case DstpDstmP0:
                    case DstpDstmP2:
                        channel_array[i] = std::make_shared<PWaveVVChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_Dstp, m_Dstp, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                        break;
                    case DstpDstmF2:
                        channel_array[i] = std::make_shared<FWaveVVChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_Dstp, m_Dstp, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                        break;
                    case DspDsstm:
                    case DsstpDsm:
                        channel_array[i] = std::make_shared<PWaveVPChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_Dsstp, m_Dsp, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                        break;
                    // case DsstpDsstmP0:
                    // case DsstpDsstmP2:
                    //     channel_array[i] = std::make_shared<PWaveVVChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_Dsstp, m_Dsstp, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                    //     break;
                    // case DsstpDsstmF2:
                    //     channel_array[i] = std::make_shared<FWaveVVChannel<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>>(channel_names[i], m_Dsstp, m_Dsstp, q[i], _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBarHRvD::nresonances>()));
                    //     break;

                    default:
                        throw InternalError("The number of requested channels (array) is larger than the number of known channel names.");
                }
            }

            K = std::shared_ptr<KMatrix<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>> (
                new KMatrix<EEToCCBarHRvD::nchannels, EEToCCBarHRvD::nresonances>(
                    channel_array,
                    resonance_array,
                    _get_bkgcst_matrix(std::make_index_sequence<EEToCCBarHRvD::nchannels>(), std::make_index_sequence<EEToCCBarHRvD::nchannels>()),
                    "e^+e^-->ccbar"
                    )
                );

            if (opt_qed.value())
            {
                eta = [this](const double & s) { return this->_eta_lo(s); };
            }
        }

        double _eta_lo(const double & s) const
        {
            // hep-ph/0107154 eq. (45)
            const double beta = sqrt(1.0 - 4.0 * m_Dp * m_Dp / s), beta2 = beta * beta,
                beta3 = beta2 * beta;
            const double X = (1.0 - beta) / (1.0 + beta);

            return (1.0 + beta2) / beta * real(
                4.0 * dilog(X) + 2.0 * dilog(-1.0 * X) + 3.0 * log(X) * log(2.0 / (1.0 + beta))
                + 2.0 * log(beta) * log(X))
                - 3.0 * log(4.0 / (1 - beta2)) - 4.0 * log(beta)
                - 1.0 / beta3 * (5.0 / 4.0 * power_of<2>(1.0 + beta2) - 2.0) * log(X)
                + 3.0 / 2.0 * (1.0 + beta2) / beta2
            ;
        }

        const IntermediateResult * prepare(const complex<double> & E)
        {
            // Amplitude on the first RS
            _intermediate_result.tmatrix_row_0 = K->tmatrix_row(0, E * E);
            // Amplitude on the second RS
            _intermediate_result.tmatrix2_row_0 = K->tmatrix_row(0, E * E, true);

            _intermediate_result.E = E;
            _intermediate_result.s = E * E;

            return &_intermediate_result;
        }


        double rho(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            return std::real(K->_channels[channel]->rho(intermediate_result->s));
        }

        complex<double> chew_mandelstam(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            return K->_channels[channel]->chew_mandelstam(intermediate_result->s);
        }

        complex<double> chew_mandelstam_II(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            complex<double> s = intermediate_result->s;
            auto get_channel = K->_channels[channel];

            const unsigned li = get_channel->_l_orbital;
            const double q0 = get_channel->_q0.evaluate();
            const complex<double> mi1_2 = power_of<2>(get_channel->_m1());
            const complex<double> mi2_2 = power_of<2>(get_channel->_m2());

            // Momentum of particles in their center-of-momentum frame
            const complex<double> q = 0.5 * sqrt(eos::lambda(s, mi1_2, mi2_2)) / sqrt(s);

            // Blatt-Weisskopf factors, cf eq. (50.26)
            const complex<double> Fi = kmatrix_utils::blatt_weisskopf_factor(li, q / q0);

            return get_channel->chew_mandelstam(s) +
                complex<double>(0.0, 2.0) * get_channel->rho(s) * power_of<2>(pow(q / q0, li) * Fi);
        }


        inline double sigma_eetomumu_leading_order(const double & E)
        {
            // Conversion factor between GeV^2 and nb
            static const double speedoflight = 299792458.; // Exact value
            const double GeVtonb = 10.0 * power_of<2>(1.0e18 * hbar * speedoflight);

            return GeVtonb * 4.0 * M_PI * alpha_em * alpha_em / (3.0 * E * E);
        }

        // amplitude of ee -> channel
        complex<double> T_eetochannel(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            // Get T-matrix[ee, channel]
            const complex<double> T1f = intermediate_result->tmatrix_row_0[channel];

            return T1f;
        }

        complex<double> T_II_eetochannel(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            // Get T-matrix[ee, channel] on the second RS
            const complex<double> T1f = intermediate_result->tmatrix2_row_0[channel];

            return T1f;
        }

        double sigma_eetochannel(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            // Conversion factor between GeV^2 and nb
            static const double speedoflight = 299792458.; // Exact value
            const double GeVtonb = 10.0 * power_of<2>(1.0e18 * hbar * speedoflight);

            // Channel properties
            const double Nf = 2.0 * K->_channels[channel]->_l_orbital + 1.0;
            const double rhof = std::abs(K->_channels[channel]->rho(intermediate_result->s));

            return GeVtonb / std::abs(intermediate_result->s) * Nf * rhof * norm(T_eetochannel(intermediate_result, channel));
        }

        // K matrix widths, they are not expected to match the experimental ones
        double res_partial_width(const Resonances & resonance, const Channels & channel)
        {
            return K->partial_width(resonance, channel);
        }

        double res_total_width(const Resonances & resonance)
        {
            return K->width(resonance);
        }

        // Ratios
        double R(const IntermediateResult * intermediate_result)
        {
            double total_sigma = 0.0;

            for (unsigned i = 1 ; i < EEToCCBarHRvD::nchannels ; i++)
            {
                total_sigma += sigma_eetochannel(intermediate_result, Channels(i));
            }

            return total_sigma / sigma_eetomumu_leading_order(abs(intermediate_result->E)) + Rconstant; // Add constant term
        }

        // Spectral function
        double spectral_function(const double E, const Resonances & resonance)
        {
            return K->spectral_function(resonance, E * E);
        }
    };

    EEToCCBarHRvD::EEToCCBarHRvD(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<EEToCCBarHRvD>(new Implementation<EEToCCBarHRvD>(parameters, options, *this))
    {
    }

    EEToCCBarHRvD::~EEToCCBarHRvD()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<EEToCCBarHRvD>::options
    {
        {"assume_isospin", { "true", "false" }, "true" },
        {"assume_uspin",   { "true", "false" }, "true" },
        {"QED",            { "true", "false" }, "false"},
    };

    const EEToCCBarHRvD::IntermediateResult *
    EEToCCBarHRvD::prepare(const double & E) const
    {
        return _imp->prepare(E);
    }

    const EEToCCBarHRvD::IntermediateResult *
    EEToCCBarHRvD::prepare_complex(const double & re_E, const double & im_E) const
    {
        return _imp->prepare(complex<double>(re_E, im_E));
    }

    using Resonances = Implementation<eos::EEToCCBarHRvD>::Resonances;
    using Channels = Implementation<eos::EEToCCBarHRvD>::Channels;

    double
    EEToCCBarHRvD::psi2S_ee_width() const
    {
        return _imp->res_partial_width(Resonances::psi2S, Channels::ee);
    }

    double
    EEToCCBarHRvD::psi2S_eff_width() const
    {
        return _imp->res_partial_width(Resonances::psi2S, Channels::eff2S);
    }

    double
    EEToCCBarHRvD::psi2S_total_width() const
    {
        return _imp->res_total_width(Resonances::psi2S);
    }

    double
    EEToCCBarHRvD::psi3770_total_width() const
    {
        return _imp->res_total_width(Resonances::psi3770);
    }

    double
    EEToCCBarHRvD::psi3770_D0Dbar0_width() const
    {
        return _imp->res_partial_width(Resonances::psi3770, Channels::D0Dbar0);
    }

    double
    EEToCCBarHRvD::psi3770_DpDm_width() const
    {
        return _imp->res_partial_width(Resonances::psi3770, Channels::DpDm);
    }

    double
    EEToCCBarHRvD::psi3770_eff_width() const
    {
        return _imp->res_partial_width(Resonances::psi3770, Channels::eff3770);
    }

    double
    EEToCCBarHRvD::psi4040_D0Dbar0_width() const
    {
        return _imp->res_partial_width(Resonances::psi4040, Channels::D0Dbar0);
    }

    double
    EEToCCBarHRvD::psi4040_DpDm_width() const
    {
        return _imp->res_partial_width(Resonances::psi4040, Channels::DpDm);
    }

    double
    EEToCCBarHRvD::psi4040_D0Dbarst0_width() const
    {
        return _imp->res_partial_width(Resonances::psi4040, Channels::D0Dbarst0)
             + _imp->res_partial_width(Resonances::psi4040, Channels::Dst0Dbar0);
    }

    double
    EEToCCBarHRvD::psi4040_DpDstm_width() const
    {
        return _imp->res_partial_width(Resonances::psi4040, Channels::DpDstm)
             + _imp->res_partial_width(Resonances::psi4040, Channels::DstpDm);
    }

    double
    EEToCCBarHRvD::psi4040_DspDsm_width() const
    {
        return _imp->res_partial_width(Resonances::psi4040, Channels::DspDsm);
    }

    double
    EEToCCBarHRvD::psi4040_Dst0Dbarst0_width() const
    {
        return _imp->res_partial_width(Resonances::psi4040, Channels::Dst0Dbarst0P0)
             + _imp->res_partial_width(Resonances::psi4040, Channels::Dst0Dbarst0P2)
             + _imp->res_partial_width(Resonances::psi4040, Channels::Dst0Dbarst0F2);
    }

    double
    EEToCCBarHRvD::psi4040_DstpDstm_width() const
    {
        return _imp->res_partial_width(Resonances::psi4040, Channels::DstpDstmP0)
             + _imp->res_partial_width(Resonances::psi4040, Channels::DstpDstmP2)
             + _imp->res_partial_width(Resonances::psi4040, Channels::DstpDstmF2);
    }

    double
    EEToCCBarHRvD::psi4040_eff_width() const
    {
        return _imp->res_partial_width(Resonances::psi4040, Channels::eff4040);
    }

    double
    EEToCCBarHRvD::psi4040_total_width() const
    {
        return _imp->res_total_width(Resonances::psi4040);
    }

    double
    EEToCCBarHRvD::psi4160_D0Dbar0_width() const
    {
        return _imp->res_partial_width(Resonances::psi4160, Channels::D0Dbar0);
    }

    double
    EEToCCBarHRvD::psi4160_DpDm_width() const
    {
        return _imp->res_partial_width(Resonances::psi4160, Channels::DpDm);
    }

    double
    EEToCCBarHRvD::psi4160_D0Dbarst0_width() const
    {
        return _imp->res_partial_width(Resonances::psi4160, Channels::D0Dbarst0)
            +  _imp->res_partial_width(Resonances::psi4160, Channels::Dst0Dbar0);
    }

    double
    EEToCCBarHRvD::psi4160_DpDstm_width() const
    {
        return _imp->res_partial_width(Resonances::psi4160, Channels::DpDstm)
             + _imp->res_partial_width(Resonances::psi4160, Channels::DstpDm);
    }

    double
    EEToCCBarHRvD::psi4160_DspDsm_width() const
    {
        return _imp->res_partial_width(Resonances::psi4160, Channels::DspDsm);
    }

    double
    EEToCCBarHRvD::psi4160_Dst0Dbarst0_width() const
    {
        return _imp->res_partial_width(Resonances::psi4160, Channels::Dst0Dbarst0P0)
             + _imp->res_partial_width(Resonances::psi4160, Channels::Dst0Dbarst0P2)
             + _imp->res_partial_width(Resonances::psi4160, Channels::Dst0Dbarst0F2);
    }

    double
    EEToCCBarHRvD::psi4160_DstpDstm_width() const
    {
        return _imp->res_partial_width(Resonances::psi4160, Channels::DstpDstmP0)
             + _imp->res_partial_width(Resonances::psi4160, Channels::DstpDstmP2)
             + _imp->res_partial_width(Resonances::psi4160, Channels::DstpDstmF2);
    }

    double
    EEToCCBarHRvD::psi4160_DspDsstm_width() const
    {
        return _imp->res_partial_width(Resonances::psi4160, Channels::DspDsstm)
              +_imp->res_partial_width(Resonances::psi4160, Channels::DsstpDsm);
    }

    double
    EEToCCBarHRvD::psi4160_eff_width() const
    {
        return _imp->res_partial_width(Resonances::psi4160, Channels::eff4160);
    }

    double
    EEToCCBarHRvD::psi4160_total_width() const
    {
        return _imp->res_total_width(Resonances::psi4160);
    }



    double
    EEToCCBarHRvD::sigma_eetoee(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return _imp->sigma_eetochannel(ir, Channels::ee);
    }

    double
    EEToCCBarHRvD::sigma_eetoeff(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return (
            + _imp->sigma_eetochannel(ir, Channels::eff2S)
            + _imp->sigma_eetochannel(ir, Channels::eff3770)
            );
    }

    double
    EEToCCBarHRvD::sigma_eetoD0Dbar0(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return _imp->sigma_eetochannel(ir, Channels::D0Dbar0);
    }

    double
    EEToCCBarHRvD::sigma_eetoDpDm(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        const double qed = 1.0 + _imp->alpha_em / M_PI * _imp->eta(real(ir->s));
        return _imp->sigma_eetochannel(ir, Channels::DpDm) * qed;
    }

    double
    EEToCCBarHRvD::rho_ee(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return _imp->rho(ir, Channels::ee);
    }

    double
    EEToCCBarHRvD::rho_eff(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return _imp->rho(ir, Channels::eff3770);
    }

    double
    EEToCCBarHRvD::rho_D0Dbar0(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return _imp->rho(ir, Channels::D0Dbar0);
    }

    double
    EEToCCBarHRvD::rho_DpDm(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return _imp->rho(ir, Channels::DpDm);
    }

    double
    EEToCCBarHRvD::re_chew_mandelstam_ee(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam(ir, Channels::ee));
    }

    double
    EEToCCBarHRvD::im_chew_mandelstam_ee(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam(ir, Channels::ee));
    }

    double
    EEToCCBarHRvD::re_chew_mandelstam_eff(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam(ir, Channels::eff3770));
    }

    double
    EEToCCBarHRvD::im_chew_mandelstam_eff(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam(ir, Channels::eff3770));
    }

    double
    EEToCCBarHRvD::re_chew_mandelstam_D0Dbar0(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBarHRvD::im_chew_mandelstam_D0Dbar0(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBarHRvD::re_chew_mandelstam_DpDm(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam(ir, Channels::DpDm));
    }

    double
    EEToCCBarHRvD::im_chew_mandelstam_DpDm(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam(ir, Channels::DpDm));
    }
    double
    EEToCCBarHRvD::re_chew_mandelstam_II_ee(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam_II(ir, Channels::ee));
    }

    double
    EEToCCBarHRvD::im_chew_mandelstam_II_ee(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam_II(ir, Channels::ee));
    }

    double
    EEToCCBarHRvD::re_chew_mandelstam_II_eff(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam_II(ir, Channels::eff3770));
    }

    double
    EEToCCBarHRvD::im_chew_mandelstam_II_eff(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam_II(ir, Channels::eff3770));
    }

    double
    EEToCCBarHRvD::re_chew_mandelstam_II_D0Dbar0(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam_II(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBarHRvD::im_chew_mandelstam_II_D0Dbar0(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam_II(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBarHRvD::re_chew_mandelstam_II_DpDm(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->chew_mandelstam_II(ir, Channels::DpDm));
    }

    double
    EEToCCBarHRvD::im_chew_mandelstam_II_DpDm(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->chew_mandelstam_II(ir, Channels::DpDm));
    }


    double
    EEToCCBarHRvD::re_T_eetoee(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->T_eetochannel(ir, Channels::ee));
    }

    double
    EEToCCBarHRvD::im_T_eetoee(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->T_eetochannel(ir, Channels::ee));
    }

    double
    EEToCCBarHRvD::re_T_eetoeff(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->T_eetochannel(ir, Channels::eff3770));
    }

    double
    EEToCCBarHRvD::im_T_eetoeff(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->T_eetochannel(ir, Channels::eff3770));
    }

    double
    EEToCCBarHRvD::re_T_eetoDpDm(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->T_eetochannel(ir, Channels::DpDm));
    }

    double
    EEToCCBarHRvD::im_T_eetoDpDm(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->T_eetochannel(ir, Channels::DpDm));
    }

    double
    EEToCCBarHRvD::re_T_eetoD0Dbar0(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->T_eetochannel(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBarHRvD::im_T_eetoD0Dbar0(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->T_eetochannel(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBarHRvD::re_T_II_eetoee(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->T_II_eetochannel(ir, Channels::ee));
    }

    double
    EEToCCBarHRvD::im_T_II_eetoee(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->T_II_eetochannel(ir, Channels::ee));
    }

    double
    EEToCCBarHRvD::re_T_II_eetoeff(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->T_II_eetochannel(ir, Channels::eff3770));
    }

    double
    EEToCCBarHRvD::im_T_II_eetoeff(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->T_II_eetochannel(ir, Channels::eff3770));
    }

    double
    EEToCCBarHRvD::re_T_II_eetoDpDm(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->T_II_eetochannel(ir, Channels::DpDm));
    }

    double
    EEToCCBarHRvD::im_T_II_eetoDpDm(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->T_II_eetochannel(ir, Channels::DpDm));
    }

    double
    EEToCCBarHRvD::re_T_II_eetoD0Dbar0(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return real(_imp->T_II_eetochannel(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBarHRvD::im_T_II_eetoD0Dbar0(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return imag(_imp->T_II_eetochannel(ir, Channels::D0Dbar0));
    }

    double
    EEToCCBarHRvD::psi3770_spectral_function(const double & E) const
    {
        return _imp->spectral_function(E, Resonances::psi3770);
    }


    double
    EEToCCBarHRvD::R(const EEToCCBarHRvD::IntermediateResult * ir) const
    {
        return _imp->R(ir);
    }

    const std::set<ReferenceName>
    EEToCCBarHRvD::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    EEToCCBarHRvD::begin_options()
    {
        return Implementation<EEToCCBarHRvD>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    EEToCCBarHRvD::end_options()
    {
        return Implementation<EEToCCBarHRvD>::options.cend();
    }
}
