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
        UsedParameter hbar;
        UsedParameter alpha_em;
        UsedParameter m_e;
        UsedParameter m_D0;
        UsedParameter m_D;
        UsedParameter m_Dst0;
        UsedParameter m_Dst;
        UsedParameter m_Ds;
        UsedParameter m_Dsst;

        static const inline std::vector<std::string> resonance_names =
        {
            "psi(2S)", "psi(3770)", "psi(4040)", "psi(4160)", "psi(4415)"
        };

        enum Resonances
        {
            psi2S = 0, psi3770, psi4040, psi4160, psi4415
        };

        static const inline std::vector<std::string> channel_names =
        {
            "e^+e^-", "eff", "D^0Dbar^0", "D^+D^-", "D^0Dbar^*0", "D^*0Dbar^0",
            "D^+D^*-", "D^*+D^-", "D_s^+D_s^-", "D^*0Dbar^*0P0", "D^*0Dbar^*0P2", "D^*0Dbar^*0F2", "D^*+D^*-P0",
            "D^*+D^*-P2", "D^*+D^*-F2", "D_s^+D_s^*-", "D_s^*+D_s^-", "D_s^*+D_s^*-P0", "D_s^*+D_s^*-P2",
            "D_s^*+D_s^*-F2"
        };

        enum Channels
        {
            ee = 0, eff, D0Dbar0, DpDm, D0Dbarst0, Dst0Dbar0,
            DpDstm, DstpDm, DspDsm, Dst0Dbarst0P0, Dst0Dbarst0P2, Dst0Dbarst0F2, DstpDstmP0,
            DstpDstmP2, DstpDstmF2, DspDsstm, DsstpDsm, DsstpDsstmP0, DsstpDsstmP2,
            DsstpDsstmF2
        };

        // Charmonium masses and effective sizes
        std::array<UsedParameter, EEToCCBar::nresonances> m;
        std::array<UsedParameter, EEToCCBar::nresonances> q;

        // Channel-Resonance couplings
        std::array<std::array<UsedParameter, EEToCCBar::nchannels>, EEToCCBar::nresonances> g0;

        // Non-cc contribution to the Rc ratio
        std::array<std::array<UsedParameter, EEToCCBar::nchannels>, EEToCCBar::nchannels> bkgcst;

        // R_uds
        UsedParameter Rconstant;

        // Normalization of the exclusive channels
        UsedParameter exclusive_norm;

        std::shared_ptr<KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances>> K;

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
                UsedParameter(p[resonance_names[indices] + "::q_R"], u)...
            }};
        }

        long unsigned _filter_channel_index(Channels channel)
        {
            switch (channel)
            {
                case DpDm:
                    return D0Dbar0;

                case Dst0Dbar0:
                case DpDstm:
                case DstpDm:
                    return D0Dbarst0;

                case Dst0Dbarst0F2:
                    return Dst0Dbarst0P2;

                case DstpDstmP0:
                    return Dst0Dbarst0P0;

                case DstpDstmP2:
                case DstpDstmF2:
                    return Dst0Dbarst0P2;

                case DsstpDsm:
                    return DspDsstm;

                case DsstpDsstmF2:
                    return DsstpDsstmP2;

                default:
                    return channel;
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
                UsedParameter(p["ee->ccbar::g0(" + resonance_names[row_index] + "," + channel_names[_filter_channel_index(Channels(column_indices))] + ")"], u)...
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
                UsedParameter(p["ee->ccbar::c(" + channel_names[_filter_channel_index(Channels(row_index))] + ","
                        + channel_names[_filter_channel_index(Channels(column_indices))] + ")"], u)...
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

            m(_resonance_masses(p, u, std::make_index_sequence<EEToCCBar::nresonances>())),
            q(_resonance_sizes(p, u, std::make_index_sequence<EEToCCBar::nresonances>())),
            g0(_g0_matrix(p, u, std::make_index_sequence<EEToCCBar::nresonances>(), std::make_index_sequence<EEToCCBar::nchannels>())),
            bkgcst(_c_matrix(p, u, std::make_index_sequence<EEToCCBar::nchannels>(), std::make_index_sequence<EEToCCBar::nchannels>())),

            Rconstant(p["ee->ccbar::Rconstant"], u),
            exclusive_norm(p["ee->ccbar::exclusive_norm"], u)
        {
            std::array<std::shared_ptr<KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances>::Resonance>, EEToCCBar::nresonances> resonance_array;
            for (unsigned i = 0; i<EEToCCBar::nresonances; i++)
            {
                resonance_array[i] = std::make_shared<CharmoniumResonance<EEToCCBar::nchannels, EEToCCBar::nresonances>>(resonance_names[i], m[i], q[i]);
            }

            std::array<std::shared_ptr<KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances>::Channel>, EEToCCBar::nchannels> channel_array;
            for (unsigned i = 0; i<EEToCCBar::nchannels; i++)
            {
                switch (Channels(i))
                {
                    case ee:
                    case eff:
                        channel_array[i] = std::make_shared<EffChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_e, m_e, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case D0Dbar0:
                        channel_array[i] = std::make_shared<PWavePPChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_D0, m_D0, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case DpDm:
                        channel_array[i] = std::make_shared<PWavePPChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_D, m_D, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case D0Dbarst0:
                    case Dst0Dbar0:
                        channel_array[i] = std::make_shared<PWaveVPChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_D0, m_Dst0, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case DpDstm:
                    case DstpDm:
                        channel_array[i] = std::make_shared<PWaveVPChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_D, m_Dst, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case DspDsm:
                        channel_array[i] = std::make_shared<PWavePPChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_Ds, m_Ds, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case Dst0Dbarst0P0:
                    case Dst0Dbarst0P2:
                        channel_array[i] = std::make_shared<PWaveVVChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_Dst0, m_Dst0, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case Dst0Dbarst0F2:
                        channel_array[i] = std::make_shared<FWaveVVChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_Dst0, m_Dst0, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case DstpDstmP0:
                    case DstpDstmP2:
                        channel_array[i] = std::make_shared<PWaveVVChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_Dst, m_Dst, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case DstpDstmF2:
                        channel_array[i] = std::make_shared<FWaveVVChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_Dst, m_Dst, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case DspDsstm:
                    case DsstpDsm:
                        channel_array[i] = std::make_shared<PWaveVPChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_Ds, m_Dsst, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case DsstpDsstmP0:
                    case DsstpDsstmP2:
                        channel_array[i] = std::make_shared<PWaveVVChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_Dsst, m_Dsst, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case DsstpDsstmF2:
                        channel_array[i] = std::make_shared<FWaveVVChannel<EEToCCBar::nchannels, EEToCCBar::nresonances>>(channel_names[i], m_Dsst, m_Dsst, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;

                    default:
                        throw InternalError("The number of requested channels (array) is larger than the number of known channel names.");
                }
            }

            K = std::shared_ptr<KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances>> (
                new KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances>(
                    channel_array,
                    resonance_array,
                    _get_bkgcst(std::make_index_sequence<EEToCCBar::nchannels>(),std::make_index_sequence<EEToCCBar::nchannels>()),
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

        double sigma_eetochannel(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            // Conversion factor between GeV^2 and nb
            const double speedoflight = 299792458.; //Exact value
            const double GeVtonb = 10 * pow( 1.e18 * hbar * speedoflight, 2.0);

            // Channel properties
            const double Nf = 2 * K->_channels[Channels(channel)]->_l_orbital + 1;
            const double rhof = real(K->_channels[Channels(channel)]->rho(intermediate_result->s));

            // Get T-matrix[ee, channel]
            const complex<double> T1f = intermediate_result->tmatrix_row_0[Channels(channel)];

            return GeVtonb * 16. * M_PI / intermediate_result->s * Nf * rhof * norm(T1f);
        }

        // Widths
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
            double total_xsec = 0.0;

            for (unsigned i = 1; i < EEToCCBar::nchannels; i++)
            {
                total_xsec += sigma_eetochannel(intermediate_result, Channels(i));
            }

            return total_xsec / sigma_eetomumu(intermediate_result->E) + Rconstant; //Add constant term
        }
    };

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

    using Resonances = Implementation<eos::EEToCCBar>::Resonances;
    using Channels = Implementation<eos::EEToCCBar>::Channels;

    double
    EEToCCBar::psi2S_ee_width() const
    {
        return _imp->res_partial_width(Resonances::psi2S, Channels::ee);
    }

    double
    EEToCCBar::psi2S_eff_width() const
    {
        return _imp->res_partial_width(Resonances::psi2S, Channels::eff);
    }

    double
    EEToCCBar::psi2S_total_width() const
    {
        return _imp->res_total_width(Resonances::psi2S);
    }

    double
    EEToCCBar::psi3770_total_width() const
    {
        return _imp->res_total_width(Resonances::psi3770);
    }

    double
    EEToCCBar::psi4040_total_width() const
    {
        return _imp->res_total_width(Resonances::psi4040);
    }

    double
    EEToCCBar::psi4160_total_width() const
    {
        return _imp->res_total_width(Resonances::psi4160);
    }

    double
    EEToCCBar::psi4415_total_width() const
    {
        return _imp->res_total_width(Resonances::psi4415);
    }


    double
    EEToCCBar::sigma_eetoee(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, Channels::ee);
    }

    double
    EEToCCBar::sigma_eetoeff(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, Channels::eff);
    }

    double
    EEToCCBar::sigma_eetoD0Dbar0(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, Channels::D0Dbar0);
    }

    double
    EEToCCBar::sigma_eetoDpDm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, Channels::DpDm);
    }

    double
    EEToCCBar::sigma_eetoD0Dbarst0(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, Channels::D0Dbarst0)
            + _imp->sigma_eetochannel(ir, Channels::Dst0Dbar0)
            );
    }

    double
    EEToCCBar::sigma_eetoDpDstm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, Channels::DpDstm)
            + _imp->sigma_eetochannel(ir, Channels::DstpDm)
            );
    }

    double
    EEToCCBar::sigma_eetoDspDsm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, Channels::DspDsm);
    }

    double
    EEToCCBar::sigma_eetoDst0Dbarst0(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, Channels::Dst0Dbarst0P0)
            + _imp->sigma_eetochannel(ir, Channels::Dst0Dbarst0P2)
            + _imp->sigma_eetochannel(ir, Channels::Dst0Dbarst0F2)
            );
    }

    double
    EEToCCBar::sigma_eetoDstpDstm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, Channels::DstpDstmP0)
            + _imp->sigma_eetochannel(ir, Channels::DstpDstmP2)
            + _imp->sigma_eetochannel(ir, Channels::DstpDstmF2)
            );
    }

    double
    EEToCCBar::sigma_eetoDspDsstm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, Channels::DspDsstm)
            + _imp->sigma_eetochannel(ir, Channels::DsstpDsm));
    }

    double
    EEToCCBar::sigma_eetoDsstpDsstm(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, Channels::DsstpDsstmP0)
            + _imp->sigma_eetochannel(ir, Channels::DsstpDsstmP2)
            + _imp->sigma_eetochannel(ir, Channels::DsstpDsstmF2)
            );
    }

    double
    EEToCCBar::R(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->R(ir);
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
