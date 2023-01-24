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
    struct Implementation<EEToCCBar>
    {
        UsedParameter hbar;
        UsedParameter alpha_em;
        UsedParameter m_e;
        UsedParameter m_D0;
        UsedParameter m_D;

        bool assume_isospin;

        static const inline std::vector<std::string> resonance_names =
        {
            "J/psi", "psi(2S)", "psi(3770)"
        };

        enum Resonances
        {
            Jpsi = 0, psi2S, psi3770
        };

        static const inline std::vector<std::string> channel_names =
        {
            "e^+e^-", "eff(Jpsi)", "Jpsipipi", "eff(2S)", "D^0Dbar^0", "D^+D^-", "eff(3770)"
        };

        enum Channels
        {
            ee = 0, effJpsi, Jpsipipi, eff2S, D0Dbar0, DpDm, eff3770
        };

        // Charmonium masses and effective sizes
        std::array<UsedParameter, EEToCCBar::nresonances> m;
        std::array<UsedParameter, EEToCCBar::nresonances> q;

        // Channel-Resonance couplings
        std::array<std::array<UsedParameter, EEToCCBar::nchannels>, EEToCCBar::nresonances> g0;

        // Non-cc contribution to the Rc ratio
        std::array<std::array<std::array<UsedParameter, EEToCCBar::nchannels>, EEToCCBar::nchannels>, EEToCCBar::order + 1> bkgpol;

        // R_uds
        UsedParameter Rconstant;

        // Normalization of the exclusive channels
        UsedParameter exclusive_norm;

        std::shared_ptr<KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances, EEToCCBar::order>> K;

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
                UsedParameter(p["q_R::" + resonance_names[indices]], u)...
            }};
        }

        long unsigned _filter_channel_index(Channels channel)
        {
            if (assume_isospin)
            {
                switch (channel)
                {
                    case DpDm:
                        return D0Dbar0;

                    default:
                        return channel;
                }
            }
            else
            {
                switch (channel)
                {
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
        auto _c_row(const Parameters & p, ParameterUser & u, T power, T row_index, std::integer_sequence<T, column_indices...>)
            -> std::array<UsedParameter, sizeof...(column_indices)>
        {
            if (sizeof...(column_indices) > channel_names.size())
                throw InternalError("The number of requested channels (column) is larger than the number of channel names.");

            return std::array<UsedParameter, sizeof...(column_indices)>
            {{
                UsedParameter(p["ee->ccbar::c_" + std::to_string(power) + "(" + channel_names[_filter_channel_index(Channels(row_index))] + ","
                        + channel_names[_filter_channel_index(Channels(column_indices))] + ")"], u)...
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
            if (sizeof...(power_indices) >EEToCCBar::order + 1)
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

            assume_isospin(destringify<bool>(o.get("assume_isospin", "false"))),

            m(_resonance_masses(p, u, std::make_index_sequence<EEToCCBar::nresonances>())),
            q(_resonance_sizes(p, u, std::make_index_sequence<EEToCCBar::nresonances>())),
            g0(_g0_matrix(p, u, std::make_index_sequence<EEToCCBar::nresonances>(), std::make_index_sequence<EEToCCBar::nchannels>())),
            bkgpol(_c_array(p, u, std::make_index_sequence<EEToCCBar::order + 1>(), std::make_index_sequence<EEToCCBar::nchannels>(), std::make_index_sequence<EEToCCBar::nchannels>())),

            Rconstant(p["ee->ccbar::Rconstant"], u),
            exclusive_norm(p["ee->ccbar::exclusive_norm"], u)
        {
            std::array<std::shared_ptr<KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances, EEToCCBar::order>::Resonance>, EEToCCBar::nresonances> resonance_array;
            for (unsigned i = 0; i<EEToCCBar::nresonances; i++)
            {
                resonance_array[i] = std::make_shared<CharmoniumResonance<EEToCCBar::nchannels, EEToCCBar::nresonances, EEToCCBar::order>>(resonance_names[i], m[i], q[i]);
            }

            std::array<std::shared_ptr<KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances, EEToCCBar::order>::Channel>, EEToCCBar::nchannels> channel_array;
            for (unsigned i = 0; i<EEToCCBar::nchannels; i++)
            {
                switch (Channels(i))
                {
                    case ee:
                    case effJpsi:
                        channel_array[i] = std::make_shared<EffChannel<EEToCCBar::nchannels, EEToCCBar::nresonances, EEToCCBar::order>>(channel_names[i], m_e, m_e, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case Jpsipipi:
                        channel_array[i] = std::make_shared<EffChannel<EEToCCBar::nchannels, EEToCCBar::nresonances, EEToCCBar::order>>(channel_names[i], m[0], m_e, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case eff2S:
                    case eff3770:
                        channel_array[i] = std::make_shared<EffChannel<EEToCCBar::nchannels, EEToCCBar::nresonances, EEToCCBar::order>>(channel_names[i], m_e, m_e, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case D0Dbar0:
                        channel_array[i] = std::make_shared<PWavePPChannel<EEToCCBar::nchannels, EEToCCBar::nresonances, EEToCCBar::order>>(channel_names[i], m_D0, m_D0, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;
                    case DpDm:
                        channel_array[i] = std::make_shared<PWavePPChannel<EEToCCBar::nchannels, EEToCCBar::nresonances, EEToCCBar::order>>(channel_names[i], m_D, m_D, _get_g0_column(_filter_channel_index(Channels(i)), std::make_index_sequence<EEToCCBar::nresonances>()));
                        break;

                    default:
                        throw InternalError("The number of requested channels (array) is larger than the number of known channel names.");
                }
            }

            K = std::shared_ptr<KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances, EEToCCBar::order>> (
                new KMatrix<EEToCCBar::nchannels, EEToCCBar::nresonances, EEToCCBar::order>(
                    channel_array,
                    resonance_array,
                    _get_bkgpol(std::make_index_sequence<EEToCCBar::order + 1>(), std::make_index_sequence<EEToCCBar::nchannels>(), std::make_index_sequence<EEToCCBar::nchannels>()),
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
        complex<double> amplitude_eetochannel(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            // Conversion factor between GeV^2 and nb
            const double speedoflight = 299792458.; //Exact value
            const double GeVtonb = 10 * pow( 1.e18 * hbar * speedoflight, 2.0);

            // Channel properties
            const double Nf = 2 * K->_channels[Channels(channel)]->_l_orbital + 1;
            const double rhof = real(K->_channels[Channels(channel)]->rho(intermediate_result->s));

            // Get T-matrix[ee, channel]
            const complex<double> T1f = intermediate_result->tmatrix_row_0[Channels(channel)];

            return sqrt(GeVtonb * 16. * M_PI / intermediate_result->s * Nf * rhof) * T1f;
        }

        double sigma_eetochannel(const IntermediateResult * intermediate_result, const Channels & channel)
        {
            return norm(amplitude_eetochannel(intermediate_result, channel));
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
        {"assume_isospin", { "true", "false" }, "false"},
    };

    const EEToCCBar::IntermediateResult *
    EEToCCBar::prepare(const double & E) const
    {
        return _imp->prepare(E);
    }

    using Resonances = Implementation<eos::EEToCCBar>::Resonances;
    using Channels = Implementation<eos::EEToCCBar>::Channels;

    double
    EEToCCBar::Jpsi_ee_width() const
    {
        return _imp->res_partial_width(Resonances::Jpsi, Channels::ee);
    }

    double
    EEToCCBar::Jpsi_eff_width() const
    {
        return _imp->res_partial_width(Resonances::Jpsi, Channels::eff2S);
    }

    double
    EEToCCBar::Jpsi_total_width() const
    {
        return _imp->res_total_width(Resonances::Jpsi);
    }

    double
    EEToCCBar::psi2S_ee_width() const
    {
        return _imp->res_partial_width(Resonances::psi2S, Channels::ee);
    }

    double
    EEToCCBar::psi2S_eff_width() const
    {
        return _imp->res_partial_width(Resonances::psi2S, Channels::eff2S);
    }

    double
    EEToCCBar::psi2S_Jpsipipi_width() const
    {
        return _imp->res_partial_width(Resonances::psi2S, Channels::Jpsipipi);
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
    EEToCCBar::psi3770_D0Dbar0_width() const
    {
        return _imp->res_partial_width(Resonances::psi3770, Channels::D0Dbar0);
    }

    double
    EEToCCBar::psi3770_DpDm_width() const
    {
        return _imp->res_partial_width(Resonances::psi3770, Channels::DpDm);
    }

    double
    EEToCCBar::psi3770_eff_width() const
    {
        return _imp->res_partial_width(Resonances::psi3770, Channels::eff3770);
    }

    double
    EEToCCBar::psi3770_Jpsipipi_width() const
    {
        return _imp->res_partial_width(Resonances::psi3770, Channels::Jpsipipi);
    }

    double
    EEToCCBar::sigma_eetoee(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, Channels::ee);
    }

    double
    EEToCCBar::sigma_eetoeff(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * (
              _imp->sigma_eetochannel(ir, Channels::effJpsi)
            + _imp->sigma_eetochannel(ir, Channels::eff2S)
            + _imp->sigma_eetochannel(ir, Channels::eff3770)
            );
    }

    double
    EEToCCBar::sigma_eetoJpsipipi(const EEToCCBar::IntermediateResult * ir) const
    {
        return _imp->exclusive_norm * _imp->sigma_eetochannel(ir, Channels::Jpsipipi);
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
