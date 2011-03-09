/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <config.h>
#include <src/observable.hh>
#include <src/utils/cartesian-product.hh>
#include <src/utils/chi-squared.hh>
#include <src/utils/destringify.hh>
#include <src/utils/instantiation_policy-impl.hh>
#include <src/utils/one-of.hh>
#include <src/utils/power_of.hh>
#include <src/utils/scan_file.hh>
#include <src/utils/stringify.hh>
#include <src/utils/thread_pool.hh>
#include <src/utils/wilson-polynomial.hh>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <list>

using namespace eos;

class DoUsage
{
    private:
        std::string _what;

    public:
        DoUsage(const std::string & what) :
            _what(what)
        {
        }

        const std::string & what() const
        {
            return _what;
        }
};

struct ObservableInput
{
    ObservablePtr observable;

    double min, central, max;
};

struct ObservableRatioInput
{
    ObservablePtr numerator, denominator;

    double min, central, max;
};

typedef OneOf<ObservableInput, ObservableRatioInput> Input;

struct ScanData
{
    std::string name;

    unsigned points;

    double min;

    double max;
};

class CommandLine :
    public InstantiationPolicy<CommandLine, Singleton>
{
    public:
        Parameters parameters;

        std::list<ScanData> scans;

        std::list<std::string> coefficients;

        std::list<std::string> variations;

        std::list<Input> inputs;

        std::string output;

        CommandLine() :
            parameters(Parameters::Defaults())
        {
        }

        void parse(int argc, char ** argv)
        {
            std::shared_ptr<Kinematics> kinematics(new Kinematics);

            for (char ** a(argv + 1), ** a_end(argv + argc) ; a != a_end ; ++a)
            {
                std::string argument(*a);
                if ("--scan-abs" == argument)
                {
                    std::string coefficient = std::string(*(++a));
                    std::string name = "Abs{" + coefficient + "}";
                    unsigned points = destringify<unsigned>(*(++a));
                    double min = destringify<double>(*(++a));
                    double max = destringify<double>(*(++a));
                    scans.push_back(ScanData{name, points, min, max});
                    if (coefficients.cend() == std::find(coefficients.cbegin(), coefficients.cend(), coefficient))
                        coefficients.push_back(coefficient);

                    continue;
                }
                if ("--scan-arg" == argument)
                {
                    std::string coefficient = std::string(*(++a));
                    std::string name = "Arg{" + coefficient + "}";
                    unsigned points = destringify<unsigned>(*(++a));
                    double min = destringify<double>(*(++a));
                    double max = destringify<double>(*(++a));
                    scans.push_back(ScanData{name, points, min, max});
                    if (coefficients.cend() == std::find(coefficients.cbegin(), coefficients.cend(), coefficient))
                        coefficients.push_back(coefficient);

                    continue;
                }

                if ("--kinematics" == argument)
                {
                    std::string name = std::string(*(++a));
                    double value = destringify<double>(*(++a));
                    kinematics->declare(name);
                    kinematics->set(name, value);

                    continue;
                }

                if ("--observable" == argument)
                {
                    std::string observable_name(*(++a));

                    ObservableInput input;
                    input.observable = Observable::make(observable_name, parameters, *kinematics, Options());
                    if (! input.observable)
                        throw DoUsage("Unknown observable '" + observable_name + "'");

                    input.min = destringify<double>(*(++a));
                    input.central = destringify<double>(*(++a));
                    input.max = destringify<double>(*(++a));

                    inputs.push_back(input);
                    kinematics.reset(new Kinematics);

                    continue;
                }

                if ("--ratio" == argument)
                {
                    std::string numerator_name(*(++a)), denominator_name(*(++a));

                    ObservableRatioInput input;
                    input.numerator = Observable::make(numerator_name, parameters, *kinematics, Options());
                    if (! input.numerator)
                        throw DoUsage("Unknown observable '" + numerator_name + "'");

                    input.denominator = Observable::make(denominator_name, parameters, *kinematics, Options());
                    if (! input.numerator)
                        throw DoUsage("Unknown observable '" + denominator_name + "'");

                    input.min = destringify<double>(*(++a));
                    input.central = destringify<double>(*(++a));
                    input.max = destringify<double>(*(++a));

                    inputs.push_back(input);
                    kinematics.reset(new Kinematics);

                    continue;
                }

                if ("--vary" == argument)
                {
                    std::string variation_name(*(++a));
                    Parameter p = parameters[variation_name];
                    variations.push_back(variation_name);

                    continue;
                }

                if ("--output" == argument)
                {
                    output = std::string(*(++a));
                    continue;
                }

                throw DoUsage("Unknown command line argument: " + argument);
            }
        }
};

class WilsonScannerPolynomial
{
    private:
        CartesianProduct<std::vector<double>> _points;

        std::vector<Parameter> _scan_parameters;

        std::vector<Parameter> _variations;

        std::vector<std::tuple<ObservablePtr, double, double, double>> _observables;

        std::vector<Ticket> _tickets;

        ScanFile _output;

        std::vector<ScanFile::DataSet> _data_sets;

    public:
        WilsonScannerPolynomial() :
            _output(ScanFile::Create(CommandLine::instance()->output, "eos-scan-polynomial"))
        {
            std::cout << std::scientific;
            std::cout << "# Scan generated by eos-scan-polynomial (" EOS_GITHEAD ")" << std::endl;
            std::cout << "# Coefficients:" << std::endl;
            for (auto c = CommandLine::instance()->coefficients.cbegin(), c_end = CommandLine::instance()->coefficients.cend() ; c != c_end ; ++c)
            {
                std::cout << "#   " << *c << std::endl;
            }

            std::cout << "# Scans:" << std::endl;
            CartesianProduct<std::vector<double>> cp;
            for (auto s = CommandLine::instance()->scans.cbegin(), s_end = CommandLine::instance()->scans.cend() ; s != s_end ; ++s)
            {
                double delta = (s->max - s->min) / s->points;

                std::vector<double> scan_set;
                for (unsigned i = 0 ; i <= s->points ; ++i)
                {
                    scan_set.push_back(s->min + delta * i);
                }

                _points.over(scan_set);
                _scan_parameters.push_back(CommandLine::instance()->parameters[s->name]);

                std::cout << "#   " << s->name << ": [" << s->min << ", " << s->max << "], increment = " << delta << std::endl;
            }

            std::cout << "# Inputs:" << std::endl;
            for (auto i = CommandLine::instance()->inputs.cbegin(), i_end = CommandLine::instance()->inputs.cend() ; i != i_end ; ++i)
            {
                i->accept(*this);
            }

            std::cout << "# Variations:" << std::endl;
            for (auto v = CommandLine::instance()->variations.cbegin(), v_end = CommandLine::instance()->variations.cend() ; v != v_end ; ++v)
            {
                Parameter p = CommandLine::instance()->parameters[*v];
                std::cout << "#   " << p.name() << ": " << p.min() << " < " << p << " < " << p.max() << std::endl;
                _variations.push_back(p);
            }
        }

        void visit(const ObservableInput & i)
        {
            Parameters parameters = CommandLine::instance()->parameters;

            std::cout << "#   " << i.observable->name() << '[' << i.observable->kinematics().as_string() << "] = (" << i.min << ", " << i.central << ", " << i.max << ")" << std::endl;

            ObservablePtr observable = make_polynomial_observable(make_polynomial(i.observable, CommandLine::instance()->coefficients), parameters);
            _observables.push_back(std::make_tuple(observable, i.min, i.central, i.max));
        }

        void visit(const ObservableRatioInput & i)
        {
            Parameters parameters = CommandLine::instance()->parameters;

            std::cout << "#   " << i.numerator->name() << "[Kinematics]" << " / " << i.denominator->name() << "[Kinematics]" << " = (" << i.min << ", " << i.central << ", " << i.max << ")" << std::endl;

            ObservablePtr observable = make_polynomial_ratio(make_polynomial(i.numerator, CommandLine::instance()->coefficients),
                    make_polynomial(i.denominator, CommandLine::instance()->coefficients),
                    parameters);
            _observables.push_back(std::make_tuple(observable, i.min, i.central, i.max));
        }

        void scan_range(const CartesianProduct<std::vector<double>>::Iterator & begin, const CartesianProduct<std::vector<double>>::Iterator & end,
                unsigned index)
        {
            // Clone
            Parameters parameters = CommandLine::instance()->parameters.clone();

            std::vector<Parameter> scan_parameters;
            for (auto p = _scan_parameters.cbegin(), p_end = _scan_parameters.cend() ; p != p_end ; ++p)
            {
                scan_parameters.push_back(parameters[p->name()]);
            }

            std::vector<std::tuple<ObservablePtr, double, double, double>> observables;
            for (auto o = _observables.cbegin(), o_end = _observables.cend() ; o != o_end ; ++o)
            {
                observables.push_back(std::make_tuple(std::get<0>(*o)->clone(parameters), std::get<1>(*o), std::get<2>(*o), std::get<3>(*o)));
            }

            std::vector<Parameter> variations;
            for (auto v = _variations.cbegin(), v_end = _variations.cend() ; v != v_end ; ++v)
            {
                variations.push_back(parameters[v->name()]);
            }

            // Scan our range
            for (auto i = begin, i_end = end ; i != i_end ; ++i)
            {
                std::vector<double> v = *i;
                std::vector<double> result;

                auto p = scan_parameters.begin();
                for (auto j = v.begin(), j_end = v.end() ; j != j_end ; ++j, ++p)
                {
                    (*p) = *j;
                    result.push_back(*j);
                }

                // Calculate chi^2
                double chi_squared = 0.0;
                for (auto o = observables.begin(), o_end = observables.end() ; o != o_end ; ++o)
                {
                    ObservablePtr observable = std::get<0>(*o);
                    double central = observable->evaluate();
                    double delta_min = 0.0, delta_max = 0.0;

                    for (auto v = variations.begin(), v_end = variations.end() ; v != v_end ; ++v)
                    {
                        double old_v = *v;
                        double max = 0.0, min = 0.0, value;

                        // Handle parameters with lowered values
                        *v = v->min();
                        value = observable->evaluate();

                        if (value > central)
                            max = value - central;

                        if (value < central)
                            min = central - value;

                        // Handle parameters with raised values
                        *v = v->max();
                        value = observable->evaluate();

                        if (value > central)
                            max = std::max(max, value - central);

                        if (value < central)
                            min = std::max(min, central -value);

                        *v = old_v;

                        delta_min += min * min;
                        delta_max += max * max;
                    }

                    chi_squared += ChiSquared::with_theory_offset(central - std::sqrt(delta_min), central, central + std::sqrt(delta_max),
                            std::get<1>(*o), std::get<2>(*o), std::get<3>(*o));
                }

                result.push_back(chi_squared);

                _data_sets[index] << result;
            }
        }

        void scan()
        {
            // Find chunk size for N chunks
            unsigned chunk_size = _points.size() / ThreadPool::instance()->number_of_threads();

            // Enqueue N - 1 jobs
            auto c = _points.begin();
            unsigned i = 0;
            for ( ; i < ThreadPool::instance()->number_of_threads() - 1 ; ++i, c += chunk_size)
            {
                auto begin = c, end = c;
                end += chunk_size;

                _data_sets.push_back(_output.add("chunk #" + stringify(i), _scan_parameters.size() + 1));
                _tickets.push_back(ThreadPool::instance()->enqueue(std::bind(&WilsonScannerPolynomial::scan_range, this, begin, end, i)));
            }

            // Enqueue an extra job for the remains
            _data_sets.push_back(_output.add("chunk #" + stringify(i), _scan_parameters.size() + 1));
            _tickets.push_back(ThreadPool::instance()->enqueue(std::bind(&WilsonScannerPolynomial::scan_range, this, c, _points.end(), i)));

            // Wait for job completion
            for (auto t = _tickets.begin(), t_end = _tickets.end() ; t != t_end ; ++t)
            {
                t->wait();
            }
        }
};

int
main(int argc, char * argv[])
{
    try
    {
        CommandLine::instance()->parse(argc, argv);

        if (CommandLine::instance()->inputs.empty())
            throw DoUsage("Need to specify at least one input!");

        if (CommandLine::instance()->output.empty())
            throw DoUsage("Need to specify output!");

        WilsonScannerPolynomial scanner;
        scanner.scan();
    }
    catch(DoUsage & e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-scan-polynomial" << std::endl;
        std::cout << "  [--vary PARAMETER]*" << std::endl;
        std::cout << "  [[--kinematics NAME VALUE]* --observable NAME MIN CENTRAL MAX]+" << std::endl;
        std::cout << "  [[--scan-abs COEFFICIENT POINTS MIN MAX] | [--scan-arg COEFFICIENT POINTS MIN MAX]]+" << std::endl;
    }
    catch(Exception & e)
    {
        std::cerr << "Caught exception: '" << e.what() << "'" << std::endl;
        return EXIT_FAILURE;
    }
    catch(...)
    {
        std::cerr << "Aborting after unknown exception" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
