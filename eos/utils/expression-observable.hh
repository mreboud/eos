#ifndef EXPRESSION_OBSERVABLE_HH
#define EXPRESSION_OBSERVABLE_HH 1

#include <eos/observable.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>

#include <string>


namespace eos
{

    class ExpressionObservableEntry : ObservableEntry
    {
    public:

        ExpressionObservableEntry(std::string input);

        ObservablePtr make(const Parameters &, const Kinematics &, const Options &) const
        {
            return nullptr;
        }

        ~ExpressionObservableEntry();

    };

}

#endif
