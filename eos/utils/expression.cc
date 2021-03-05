#include <iostream>
#include "expression.hh"

namespace eos::ex
{

    double BinaryExpression::Add(double a, double b)      { return a+b; }
    double BinaryExpression::Subtract(double a, double b) { return a-b; }
    double BinaryExpression::Multuply(double a, double b) { return a*b; }
    double BinaryExpression::Divide(double a, double b)   { return a/b; }

    BinaryExpression::func
    BinaryExpression::Method(char op)
    {
        switch(op) {
            case '+': return BinaryExpression::Add;
            case '-': return BinaryExpression::Subtract;
            case '*': return BinaryExpression::Multuply;
            case '/': return BinaryExpression::Divide;
            default:  return nullptr;
        }
    }

    double BinaryExpression::Evaluate() const
    {
            func f = BinaryExpression::Method(_op);
            assert(f && _lhs && _rhs);
            return f(_lhs->Evaluate(), _rhs->Evaluate());
    }

    double ConstantExpression::Evaluate() const
    {
         return _value;
    }

    double VariableExpression::Evaluate() const
    {
        return get(_name);
    }


    //Print methods
    std::ostream&
    BinaryExpression::Print(std::ostream& os) const
    {
        return os << "BinaryExpression(" << *_lhs << " " << _op << " " << *_rhs << ")";
    }

    std::ostream&
    ConstantExpression::Print(std::ostream& os) const
    {
        return os << "ConstantExpression(" << _value << ")";
    }

    std::ostream&
    VariableExpression::Print(std::ostream& os) const
    {
        return os << "VariableExpression('" << _name << "')";
    }

    std::ostream&
    Expression::Print(std::ostream& os) const
    {
        return os << "Expression(" << *_e << ")";
    }

}
