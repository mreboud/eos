#ifndef EXPRESSION_HH
#define EXPRESSION_HH 1

#include <cassert>
#include <memory>
#include <iostream>
#include <map>

namespace eos::ex
{

    struct AbstractExpression;
    typedef std::shared_ptr<AbstractExpression> Ptr;

    struct AbstractExpression {
        virtual ~AbstractExpression() {}
        virtual double Evaluate() const = 0;
        virtual std::ostream& Print(std::ostream& os) const = 0;

        friend std::ostream& operator<<(std::ostream& os, AbstractExpression const& e)
            { return e.Print(os); }

        protected: AbstractExpression() {}
    };



    template <typename Expr> // general purpose, static Expression cloner
        static Ptr make_from(Expr const& t) { return std::make_shared<Expr>(t); }



    struct BinaryExpression : AbstractExpression
    {
        BinaryExpression() {}

        template<typename L, typename R>
        BinaryExpression(char op, L const& l, R const& r)
            : _op(op), _lhs(make_from(l)), _rhs(make_from(r))
        {}

        double Evaluate() const;

      private:
        char _op;
        Ptr _lhs, _rhs;

        typedef double(*func)(double, double);

        static double Add(double a, double b);
        static double Subtract(double a, double b);
        static double Multuply(double a, double b);
        static double Divide(double a, double b);

        static BinaryExpression::func Method(char op);

        std::ostream& Print(std::ostream& os) const;
    };



    struct ConstantExpression : AbstractExpression {

        double _value;

        ConstantExpression(double v = 0) : _value(v) {}

        double Evaluate() const;

        virtual std::ostream& Print(std::ostream& os) const;
    };



    struct VariableExpression : AbstractExpression {

        std::string _name;

        static double& get(std::string const& name) {
            static std::map<std::string, double> _symbols;
            return _symbols[name];
            /*switch(name) {
             *    case 'a': static double a; return a;
             *    case 'b': static double b; return b;
             *    default:  throw "undefined variable";
             *}
             */
        }

        double Evaluate() const;

        virtual std::ostream& Print(std::ostream& os) const;
    };



    struct Expression : AbstractExpression
    {
        Expression() { }

        template <typename E>
        Expression(E const& e) : _e(make_from(e)) { } // cloning the expression

        double Evaluate() const { assert(_e); return _e->Evaluate(); }

        // special purpose overload to avoid unnecessary wrapping
        friend Ptr make_from(Expression const& t) { return t._e; }
      private:
        Ptr _e;
        virtual std::ostream& Print(std::ostream& os) const;
    };
}

#endif
