#ifndef EXPRESSION_PARSER_HH
#define EXPRESSION_PARSER_HH 1

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/adapted.hpp>
#include "expression.hh"

namespace qi    = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phx   = boost::phoenix;

namespace eos
{

    template <typename Iterator>
    struct ExpressionParser : qi::grammar<Iterator, eos::ex::Expression(), ascii::space_type>
    {
        struct MakeBinaryExpression {
            template<typename,typename,typename> struct result { typedef eos::ex::BinaryExpression type; };

            template<typename C, typename L, typename R>
                eos::ex::BinaryExpression operator()(C op, L const& lhs, R const& rhs) const
                { return eos::ex::BinaryExpression(op, lhs, rhs); }
        };

        phx::function<MakeBinaryExpression> makebinary;

        // Constructor
        ExpressionParser();

        qi::rule<Iterator, eos::ex::Expression()        , ascii::space_type> expression;
        qi::rule<Iterator, eos::ex::Expression()        , ascii::space_type> additive_expr;
        qi::rule<Iterator, eos::ex::Expression()        , ascii::space_type> non_additive_expr;

        qi::rule<Iterator, eos::ex::Expression()        , ascii::space_type> primary_expr;
        qi::rule<Iterator, eos::ex::ConstantExpression(), ascii::space_type> constant;
        qi::rule<Iterator, eos::ex::VariableExpression(), ascii::space_type> variable;
        qi::rule<Iterator, std::string()                , ascii::space_type> string;

        // Destuctor
        ~ExpressionParser();
    };

}
#endif