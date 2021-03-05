#ifndef EXPRESSION_PARSER_IMPL_HH
#define EXPRESSION_PARSER_IMPL_HH 1

// #define BOOST_SPIRIT_DEBUG
// #define BOOST_RESULT_OF_USE_DECLTYPE
// #define BOOST_SPIRIT_USE_PHOENIX_V3

#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/adapted.hpp>
#include "expression.hh"
#include "expression-parser.hh"

BOOST_FUSION_ADAPT_STRUCT(eos::ex::VariableExpression, (std::string, _name))

namespace eos
{

    template <typename Iterator>
    ExpressionParser<Iterator>::ExpressionParser() : ExpressionParser::base_type(expression)
        {
            using namespace boost::spirit::qi;

            expression =
                additive_expr                                [ _val = _1]
                ;

            additive_expr =
                non_additive_expr                            [ _val = _1]
                >> *(  (char_("+") >> non_additive_expr)     [ _val = makebinary(_1, _val, _2)]
                   |(  char_("-") >> non_additive_expr)      [ _val = makebinary(_1, _val, _2)]
                   )
                ;

            non_additive_expr =
                primary_expr                                 [ _val = _1]
                >> *( (char_('*') >> primary_expr)           [ _val = makebinary(_1, _val, _2)]
                   |( char_('/') >> primary_expr)            [ _val = makebinary(_1, _val, _2)]
                   )
                ;


            primary_expr =
                  ( '(' >> expression >> ')' )               [ _val = _1 ]
                | constant                                   [ _val = _1 ]
                | variable                                   [ _val = _1 ]
                ;

            constant = lexeme ["0x" >> hex] | double_ | int_;
            string   = '{' >> lexeme [ *~char_("}") ] > '}';
            variable = string | as_string [ alpha ];

            BOOST_SPIRIT_DEBUG_NODE(expression);
            BOOST_SPIRIT_DEBUG_NODE(additive_expr);
            BOOST_SPIRIT_DEBUG_NODE(non_additive_expr);
            BOOST_SPIRIT_DEBUG_NODE(primary_expr);
            BOOST_SPIRIT_DEBUG_NODE(constant);
            BOOST_SPIRIT_DEBUG_NODE(variable);
            BOOST_SPIRIT_DEBUG_NODE(string);
        }

    template <typename Iterator>
    ExpressionParser<Iterator>::~ExpressionParser(){}
}

#endif
