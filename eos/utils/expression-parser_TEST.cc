#include <test/test.hh>
#include <eos/utils/expression-parser-impl.hh>


using namespace test;
using namespace eos::ex;
using namespace eos;


struct ExpressionTest
{
    const std::string _input;
    const double _a, _b;

    Expression _e;
    bool _completed;
    double _result;

    typedef std::string::const_iterator It;
    ExpressionParser<It> p;

    ExpressionTest(const std::string input, double a=0, double b=0) :
    _input(input), _a(a), _b(b)
    {
        It f(_input.begin()), l(_input.end());
        _completed = qi::phrase_parse(f, l, p, ascii::space, _e);
        VariableExpression::get("a") = _a;
        VariableExpression::get("b") = _b;
        _result = _e.Evaluate();
    }
};


class ExpressionParserTest :
public TestCase
{
public:
    ExpressionParserTest() :
    TestCase("Expression parser tests")
    {
    }

    virtual void run() const
    {

        ExpressionTest test("1+2*3");

        TEST_CHECK(test._completed);
        TEST_CHECK_EQUAL(test._result, 7);

        ExpressionTest test2("(2 * a) + b", -3, 5);

        TEST_CHECK(test2._completed);
        TEST_CHECK_EQUAL(test2._result, -1);



    }
} expression_parser_test;
