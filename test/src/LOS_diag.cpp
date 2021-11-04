#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include "LOS.hpp"

#include <boost/test/included/unit_test.hpp>

#define BOOST_TEST_MODULE LOS_test_module

BOOST_AUTO_TEST_SUITE(LOS_DIAG)
BOOST_AUTO_TEST_CASE(MATRIX_10x10)
{
    LOS<double> a("test/file/10");
    a.solve(Conditional::DIAGONAL, false);

    std::vector<double> _Expected(10);
    std::iota (_Expected.begin(), _Expected.end(), 1);
    std::vector<double> _Actual = a.getX();

    for (size_t _Pos = 0; _Pos < _Expected.size(); _Pos++)
        BOOST_CHECK(fabs(_Expected[_Pos] - _Actual[_Pos]) < 0.000001);
}

BOOST_AUTO_TEST_CASE(MATRIX_20x20)
{
    LOS<double> a("test/file/20");
    a.solve(Conditional::DIAGONAL, false);

    std::vector<double> _Expected(20);
    std::iota (_Expected.begin(), _Expected.end(), 1);
    std::vector<double> _Actual = a.getX();

    for (size_t _Pos = 0; _Pos < _Expected.size(); _Pos++)
        BOOST_CHECK(fabs(_Expected[_Pos] - _Actual[_Pos]) < 0.000001);
}

BOOST_AUTO_TEST_CASE(MATRIX_50x50)
{
    LOS<double> a("test/file/50");
    a.solve(Conditional::DIAGONAL, false);

    std::vector<double> _Expected(50);
    std::iota (_Expected.begin(), _Expected.end(), 1);
    std::vector<double> _Actual   = a.getX();

    for (size_t _Pos = 0; _Pos < _Expected.size(); _Pos++)
        BOOST_CHECK(fabs(_Expected[_Pos] - _Actual[_Pos]) < 0.000001);
}
BOOST_AUTO_TEST_SUITE_END()