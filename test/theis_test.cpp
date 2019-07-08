/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <utility>
#include <random>
#include <boost/core/demangle.hpp>
#include <boost/math/interpolators/theis.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

using boost::math::interpolators::theis;


template<class Real>
void test_knots()
{
    Real t0 = 0;
    Real h = Real(1)/Real(16);
    size_t n = 512;
    std::vector<Real> v(n);
    std::mt19937 gen(323723);
    std::uniform_real_distribution<Real> dis(1.0, 2.0);

    for(size_t i = 0;  i < n; ++i) {
      v[i] = static_cast<Real>(dis(gen));
    }
    auto s = theis<decltype(v)>(std::move(v), t0, h);

    size_t i = 0;
    while (i < n) {
      Real t = t0 + i*h;
      Real expected = s[i];
      Real computed = s(t);
      CHECK_ULP_CLOSE(expected, computed, 16);
      ++i;
    }
}

template<class Real>
void test_bump()
{
    using std::exp;
    using std::abs;
    using std::sqrt;
    auto bump = [](Real x) { if (abs(x) >= 1) { return Real(0); } return exp(-Real(1)/(Real(1)-x*x)); };

    Real t0 = -1;
    size_t n = 8192+1;
    Real h = Real(2)/Real(n-1);

    std::vector<Real> v(n);
    for(size_t i = 0; i < n; ++i) {
        Real t = t0 + i*h;
        v[i] = bump(t);
    }


    std::vector<Real> v_copy = v;
    auto s = theis<decltype(v)>(std::move(v), t0, h);

    // Test the knots:
    for(size_t i = v_copy.size()/4; i < 3*v_copy.size()/4; ++i) {
        Real t = t0 + i*h;
        Real expected = v_copy[i];
        Real computed = s(t);
        if(!CHECK_MOLLIFIED_CLOSE(expected, computed, 10*std::numeric_limits<Real>::epsilon())) {
            std::cerr << "  Problem occured at abscissa " << t << "\n";
        }
    }

    std::mt19937 gen(323723);
    std::uniform_real_distribution<long double> dis(-0.85, 0.85);

    size_t i = 0;
    while (i++ < 1000)
    {
        Real t = static_cast<Real>(dis(gen));
        Real expected = bump(t);
        Real computed = s(t);
        if(!CHECK_MOLLIFIED_CLOSE(expected, computed, 10*std::numeric_limits<Real>::epsilon())) {
            std::cerr << "  Problem occured at abscissa " << t << "\n";
        }
    }
}


int main()
{
    test_knots<float>();
    test_knots<double>();
    test_knots<long double>();

    test_bump<double>();
    test_bump<long double>();

    return boost::math::test::report_errors();
}
