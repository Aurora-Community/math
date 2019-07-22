/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <vector>
#include <random>
#include <boost/math/constants/constants.hpp>
#include <boost/math/interpolators/cardinal_trigonometric.hpp>

using std::sin;
using boost::math::constants::two_pi;
using boost::math::interpolators::cardinal_trigonometric;

template<class Real>
void test_constant()
{
    Real t0 = 0;
    Real h = 1;
    for(size_t n = 1; n < 20; ++n)
    {
      Real c = 7.5;
      std::vector<Real> v(2*n+1, c);
      auto ct = cardinal_trigonometric(v, t0, h);
      CHECK_ULP_CLOSE(c, ct(0.3), 3);
    }
}

template<class Real>
void test_sampled_sine()
{
    using std::sin;
    unsigned n = 12;
    Real t0 = 0;
    Real T = 1;
    Real h = T/(2*n+1);

    std::vector<Real> v(2*n+1);
    for(size_t j = 0; j < v.size(); ++j) {
      Real t = t0 + j*h;
      v[j] = sin(two_pi<Real>()*(t-t0)/T);
    }

    auto ct = cardinal_trigonometric(v, t0, h);
}

int main()
{
    test_constant<double>();
    test_sampled_sine<double>();


    return boost::math::test::report_errors();
}
