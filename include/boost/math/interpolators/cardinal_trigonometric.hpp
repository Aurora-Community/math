//  (C) Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_INTERPOLATORS_CARDINAL_TRIGONOMETRIC_HPP
#define BOOST_MATH_INTERPOLATORS_CARDINAL_TRIGONOMETRIC_HPP
#include <cmath>
#include <type_traits>
#include <fftw3.h>
#include <boost/math/constants/constants.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <quadmath.h>
#endif

namespace boost { namespace math { namespace interpolators {

template<class RandomAccessContainer>
class cardinal_trigonometric
{
public:
    using Real = typename RandomAccessContainer::value_type;
    cardinal_trigonometric(RandomAccessContainer const & v, Real t0, Real h) : m_t0{t0}, m_h{h}
    {
        if (v.size() % 2 == 0)
        {
            throw std::logic_error("Even length has not yet been implemented.");
        }
        if constexpr(std::is_same_v<Real, double>)
        {

            size_t complex_vector_size = v.size()/2 + 1;
            fftw_complex* gamma = fftw_alloc_complex(v.size()/2 + 1);
            fftw_plan plan = fftw_plan_dft_r2c_1d(v.size(), const_cast<double*>(v.data()), gamma, FFTW_ESTIMATE);

            fftw_execute(plan);
            fftw_destroy_plan(plan);

            for (size_t k = 0; k < complex_vector_size; ++k) {
              std::cout << "gamma_" << k << " = {" << gamma[k][0] << ", " << gamma[k][1] << "}\n";
            }

            m_a.resize(complex_vector_size, std::numeric_limits<Real>::quiet_NaN());
            m_b.resize(complex_vector_size, std::numeric_limits<Real>::quiet_NaN());
            for (size_t k = 0; k < complex_vector_size; ++k) {
              m_a[k] = 2*gamma[k][0]/v.size();
            }

            for (size_t k = 0; k < complex_vector_size; ++k) {
              std::cout << "a_" << k << " = " << m_a[k] << "\n";
            }
            std::cout << "a0/2 = " << m_a[0]/2 << "\n";

            m_b[0] = std::numeric_limits<Real>::quiet_NaN();
            for (size_t k = 1; k < complex_vector_size; ++k) {
              m_b[k] = -2*gamma[k][1]/v.size();
            }

            for (size_t k = 0; k < m_b.size(); ++k) {
              std::cout << "b_" << k << " = " << m_b[k] << "\n";
            }


            fftw_free(gamma);
        }

    }

    Real operator()(Real t) const
    {
        using std::sin;
        using std::cos;
        using boost::math::constants::two_pi;
        Real s =  m_a[0]/2;
        Real x = two_pi<Real>()*(t - m_t0)/(m_h*(2*m_a.size()+1));
        for (size_t k = 1; k < m_a.size(); ++k) {
          s += m_a[k]*cos(k*x);
          s += m_b[k]*sin(k*x);
        }
        return s;
    }

private:
    Real m_h;
    Real m_t0;
    RandomAccessContainer m_a;
    RandomAccessContainer m_b;
};

}}}
#endif
