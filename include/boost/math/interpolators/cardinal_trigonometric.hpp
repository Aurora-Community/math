//  (C) Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_INTERPOLATORS_CARDINAL_TRIGONOMETRIC_HPP
#define BOOST_MATH_INTERPOLATORS_CARDINAL_TRIGONOMETRIC_HPP
#include <cmath>
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
        if (v.size() == 0)
        {
            throw std::logic_error("At least one sample is required.");
        }
        // The period sadly must be stored, since the complex vector has length that cannot be used to recover the period:
        m_T = m_h*v.size();
        if constexpr(std::is_same_v<Real, double>)
        {
            complex_vector_size = v.size()/2 + 1;
            gamma = fftw_alloc_complex(complex_vector_size);
            // The const_cast is legitimate: FFTW does not change the data as long as FFTW_ESTIMATE is provided.
            fftw_plan plan = fftw_plan_dft_r2c_1d(v.size(), const_cast<double*>(v.data()), gamma, FFTW_ESTIMATE);
            // FFTW says a null plan is impossible with the basic interface we are using, and I have no reason to doubt them.
            // But it just feels weird not to check this:
            if (!plan)
            {
                throw std::logic_error("A null fftw plan was created.");
            }

            fftw_execute(plan);
            fftw_destroy_plan(plan);

            Real denom = v.size();
            for (size_t k = 0; k < complex_vector_size; ++k)
            {
                gamma[k][0] /= denom;
                gamma[k][1] /= denom;
            }
        }
        else
        {
          throw std::logic_error("Cardinal trigonometric interpolation not yet implemented on this type.");
        }
    }

    Real operator()(Real t) const
    {
        using std::sin;
        using std::cos;
        using boost::math::constants::two_pi;
        using std::exp;
        Real s = gamma[0][0];
        Real x = two_pi<Real>()*(t - m_t0)/m_T;
        fftw_complex z;
        z[0] = cos(x);
        z[1] = sin(x);
        fftw_complex b{0, 0};
        // u = b*z
        fftw_complex u;
        for (size_t k = complex_vector_size - 1; k >= 1; --k) {
          u[0] = b[0]*z[0] - b[1]*z[1];
          u[1] = b[0]*z[1] + b[1]*z[0];
          b[0] = gamma[k][0] + u[0];
          b[1] = gamma[k][1] + u[1];
        }

        s += 2*(b[0]*z[0] - b[1]*z[1]);
        return s;
    }

    Real period() const {
      return m_T;
    }

    ~cardinal_trigonometric() {
      fftw_free(gamma);
    }

private:
    Real m_t0;
    Real m_h;
    Real m_T;
    fftw_complex* gamma;
    size_t complex_vector_size;
};

}}}
#endif
