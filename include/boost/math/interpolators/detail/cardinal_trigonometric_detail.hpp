//  (C) Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_INTERPOLATORS_DETAIL_CARDINAL_TRIGONOMETRIC_HPP
#define BOOST_MATH_INTERPOLATORS_DETAIL_CARDINAL_TRIGONOMETRIC_HPP
#include <cmath>
#include <fftw3.h>
#include <boost/math/constants/constants.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <quadmath.h>
#endif

namespace boost { namespace math { namespace interpolators { namespace detail {

template<typename Real>
class cardinal_trigonometric_detail {
public:
  cardinal_trigonometric_detail(const double* data, size_t length, double t0, double h)
  {
      throw std::domain_error("Not implemented.");
  }
};

template<>
class cardinal_trigonometric_detail<double> {
public:
  cardinal_trigonometric_detail(const double* data, size_t length, double t0, double h) : m_t0{t0}, m_h{h}
  {
    if (length == 0)
    {
        throw std::logic_error("At least one sample is required.");
    }
    if (h <= 0)
    {
        throw std::logic_error("The step size must be > 0");
    }
    // The period sadly must be stored, since the complex vector has length that cannot be used to recover the period:
    m_T = m_h*length;
    m_complex_vector_size = length/2 + 1;
    m_gamma = fftw_alloc_complex(m_complex_vector_size);
    // The const_cast is legitimate: FFTW does not change the data as long as FFTW_ESTIMATE is provided.
    fftw_plan plan = fftw_plan_dft_r2c_1d(length, const_cast<double*>(data), m_gamma, FFTW_ESTIMATE);
    // FFTW says a null plan is impossible with the basic interface we are using, and I have no reason to doubt them.
    // But it just feels weird not to check this:
    if (!plan)
    {
        throw std::logic_error("A null fftw plan was created.");
    }

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    double denom = length;
    for (size_t k = 0; k < m_complex_vector_size; ++k)
    {
        m_gamma[k][0] /= denom;
        m_gamma[k][1] /= denom;
    }
  }

  cardinal_trigonometric_detail(const cardinal_trigonometric_detail& old)  = delete;

  cardinal_trigonometric_detail& operator=(const cardinal_trigonometric_detail&) = delete;

  cardinal_trigonometric_detail(cardinal_trigonometric_detail &&) = delete;

  /*cardinal_trigonometric_detail(cardinal_trigonometric_detail&& other) {
    this->m_t0 = other.m_t0;
    this->m_h = other.m_h;
    this->m_T = other.m_T;
    this->m_gamma = other.m_gamma;
    //other.m_gamma = nullptr;
    this->m_complex_vector_size = m_complex_vector_size;

    std::cout << "Calling move constructor\n";
    std::cout << "Other pointer = " << other.m_gamma << "\n";
    std::cout << "this gamma pointer = " << this->m_gamma << "\n";
    std::cout << "this->m_T = " << this->m_T << "\n";
  }*/

  double operator()(double t) const
  {
      using std::sin;
      using std::cos;
      using boost::math::constants::two_pi;
      using std::exp;
      double s = m_gamma[0][0];
      double x = two_pi<double>()*(t - m_t0)/m_T;
      fftw_complex z;
      z[0] = cos(x);
      z[1] = sin(x);
      fftw_complex b{0, 0};
      // u = b*z
      fftw_complex u;
      for (size_t k = m_complex_vector_size - 1; k >= 1; --k) {
        u[0] = b[0]*z[0] - b[1]*z[1];
        u[1] = b[0]*z[1] + b[1]*z[0];
        b[0] = m_gamma[k][0] + u[0];
        b[1] = m_gamma[k][1] + u[1];
      }

      s += 2*(b[0]*z[0] - b[1]*z[1]);
      return s;
  }

  double period() const {
    return m_T;
  }

  ~cardinal_trigonometric_detail() {
    if (m_gamma) {
      fftw_free(m_gamma);
    }
  }


private:
  double m_t0;
  double m_h;
  double m_T;
  fftw_complex* m_gamma;
  size_t m_complex_vector_size;
};

}}}}
#endif
