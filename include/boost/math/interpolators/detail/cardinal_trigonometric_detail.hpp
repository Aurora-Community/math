//  (C) Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_INTERPOLATORS_DETAIL_CARDINAL_TRIGONOMETRIC_HPP
#define BOOST_MATH_INTERPOLATORS_DETAIL_CARDINAL_TRIGONOMETRIC_HPP
#include <cmath>
#include <stdexcept>
#include <fftw3.h>
#include <boost/math/constants/constants.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <quadmath.h>
#endif

namespace boost { namespace math { namespace interpolators { namespace detail {

template<typename Real>
class cardinal_trigonometric_detail {
public:
  cardinal_trigonometric_detail(const Real* data, size_t length, Real t0, Real h)
  {
      m_data = data;
      m_length = length;
      m_t0 = t0;
      m_h = h;
      throw std::domain_error("Not implemented.");
  }
private:
    size_t m_length;
    Real m_t0;
    Real m_h;
    Real* m_data;
};

template<>
class cardinal_trigonometric_detail<float> {
public:
  cardinal_trigonometric_detail(const float* data, size_t length, float t0, float h) : m_t0{t0}, m_h{h}
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
    m_gamma = fftwf_alloc_complex(m_complex_vector_size);
    // The const_cast is legitimate: FFTW does not change the data as long as FFTW_ESTIMATE is provided.
    fftwf_plan plan = fftwf_plan_dft_r2c_1d(length, const_cast<float*>(data), m_gamma, FFTW_ESTIMATE);
    // FFTW says a null plan is impossible with the basic interface we are using, and I have no reason to doubt them.
    // But it just feels weird not to check this:
    if (!plan)
    {
      throw std::logic_error("A null fftw plan was created.");
    }

    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    float denom = length;
    for (size_t k = 0; k < m_complex_vector_size; ++k)
    {
        m_gamma[k][0] /= denom;
        m_gamma[k][1] /= denom;
    }
  }

  cardinal_trigonometric_detail(const cardinal_trigonometric_detail& old)  = delete;

  cardinal_trigonometric_detail& operator=(const cardinal_trigonometric_detail&) = delete;

  cardinal_trigonometric_detail(cardinal_trigonometric_detail &&) = delete;

  float operator()(float t) const
  {
      using std::sin;
      using std::cos;
      using boost::math::constants::two_pi;
      using std::exp;
      float s = m_gamma[0][0];
      float x = two_pi<float>()*(t - m_t0)/m_T;
      fftwf_complex z;
      z[0] = cos(x);
      z[1] = sin(x);
      fftwf_complex b{0, 0};
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

  float period() const
  {
    return m_T;
  }

  float integrate() const {
      return m_T*m_gamma[0][0];
  }

  ~cardinal_trigonometric_detail()
  {
    if (m_gamma)
    {
      fftwf_free(m_gamma);
      m_gamma = nullptr;
    }
  }


private:
  float m_t0;
  float m_h;
  float m_T;
  fftwf_complex* m_gamma;
  size_t m_complex_vector_size;
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
    m_T = m_h*length;
    m_complex_vector_size = length/2 + 1;
    m_gamma = fftw_alloc_complex(m_complex_vector_size);
    fftw_plan plan = fftw_plan_dft_r2c_1d(length, const_cast<double*>(data), m_gamma, FFTW_ESTIMATE);
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
    for (size_t k = m_complex_vector_size - 1; k >= 1; --k)
    {
      u[0] = b[0]*z[0] - b[1]*z[1];
      u[1] = b[0]*z[1] + b[1]*z[0];
      b[0] = m_gamma[k][0] + u[0];
      b[1] = m_gamma[k][1] + u[1];
    }

    s += 2*(b[0]*z[0] - b[1]*z[1]);
    return s;
  }

  double period() const
  {
    return m_T;
  }

  double integrate() const
  {
    return m_T*m_gamma[0][0];
  }

  ~cardinal_trigonometric_detail()
  {
    if (m_gamma)
    {
      fftw_free(m_gamma);
      m_gamma = nullptr;
    }
  }

private:
  double m_t0;
  double m_h;
  double m_T;
  fftw_complex* m_gamma;
  size_t m_complex_vector_size;
};


template<>
class cardinal_trigonometric_detail<long double> {
public:
  cardinal_trigonometric_detail(const long double* data, size_t length, long double t0, long double h) : m_t0{t0}, m_h{h}
  {
    if (length == 0)
    {
      throw std::logic_error("At least one sample is required.");
    }
    if (h <= 0)
    {
      throw std::logic_error("The step size must be > 0");
    }
    m_T = m_h*length;
    m_complex_vector_size = length/2 + 1;
    m_gamma = fftwl_alloc_complex(m_complex_vector_size);
    fftwl_plan plan = fftwl_plan_dft_r2c_1d(length, const_cast<long double*>(data), m_gamma, FFTW_ESTIMATE);
    if (!plan)
    {
      throw std::logic_error("A null fftw plan was created.");
    }

    fftwl_execute(plan);
    fftwl_destroy_plan(plan);

    long double denom = length;
    for (size_t k = 0; k < m_complex_vector_size; ++k)
    {
      m_gamma[k][0] /= denom;
      m_gamma[k][1] /= denom;
    }
  }

  cardinal_trigonometric_detail(const cardinal_trigonometric_detail& old)  = delete;

  cardinal_trigonometric_detail& operator=(const cardinal_trigonometric_detail&) = delete;

  cardinal_trigonometric_detail(cardinal_trigonometric_detail &&) = delete;

  long double operator()(long double t) const
  {
    using std::sin;
    using std::cos;
    using boost::math::constants::two_pi;
    using std::exp;
    long double s = m_gamma[0][0];
    long double x = two_pi<long double>()*(t - m_t0)/m_T;
    fftwl_complex z;
    z[0] = cos(x);
    z[1] = sin(x);
    fftwl_complex b{0, 0};
    fftwl_complex u;
    for (size_t k = m_complex_vector_size - 1; k >= 1; --k)
    {
      u[0] = b[0]*z[0] - b[1]*z[1];
      u[1] = b[0]*z[1] + b[1]*z[0];
      b[0] = m_gamma[k][0] + u[0];
      b[1] = m_gamma[k][1] + u[1];
    }

    s += 2*(b[0]*z[0] - b[1]*z[1]);
    return s;
  }

  long double period() const
  {
    return m_T;
  }

  double integrate() const
  {
    return m_T*m_gamma[0][0];
  }

  ~cardinal_trigonometric_detail()
  {
    if (m_gamma)
    {
      fftwl_free(m_gamma);
      m_gamma = nullptr;
    }
  }

private:
  long double m_t0;
  long double m_h;
  long double m_T;
  fftwl_complex* m_gamma;
  size_t m_complex_vector_size;
};

#ifdef BOOST_HAS_FLOAT128
template<>
class cardinal_trigonometric_detail<__float128> {
public:
  cardinal_trigonometric_detail(const __float128* data, size_t length, __float128 t0, __float128 h) : m_t0{t0}, m_h{h}
  {
    if (length == 0)
    {
      throw std::logic_error("At least one sample is required.");
    }
    if (h <= 0)
    {
      throw std::logic_error("The step size must be > 0");
    }
    m_T = m_h*length;
    m_complex_vector_size = length/2 + 1;
    m_gamma = fftwq_alloc_complex(m_complex_vector_size);
    fftwq_plan plan = fftwq_plan_dft_r2c_1d(length, reinterpret_cast<__float128*>(const_cast<__float128*>(data)), m_gamma, FFTW_ESTIMATE);
    if (!plan)
    {
      throw std::logic_error("A null fftw plan was created.");
    }

    fftwq_execute(plan);
    fftwq_destroy_plan(plan);

    __float128 denom = length;
    for (size_t k = 0; k < m_complex_vector_size; ++k)
    {
      m_gamma[k][0] /= denom;
      m_gamma[k][1] /= denom;
    }
  }

  cardinal_trigonometric_detail(const cardinal_trigonometric_detail& old)  = delete;

  cardinal_trigonometric_detail& operator=(const cardinal_trigonometric_detail&) = delete;

  cardinal_trigonometric_detail(cardinal_trigonometric_detail &&) = delete;

  __float128 operator()(__float128 t) const
  {
    using std::sin;
    using std::cos;
    using boost::math::constants::two_pi;
    using std::exp;
    __float128 s = m_gamma[0][0];
    __float128 x = two_pi<__float128>()*(t - m_t0)/m_T;
    fftwq_complex z;
    z[0] = cosq(x);
    z[1] = sinq(x);
    fftwq_complex b{0, 0};
    fftwq_complex u;
    for (size_t k = m_complex_vector_size - 1; k >= 1; --k)
    {
      u[0] = b[0]*z[0] - b[1]*z[1];
      u[1] = b[0]*z[1] + b[1]*z[0];
      b[0] = m_gamma[k][0] + u[0];
      b[1] = m_gamma[k][1] + u[1];
    }

    s += 2*(b[0]*z[0] - b[1]*z[1]);
    return s;
  }

  __float128 period() const
  {
    return m_T;
  }

  __float128 integrate() const
  {
    return m_T*m_gamma[0][0];
  }


  ~cardinal_trigonometric_detail()
  {
    if (m_gamma)
    {
      fftwq_free(m_gamma);
      m_gamma = nullptr;
    }
  }


private:
  __float128 m_t0;
  __float128 m_h;
  __float128 m_T;
  fftwq_complex* m_gamma;
  size_t m_complex_vector_size;
};
#endif

}}}}
#endif
