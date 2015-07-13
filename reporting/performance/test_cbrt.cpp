//  Copyright John Maddock 2015.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#  pragma warning (disable : 4224)
#endif

#include <boost/math/special_functions/cbrt.hpp>
#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include "../../test/table_type.hpp"
#include "table_helper.hpp"
#include "performance.hpp"
#include <iostream>

int main()
{
   typedef double T;
#define SC_(x) static_cast<double>(x)
#  include "../../test/cbrt_data.ipp"

   add_data(cbrt_data);

   unsigned data_total = data.size();

   screen_data([](const std::vector<double>& v){  return boost::math::cbrt(v[1]);  }, [](const std::vector<double>& v){ return v[0];  });

#if defined(TEST_C99) && !defined(COMPILER_COMPARISON_TABLES)
   screen_data([](const std::vector<double>& v){  return ::cbrt(v[1]);  }, [](const std::vector<double>& v){ return v[0];  });
#endif

   unsigned data_used = data.size();
   std::string function = "cbrt[br](" + boost::lexical_cast<std::string>(data_used) + "/" + boost::lexical_cast<std::string>(data_total) + " tests selected)";

   double time = exec_timed_test([](const std::vector<double>& v){  return boost::math::cbrt(v[1]);  });
   std::cout << time << std::endl;
#if defined(COMPILER_COMPARISON_TABLES)
   report_execution_time(time, std::string("Compiler Option Comparison on ") + BOOST_PLATFORM, "boost::math::cbrt", get_compiler_options_name());
#else
   report_execution_time(time, std::string("Library Comparison on ") + BOOST_PLATFORM, function, "Boost");
#endif
   //
   // Boost again, but with promotion to long double turned off:
   //
#if !defined(COMPILER_COMPARISON_TABLES)
   if(sizeof(long double) != sizeof(double))
   {
      double time = exec_timed_test([](const std::vector<double>& v){  return boost::math::cbrt(v[1], boost::math::policies::make_policy(boost::math::policies::promote_double<false>()));  });
      std::cout << time << std::endl;
      report_execution_time(time, std::string("Library Comparison on ") + BOOST_PLATFORM, function, "Boost[br](no internal promotion to long double)");
   }
#endif


#if defined(TEST_C99) && !defined(COMPILER_COMPARISON_TABLES)
   time = exec_timed_test([](const std::vector<double>& v){  return ::cbrt(v[1]);  });
   std::cout << time << std::endl;
   report_execution_time(time, std::string("Library Comparison on ") + BOOST_PLATFORM, function, "math.h");
#endif


   return 0;
}

