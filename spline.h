/*
  Abraham Flores
  Notre Dame Physics REU 2016
  6/9/2016
  Language C++


  SPLINE INTERPOLATION HEADER
  Restrictions:
  *Take vector or array of data points separate into groups of four
  *Approximate each group with a cubic function
  *Model given function with piecewise cubics
  *Match the first and second derivative of each point.
  *first and last endpoint have their second derivative set to zero

  -These requirements will be met by the gsl functions

  USE:
  give given_func a function of r
  Be able to use the spline to approximate integrals over infinite ranges that undergo a change
  of variables so the they range from zero to one.
  *To a reasonable degree of accuracy

  ASSUMPTIONS:
  *Range from 0 to 1
  *1000 data points //can be adjusted for more or less. Also could be user input

  Resources to gain to better understanding of the math and routines used;
  https://en.wikipedia.org/wiki/Spline_interpolation
  https://www.physics.utah.edu/~detar/phys6720/handouts/cubic_spline/cubic_spline/node1.html
  http://apps.nrbook.com/empanel/index.html#
  - Chapter 3 , specifically 3.3
  https://www.gnu.org/software/gsl/manual/gsl-ref.pdf
  - Chaper 28, Chapter 40

  + 9/28/16 (mac): Add preprocessor flag SPLINE_NO_FANCY_INTEGRATION to
  disable unused Steffen and Akima integration, for use with older
  GSL.
  + 3/26/17 (mac): Reformat code.
  + 03/03/22 (pjf):
    + Fix const-correctness of CubicIntegrate.
    + Add template overload of CubicIntegrate for STL-like containers (with size() and data)
*/

#ifndef CUBIC_SPLINE_INTERPOLATION_
#define CUBIC_SPLINE_INTERPOLATION_

#include <stdexcept>
namespace spline{

  double CubicIntegrate(const double x[], const double y[], int n);
  /*
    x: Array of x axis data points
    y: Array of y axis data points
    n: Number of data points

    Core Function:
    -Initializes types for gsl functions
    -Initializes spline
    -Calls gsl to evaluate integral from lower bound to upper bound
    -returns Integral value
  */

  template<
      typename T,
      decltype(std::declval<T>().size())* = nullptr,
      decltype(std::declval<T>().data())* = nullptr
    >
  inline double CubicIntegrate(const T& x, const T& y)
  {
    if (x.size() != y.size())
      throw std::invalid_argument("unequal length containers");
    return CubicIntegrate(x.data(), y.data(), x.size());
  }

#ifndef SPLINE_NO_FANCY_INTEGRATION
  double AkimaIntegrate(double *x,double *y, int n);
  /*
    x: Array of x axis data points
    y: Array of y axis data points
    n: Number of data points

    Core Function:
    -Initializes types for gsl functions
    -Initializes spline
    -Calls gsl to evaluate integral from lower bound to upper bound
    -returns Integral value
  */

  double SteffenIntegrate(double *x,double *y, int n);
  /*
    x: Array of x axis data points
    y: Array of y axis data points
    n: Number of data points

    Core Function:
    -Initializes types for gsl functions
    -Initializes spline
    -Calls gsl to evaluate integral from lower bound to upper bound
    -returns Integral value
  */
#endif

}
#endif

