/*
  Abraham Flores
  Notre Dame Physics REU 2016
  6/9/2016
  Language C++

  + 3/26/17 (mac):
    - Reformat code.
    - Add dispatch function WaveFunctionValue for evaluating wave functions from any basis.
    - Eliminate cut-and-paste duplicate (r,n,l,b) and (r,wf) versions of each radial function.
    - Allow for shift in power of r multiplying wave function to facilitate avoiding
      divide-by-zero on operators with small negative powers of r.
  + 10/17/20 (pjf):
    - Prevent overflows by evaluating harmonic oscillator norms with sum of
      log(gamma()) instead of product of gamma().
    - Remove Factorial(), LaguerrePolynomial(), and JacobiPolynomial() from header.
*/

#ifndef WAVEFUNCTION_BASIS_FUNCTIONS_
#define WAVEFUNCTION_BASIS_FUNCTIONS_

#include "wavefunction_class.h"
namespace spline {
  namespace basis{

    double HarmonicCoordinate(double r, int n, int l, double b, double power_shift=0);
    //Evaluates Radial wavefunctions in coordinate space at point r for the Harmonic Oscillator basis
    //Defined in Coulomb-Sturmian(Laguerre functions) basis for the nuclear many-body problem. M.A.CAPRIO,P.MARIS AND J.P.VARY
    // Equation 4
    /*
      r: point at which to evaluate the function
      b: length parameter
      n: quantum number n
      l: quantum number l
    */

    double HarmonicMomentum(double k, int n, int l, double b, double power_shift=0);
    //Evaluates Radial wavefunctions in momentum eigenspace at momentum k for the Harmonic Oscillator basis
    //Defined in Coulomb-Sturmian(Laguerre functions) basis for the nuclear many-body problem. M.A.CAPRIO,P.MARIS AND J.P.VARY
    // Equation 8
    /*
      k: momentum at which to evaluate the function
      b: length parameter
      n: quantum number n
      l: quantum number l
      hbar = 1
    */

    double LaguerreCoordinate(double r, int n, int l, double b, double power_shift=0);
    //Evaluates Radial wavefunctions in coordinate space at point r for the cs basis
    //Defined in Coulomb-Sturmian(Laguerre functions) basis for the nuclear many-body problem. M.A.CAPRIO,P.MARIS AND J.P.VARY
    // Equation 34
    /*
      r: point at which to evaluate the function
      b: length parameter
      n: quatnum number n
      l: quantum number l
    */

    double LaguerreMomentum(double k, int n, int l, double b, double power_shift=0);
    //Evaluates Radial wavefunctions in momentum eigenspace at momentum k for the cs basis
    //Defined in Coulomb-Sturmian(Laguerre functions) basis for the nuclear many-body problem. M.A.CAPRIO,P.MARIS AND J.P.VARY
    // Equation 35
    /*
      k: momentum at which to evaluate the function
      b: length parameter
      n: energy level, integer >= 0
      l: quantum number l
    */


    double WaveFunctionValue(double x, WaveFunction wf, double power_shift=0);
    // Calculate value of wave function according to stored
    // information as to what basis it is a part of.
    //
    // Serves as dispatch function to above specialized wave
    // functions.
    //
    // Note: When the power is shifted, the dimensional factor is
    // compensated, so only the power of x itself changes.
    //
    // Arguments:
    //   x: coordinate or momentum
    //   wf: wave function labels
    //   power_shift: exponent on radial monomial factor exponent



  }//Close Namespace
}  // namespace spline
#endif
