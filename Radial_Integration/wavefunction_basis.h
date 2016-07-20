/*
Abraham Flores
Notre Dame Physics REU 2016
6/9/2016
Language C++

CS_BASIS HEADER
*/

#ifndef WAVEFUNCTION_BASIS_FUNCTIONS
#define WAVEFUNCTION_BASIS_FUNCTIONS

#include "wavefunction_class.h"
namespace basis{

double Factorial(int n);
//Calulates the factiorial of n = n!

double LaguerrePolynomial(int n, int k, double x);
//More information https://en.wikipedia.org/wiki/Laguerre_polynomials
//Returns the associated Laguerre polynomial of order n evaluated at x.
//n: polynomial order, must be an integer >= 0
//x: point at which to evaluate the function

double JacobiPolynomial(int n, double a, double b, double x);
//More information https://en.wikipedia.org/wiki/Jacobi_polynomials
//Returns the Hermite polynomial of order n evaluated at x.
//n: polynomial order, must be an integer >= 0
//x: point at which to evaluate the function

double HarmonicCoordinate(double r, int n, int l, double b);
//Evaluates Radial wavefunctions in coordinate space at point r for the Harmonic Oscillator basis
//Defined in Coulomb-Sturmian(Laguerre functions) basis for the nuclear many-body problem. M.A.CAPRIO,P.MARIS AND J.P.VARY
// Equation 4
/*
r: point at which to evaluate the function
b: length parameter
n: quantum number n
l: quantum number l
*/
double HarmonicCoordinate(double r, WaveFunction wf);

double HarmonicMomentum(double k, int n, int l, double b);
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
double HarmonicMomentum(double k, WaveFunction wf);

double LaguerreCoordinate(double r, int n, int l, double b);
//Evaluates Radial wavefunctions in coordinate space at point r for the cs basis
//Defined in Coulomb-Sturmian(Laguerre functions) basis for the nuclear many-body problem. M.A.CAPRIO,P.MARIS AND J.P.VARY
// Equation 34
/*
r: point at which to evaluate the function
b: length parameter
n: quatnum number n
l: quantum number l
*/

double LaguerreCoordinate(double r, WaveFunction wf);

double LaguerreMomentum(double k, int n, int l, double b);
//Evaluates Radial wavefunctions in momentum eigenspace at momentum k for the cs basis
//Defined in Coulomb-Sturmian(Laguerre functions) basis for the nuclear many-body problem. M.A.CAPRIO,P.MARIS AND J.P.VARY
// Equation 35
/*
k: momentum at which to evaluate the function
b: length parameter
n: energy level, integer >= 0
l: quantum number l
*/
double LaguerreMomentum(double k, WaveFunction wf);

}//Close Namespace
#endif

