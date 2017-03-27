/*
  Abraham Flores
  Notre Dame Physics REU 2016
  6/9/2016
  Language C++

  Algorithm designed to generate C_S basis Wave Functions
  See header for details
*/

#include <stdexcept>
#include <limits>
#include <cmath>
#include <vector>
#include <gsl/gsl_sf_laguerre.h>

#include "wavefunction_basis.h"
#include "wavefunction_class.h"

namespace spline {
  namespace basis{

    double Factorial(int n){
      double fact=1;
      if(n!=0 && n>0){for(int i=n;i>1;i--){fact*=i;}}
      return fact;
    }

    double LaguerrePolynomial(int n, int k, double x){//Breaks down for Large x: DO NOT USE!
      std::vector<double> L = {1.0,1+k-x};
      double next;
      for(int i=1;i<n;i++){
        next = ((2.0*i+1.0+k-x)*L[i]-(i+k)*L[i-1])/(i+1.0);
        // Recursive definition of the Laguerre Polynomials
        L.push_back(next);
      }
      return L[n];
    }

    double JacobiPolynomial(int n, double a, double b, double x){
      std::vector<double> J = {1.0,0.5*(x*(a+b+2)+a-b)};
      double next;
      for(int i=2;i<=n;i++){
        next = ((2.0*i+a+b-1)*((2.0*i+a+b)*(2.0*i+a+b-2)*x+pow(a,2)-pow(b,2))*J[i-1]\
                -2.0*(i+a-1)*(i+b-1)*(2.0*i+a+b)*J[i-2])/(2.0*i*(i+a+b)*(2.0*i+a+b-2));
        // Recursive definition of the Jacobi Polynomials
        J.push_back(next);
      }
      return J[n];
    }

    double HarmonicCoordinate(double r, int n, int l, double b, double power_shift){

      double laguerre = gsl_sf_laguerre_n(n,l+0.5,pow(r/b,2));
      double norm = sqrt(2*Factorial(n)/tgamma(l+n+1.5))*pow(b,-1.5);
      double value = b*pow(1/b,-power_shift)*pow(r/b,l+1.0+power_shift)*exp(-pow(r/b,2)/2.0);

      return value*norm*laguerre;
    }

    double HarmonicMomentum(double k, int n, int l, double b, double power_shift){
      double laguerre = gsl_sf_laguerre_n(n,l+0.5,pow(b*k,2));
      double norm = sqrt(2.0*Factorial(n)/tgamma(l+n+1.5))*pow(b,1.5);
      double value = pow(-1.0,n)/b*pow(b,-power_shift)*pow(k*b,l+1.0+power_shift)*exp(-pow(k*b,2)/2.0);

      return value*norm*laguerre;
    }

    double LaguerreCoordinate(double r, int n, int l, double b, double power_shift){
      double laguerre = gsl_sf_laguerre_n(n,2*l+2,2*r/b);
      double value = b/2*pow(2/b,-power_shift)*pow(2*r/b,l+1+power_shift)*exp(-r/b);
      double norm = pow(2/b,1.5)*sqrt(Factorial(n)/Factorial(n+2*l+2));

      return value*norm*laguerre;
    }

    double LaguerreMomentum(double k, int n, int l, double b, double power_shift){
      double j = (pow(b*k,2)-1)/(pow(b*k,2)+1);
      double jacobi = JacobiPolynomial(n,l+1.5,l+0.5,j);
      double value = pow(b*k,-power_shift)*pow(b*k,l+1+power_shift)/pow(pow(b*k,2)+1,l+2)/b;
      double norm = 2*pow(b,1.5)*sqrt(Factorial(n)*Factorial(n+2*l+2))/tgamma(n+l+1.5);

      return value*norm*jacobi;
    }

    double WaveFunctionValue(double x, WaveFunction wf, double power_shift)
    {
      double value;
      switch (wf.basis())
        {
        case Basis::HC:
          value = basis::HarmonicCoordinate(x,wf.n(),wf.l(),wf.b(),power_shift);//Harmonic Oscillator Coordinate Radial Wavefunctions
          break;
        case Basis::HM:
          value = basis::HarmonicMomentum(x,wf.n(),wf.l(),wf.b(),power_shift);//Harmonic Oscillator Momentum Radial Wavefunctions
          break;
        case Basis::LC:
          value = basis::LaguerreCoordinate(x,wf.n(),wf.l(),wf.b(),power_shift);//Laguerre Coordinate Radial Wavefunctions
          break;
        case Basis::LM:
          value = basis::LaguerreMomentum(x,wf.n(),wf.l(),wf.b(),power_shift);//Laguerre Momentum Radial Wavefunctions
          break;
        }
      return value;
    }

  }//close namespace
}  // namespace spline
