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
#include <math.h>
#include <vector>
#include <gsl/gsl_sf_laguerre.h>

#include "wavefunction_basis.h"
#include "wavefunction_class.h"

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

double HarmonicCoordinate(double r, int n, int l, double b){

    double laguerre = gsl_sf_laguerre_n(n,l+0.5,pow(r/b,2));
    double norm = sqrt(2*Factorial(n)/tgamma(l+n+1.5))*pow(b,-1.5);
    double value = b*pow(r/b,l+1.0)*exp(-pow(r/b,2)/2.0);

    return value*norm*laguerre;
}

double HarmonicCoordinate(double r, WaveFunction wf){

    if(wf.basis() != Basis::HC){throw std::invalid_argument("INVALID BASIS");}

    double laguerre = gsl_sf_laguerre_n(wf.n(),wf.l()+0.5,pow(r/wf.b(),2));
    double norm = sqrt(2*Factorial(wf.n())/tgamma(wf.l()+wf.n()+1.5))*pow(wf.b(),-1.5);
    double value = wf.b()*pow(r/wf.b(),wf.l()+1.0)*exp(-pow(r/wf.b(),2)/2.0);

    return value*norm*laguerre;
}

double HarmonicMomentum(double k, int n, int l, double b){
    double laguerre = gsl_sf_laguerre_n(n,l+0.5,pow(b*k,2));
    double norm = sqrt(2.0*Factorial(n)/tgamma(l+n+1.5))*pow(b,1.5);
    double value = pow(-1.0,n)/b*pow(k*b,l+1.0)*exp(-pow(k*b,2)/2.0);

    return value*norm*laguerre;
}

double HarmonicMomentum(double k, WaveFunction wf){

    if(wf.basis() != Basis::HM){throw std::invalid_argument("INVALID BASIS");}

    double laguerre = gsl_sf_laguerre_n(wf.n(),wf.l()+0.5,pow(wf.b()*k,2));
    double norm = sqrt(2.0*Factorial(wf.n())/tgamma(wf.l()+wf.n()+1.5))*pow(wf.b(),1.5);
    double value = pow(-1.0,wf.n())/wf.b()*pow(k*wf.b(),wf.l()+1.0)*exp(-pow(k*wf.b(),2)/2.0);

    return value*norm*laguerre;
}

double LaguerreCoordinate(double r, int n, int l, double b){
    double laguerre = gsl_sf_laguerre_n(n,2*l+2,2*r/b);
    double value = b/2*pow(2*r/b,l+1)*exp(-r/b);
    double norm = pow(2/b,1.5)*sqrt(Factorial(n)/Factorial(n+2*l+2));

    return value*norm*laguerre;
}

double LaguerreCoordinate(double r, WaveFunction wf){

    if(wf.basis() != Basis::LC){throw std::invalid_argument("INVALID BASIS");}

    double laguerre = gsl_sf_laguerre_n(wf.n(),2*wf.l()+2,2*r/wf.b());
    double value = wf.b()/2*pow(2*r/wf.b(),wf.l()+1)*exp(-r/wf.b());
    double norm = pow(2/wf.b(),1.5)*sqrt(Factorial(wf.n())/Factorial(wf.n()+2*wf.l()+2));

    return value*norm*laguerre;
}

double LaguerreMomentum(double k, int n, int l, double b){
    double j = (pow(b*k,2)-1)/(pow(b*k,2)+1);
    double jacobi = JacobiPolynomial(n,l+1.5,l+0.5,j);
    double value = pow(b*k,l+1)/pow(pow(b*k,2)+1,l+2)/b;
    double norm = 2*pow(b,1.5)*sqrt(Factorial(n)*Factorial(n+2*l+2))/tgamma(n+l+1.5);

    return value*norm*jacobi;
}

double LaguerreMomentum(double k, WaveFunction wf){

    if(wf.basis() != Basis::LM){throw std::invalid_argument("INVALID BASIS");}

    double j = (pow(wf.b()*k,2)-1)/(pow(wf.b()*k,2)+1);
    double jacobi = JacobiPolynomial(wf.n(),wf.l()+1.5,wf.l()+0.5,j);
    double value = pow(wf.b()*k,wf.l()+1)/pow(pow(wf.b()*k,2)+1,wf.l()+2)/wf.b();
    double norm = 2*pow(wf.b(),1.5)*sqrt(Factorial(wf.n())*Factorial(wf.n()+2*wf.l()+2))/tgamma(wf.n()+wf.l()+1.5);

    return value*norm*jacobi;
}

}//close namespace
