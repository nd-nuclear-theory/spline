/*
  WaveFunction Class:
  Generates a Wavefunction with the relevant data members

  Abraham Flores
  Notre Dame Physics REU 2016
  6/27/2016
  Language C++

  + 3/26/17 (mac):
    - Reformat code.
    - Make use of updated dispatch function (WaveFunctionValue)
      interface to wavefunction_basis.
    - Use shift in wave function monomial power to evaluate matrix elements
      of 1/r without divide-by-zero.
*/

#ifndef RADIAL_WAVEFUNCTIONS_
#define RADIAL_WAVEFUNCTIONS_

#include <iostream>
#include <vector>
#include <utility>
#include <cmath>

namespace spline {

  enum class Basis {HC=0,HM=1,LC=2,LM=3};

  class WaveFunction{
    private:
    //Data Members
    int n_;
    int l_;
    double b_;
    Basis basis_;
    public:
    //Member Functions
    WaveFunction() : n_(1), l_(0), b_(1.0), basis_(Basis::HC) {};
    WaveFunction(int n, int l, double b, Basis basis) : n_(n), l_(l), b_(b), basis_(basis) {};
    WaveFunction(int n, int l, Basis basis) : n_(n), l_(l), b_(sqrt(2/(l+3))), basis_(basis) {};
    WaveFunction(const WaveFunction& wf): n_(wf.n_), l_(wf.l_), b_(wf.b_), basis_(wf.basis_) {};
    ~WaveFunction(){};
    int n(){return n_;};
    int l(){return l_;};
    double b(){return b_;};
    Basis basis(){return basis_;};
    std::vector<double> FindRoots();
    std::vector<std::pair<double,double>> MakeBounds(WaveFunction);
    double MatrixElement(int num_size, WaveFunction wf, int order);

    friend std::ostream& operator<<(std::ostream& out, WaveFunction& wf){
      out<<"Quantum Number n: "<<wf.n_;
      out<<" Quantum Number l: "<<wf.l_;
      out<<" Length Parameter : "<<wf.b_;
      switch (wf.basis_){
      case Basis::HC:
        out<<" Basis: Harmonic Oscillator In Coordinate Space"<<"\n";
        break;
      case Basis::HM:
        out<<" Basis: Harmonic Oscillator In Momentum Space"<<"\n";
        break;
      case Basis::LC:
        out<<" Basis: Laguerre In Coordinate Space"<<"\n";
        break;
      case Basis::LM:
        out<<" Basis: Laguerre In Momentum Space"<<"\n";
        break;

      }
      return out;
    };
  };

	// V. Constantinou (04/06/17)

	// Add in header file wavefunction_class.h so RadialIntegrand can be used outside wavefunction_class.cpp

	double RadialIntegrand(double z, WaveFunction wf_1, WaveFunction wf_2, int order);

	double MatrixElement(int num_size, WaveFunction wf1, WaveFunction wf2, int order);


}  // namespace spline
#endif
