/*
Abraham Flores
Notre Dame Physics REU 2016
6/27/2016
Language C++

WaveFunction Class:
    Generates a Wavefunction with the relevant data members

*/
#ifndef RADIAL_WAVEFUNCTIONS
#define RADIAL_WAVEFUNCTIONS

#include <iostream>
#include <vector>
#include <utility>
#include <math.h>

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
        double MatrixElement(int,WaveFunction,int);

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

#endif

