/*
Abraham Flores
Notre Dame Physics REU 2016
6/27/2016
Language C++

WaveFunction Class:
    Generates a Wavefunction with the relevant data members
*/

#include <vector>
#include <utility>
#include <algorithm>
#include <math.h>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "wavefunction_class.h"
#include "wavefunction_basis.h"
#include "spline.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;

void BuildLaguerre(MatrixXd &T, int n, double a){
    double on,off;//On or off diagonal

    for(auto i=0;i<n;i++){
        on = 2*i+a+1;
        off = sqrt(i*(i+a));

        T(i,i) = on;
        if(i!=0){
            T(i-1,i) = off;
            T(i,i-1) = off;
        }
    }
}

void BuildJacobi(MatrixXd &T, int n, double a, double b){
    double on,off;//On or off diagonal

    for(auto i=0;i<n;i++){
        on = (b*b-a*a)/((2*i+a+b)*(2*i+a+b+2));
        off = sqrt((4*i*(i+a)*(i+b)*(i+a+b))/(pow(2*i+a+b,2)*(2*i+a+b+1)*(2*i+a+b-1)));

        T(i,i) = on;
        if(i!=0){
            T(i-1,i) = off;
            T(i,i-1) = off;
        }
    }
}

std::vector<double> UniqueEigenValues(MatrixXd M){
    std::vector<double> bounds;
    Eigen::ComplexEigenSolver<MatrixXd> eigen_values(M,false);
    double real;
    std::vector<double>::iterator it;
    for(int i=0;i<M.cols();i++){
        real = eigen_values.eigenvalues()[i].real();
        it = std::find (bounds.begin(), bounds.end(), real);

        if(it == bounds.end()){
            bounds.push_back(real);
        }}

    return bounds;
}

std::vector<double> LaguerreRoots(int n, double a){
    std::vector<double> roots;
    MatrixXd T(n,n);
    T.setZero();
    std::vector<double> eigen_values;

    BuildLaguerre(T,n,a);
    eigen_values = UniqueEigenValues(T);
    for(auto ev: eigen_values){
        if(ev<0){eigen_values.erase(std::remove(eigen_values.begin(), eigen_values.end(), ev), eigen_values.end());}}

    int size_ = eigen_values.size();
    for(int i=0;i<size_;i++){
        roots.push_back(eigen_values[i]);}

    return roots;
}

std::vector<double> LaguerreSquaredRoots(int n, double a){
    std::vector<double> roots;
    MatrixXd T(n,n);
    T.setZero();
    std::vector<double> eigen_values;

    BuildLaguerre(T,n,a);
    eigen_values = UniqueEigenValues(T);

    for(auto ev: eigen_values){
        if(ev<0){eigen_values.erase(std::remove(eigen_values.begin(), eigen_values.end(), ev), eigen_values.end());}}

    for(auto ele: eigen_values){
        if(ele<0){ele*=-1;}}

    int size_ = eigen_values.size();
    for(int i=0;i<size_;i++){
        roots.push_back(eigen_values[i]);
    }

    return roots;
}

std::vector<double> JacobiRoots(int n, double a, double b){
    std::vector<double> roots;
    MatrixXd T(n,n);
    T.setZero();
    std::vector<double> eigen_values;

    BuildJacobi(T,n,a,b);
    eigen_values = UniqueEigenValues(T);

    int size_ = eigen_values.size();
    for(int i=0;i<size_;i++){
        roots.push_back(eigen_values[i]);}

    return roots;
}

std::vector<double> WaveFunction::FindRoots(){
    std::vector<double> roots;

    switch (basis_)
    {
        case Basis::HC:
            roots = LaguerreSquaredRoots(n_,l_+0.5);//Harmonic Oscillator Coordinate Radial Wavefunctions
            std::for_each(begin(roots), end(roots), [&](double& r) { r = sqrt(pow(this->b(),2)*r); });
            break;
        case Basis::HM:
            roots = LaguerreSquaredRoots(n_,l_+0.5);//Harmonic Oscillator Momentum Radial Wavefunctions
            std::for_each(begin(roots), end(roots), [&](double& r) { r = sqrt(r/pow(this->b(),2)); });
            break;
        case Basis::LC:
            roots = LaguerreRoots(n_,2*l_+2);//Laguerre Coordinate Radial Wavefunctions
            std::for_each(begin(roots), end(roots), [&](double& r) { r *= this->b(); });
            break;
        case Basis::LM:
            roots = JacobiRoots(n_,l_+1.5,l_+0.5);//Laguerre Momentum Radial Wavefunctions
            std::for_each(begin(roots), end(roots), [&](double& j) { j = (1.0/this->b())*sqrt((1+j)/(1-j)); });
            break;
    }

    std::sort (roots.begin(), roots.end());

    return roots;

}

std::vector<std::pair<double,double>> WaveFunction::MakeBounds(WaveFunction wf){
    std::vector<std::pair<double,double>> integral_bounds;
    std::vector<double> roots,roots1,roots2;

    roots1 = this->FindRoots();
    roots2 = wf.FindRoots();

    std::vector<double>::iterator it;

    int size_1 = roots1.size();
    int size_2 = roots2.size();

    for(int i=0;i<size_1;i++){roots.push_back(roots1[i]);}

    for(int i=0;i<size_2;i++){
        it = std::find (roots.begin(), roots.end(), roots2[i]);
        if(it == roots2.end()){roots.push_back(roots2[i]);}
    }

    std::sort(roots.begin(), roots.end());

    int size_ = roots.size();

    integral_bounds.push_back(std::make_pair(0,roots[0]));
    for(int i=0;i<size_-1;i++){integral_bounds.push_back(std::make_pair(roots[i],roots[i+1]));}
    integral_bounds.push_back(std::make_pair(roots.back(),pow(10,10)));

    return integral_bounds;
}

double RadialIntegrand(double z ,WaveFunction wf_1, WaveFunction wf_2 ,int order){
//Change of Variable
    double r = z/(1.0-z);
    double jacob = pow(1.0-z,-2);
//********************************************
    double v1,v2;
    switch (wf_1.basis())
    {
        case Basis::HC:
            v1 = basis::HarmonicCoordinate(r,wf_1);//Harmonic Oscillator Coordinate Radial Wavefunctions
            break;
        case Basis::HM:
            v1 = basis::HarmonicMomentum(r,wf_1);//Harmonic Oscillator Momentum Radial Wavefunctions
            break;
        case Basis::LC:
            v1 = basis::LaguerreCoordinate(r,wf_1);//Laguerre Coordinate Radial Wavefunctions
            break;
        case Basis::LM:
            v1 = basis::LaguerreMomentum(r,wf_1);//Laguerre Momentum Radial Wavefunctions
            break;
    }

    switch (wf_2.basis())
    {
        case Basis::HC:
            v2 = basis::HarmonicCoordinate(r,wf_2);//Harmonic Oscillator Coordinate Radial Wavefunctions
            break;
        case Basis::HM:
            v2 = basis::HarmonicMomentum(r,wf_2);//Harmonic Oscillator Momentum Radial Wavefunctions
            break;
        case Basis::LC:
            v2 = basis::LaguerreCoordinate(r,wf_2);//Laguerre Coordinate Radial Wavefunctions
            break;
        case Basis::LM:
            v2 = basis::LaguerreMomentum(r,wf_2);//Laguerre Momentum Radial Wavefunctions
            break;
    }
    return v1*v2*jacob*pow(r,order);
}

void BuildArrays(double *x,double *y, int n,WaveFunction wf_1,WaveFunction wf_2,int order){
    double k = n;
    double point = 0;
    double step = 1.0/k;
    for(int i=0; i<=n; i++,point += step)
    {
        x[i] = point;
        y[i] = RadialIntegrand(point,wf_1,wf_2,order);
    }
}

double WaveFunction::MatrixElement(int n,WaveFunction wf, int order){

    double x[n];
    double y[n];
    BuildArrays(x,y,n,*this,wf,order);

    return spline::CubicIntegrate(x,y,n);

}









