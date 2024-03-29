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
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "wavefunction_class.h"
#include "wavefunction_basis.h"
#include "spline.h"

namespace spline {
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

  double RadialIntegrand(double z, WaveFunction wf_1, WaveFunction wf_2, int order){

    // change of variable
    double r = z/(1.0-z);
    double jacob = pow(1.0-z,-2);

    // evaluate integrand
    double value;
    if (order >= 0)
      {
        double value1 = basis::WaveFunctionValue(r,wf_1);
        double value2 = basis::WaveFunctionValue(r,wf_2);
        value = value1*value2*jacob*pow(r,order);
      }
    else if (order == -1)
      // order negative
      //
      // factor an r out of each wave function to avoid divide by zero
      // if order -1 or -2 (operators r^-1 and r^-2)
      {
        double value1 = basis::WaveFunctionValue(r,wf_1,-0.5);
        double value2 = basis::WaveFunctionValue(r,wf_2,-0.5);
        value = value1*value2*jacob;
      }
    else if (order == -2)
      {
        double value1 = basis::WaveFunctionValue(r,wf_1,-1);
        double value2 = basis::WaveFunctionValue(r,wf_2,-1);
        value = value1*value2*jacob;
      }
    return value;

  }

  void BuildArrays(double x[], double y[], int num_size, WaveFunction wf_1, WaveFunction wf_2, int order){
    double num_steps = (num_size - 1);
    double interval_width = 1.;
    double step = interval_width/num_steps;

    for (int i=0; i < (num_size-1); ++i)
      {
        x[i] = i*step;
        y[i] = RadialIntegrand(x[i], wf_1, wf_2, order);
      }

    // handle last point separately (to avoid divide by zero in conversion from z -> r)
    x[num_size-1] = 1.0;
    y[num_size-1] = 0.0;
  }

  double WaveFunction::MatrixElement(int num_size, WaveFunction wf, int order){

    // (mac): "n" is apparently the number of *points*, not *steps*, following gsl_interp
    // conventions; rename to num_size (as in BuildArrays)

    double x[num_size];
    double y[num_size];
    BuildArrays(x,y,num_size,*this,wf,order);

    return CubicIntegrate(x,y,num_size);

  }


  double /*WaveFunction::*/MatrixElement(int num_size, WaveFunction wf1, WaveFunction wf2, int order){

    // (mac): "n" is apparently the number of *points*, not *steps*, following gsl_interp
    // conventions; rename to num_size (as in BuildArrays)

    double x[num_size];
    double y[num_size];
    BuildArrays(x,y,num_size,wf1/**this*/,wf2,order);

    return CubicIntegrate(x,y,num_size);

  }

}  // namespace spline
