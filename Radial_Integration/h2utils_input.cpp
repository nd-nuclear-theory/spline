/*
Abraham Flores
Notre Dame Physics REU 2016
6/9/2016
Language C++

*/
#include <iostream>
#include <string>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <math.h>
#include <vector>
#include <algorithm>
#include <utility>

#include "spline.h"
#include "wavefunction_class.h"


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;

void WriteFile(std::string file,std::vector<int> parm, std::vector<MatrixXd> data){

    std::ofstream myfile;
    myfile.open(file);

    for(auto ele:parm){myfile<<ele<<" ";}
    myfile<<"\n";
    myfile.precision(17);
    for(auto M:data){myfile<<M<<"\n";}
    myfile.close();
}

void RadialMatrixElement(int delta_l){

    int n = 1000;

    std::vector<std::string> hor = {"radial-me-HO-r0.dat","radial-me-HO-r1.dat","radial-me-HO-r2.dat"};
    std::vector<std::string> csr = {"radial-me-CSFN-r0.dat","radial-me-CSFN-r1.dat","radial-me-CSFN-r2.dat"};
    std::vector<std::string> hok = {"radial-me-HO-k1.dat","radial-me-HO-k2.dat"};
    std::vector<std::string> csk = {"radial-me-CSFN-k1.dat","radial-me-CSFN-k2.dat"};
    std::vector<std::vector<std::string>> all_files = {hor,hok,csr,csk};

    std::vector<MatrixXd> data;
    WaveFunction wf_1;
    WaveFunction wf_2;
    Basis basis1 = Basis::HC;
    Basis basis2 = Basis::HC;

    int l_max = 20;
    int M = 11;
    int M_ = M;
    MatrixXd T(11,11);

    int count_ = 0;
    int order = 0;

    for(auto files:all_files){
        for(auto f:files){
            order = 0;
            if(count_==1){
                order = 1;
                basis1 = Basis::HM;
                basis2 = Basis::HM;
            }
            else if(count_==2){
                basis1 = Basis::HM;
                basis2 = Basis::LC;
            }
            else if(count_==3){
                order = 1;
                basis1 = Basis::LM;
                basis2 = Basis::LM;
            }
            data.clear();
            for(int l_a = 0;l_a<=l_max;l_a++){
                for(int l_b = l_a+(delta_l%2);l_b<=std::min(l_a+delta_l,l_max);l_b++){
                    //New Blocks
                    if(l_a !=0 ){data.push_back(T);}
                    T.setZero();
                    for(int n_a = 0;n_a<M;n_a++){
                        for(int n_b = 0;n_b<M_;n_b++){
                            wf_1 = WaveFunction(n_a,l_a,1.0,basis1);
                            wf_2 = WaveFunction(n_b,l_b,1.0,basis2);
                            T(n_a,n_b) = wf_1.MatrixElement(n,wf_2,order);
                        }
                    }
                }
            }
            data.push_back(T);
            std::vector<int> parm = {l_max,delta_l,M,M_};
            WriteFile(f,parm,data);

            order++;
        }
        count_++;
    }

}

void XForm(){

    int n=100;


    int order = 0;
    std::vector<std::pair<std::string,double>> files =
    {std::make_pair("radial-xform-HO-2.000.dat", 2.000),std::make_pair("radial-xform-HO-1.915.dat", 1.915),std::make_pair("radial-xform-HO-1.834.dat",1.834),std::make_pair("radial-xform-HO-1.756.dat",1.756),\
    std::make_pair("radial-xform-HO-1.682.dat",1.682),std::make_pair("radial-xform-HO-1.610.dat",1.610),std::make_pair("radial-xform-HO-1.542.dat",1.542),std::make_pair("radial-xform-HO-1.477.dat",1.477),\
    std::make_pair("radial-xform-HO-1.414.dat",1.414),std::make_pair("radial-xform-HO-1.354.dat",1.354),std::make_pair("radial-xform-HO-1.297.dat",1.297),std::make_pair("radial-xform-HO-1.242.dat",1.242),\
    std::make_pair("radial-xform-HO-1.189.dat",1.189),std::make_pair("radial-xform-HO-1.139.dat",1.139),std::make_pair("radial-xform-HO-1.091.dat",1.091),std::make_pair("radial-xform-HO-1.044.dat",1.044),\
    std::make_pair("radial-xform-HO-1.000.dat",1.000)};

    std::vector<MatrixXd> data;
    WaveFunction wf_1;
    WaveFunction wf_2;
    Basis basis1 = Basis::HC;
    Basis basis2 = Basis::HC;

    int l_max = 20;
    int delta_l = 0;
    int M = 11;
    int M_ = M;
    MatrixXd T(11,11);

    for(auto ele:files){
            data.clear();
            for(int l_a = 0;l_a<=l_max;l_a++){
                for(int l_b = l_a+(delta_l%2);l_b<=std::min(l_a+delta_l,l_max);l_b++){
                    //New Blocks
                    if(l_a !=0 ){data.push_back(T);}
                    T.setZero();
                    for(int n_a = 0;n_a<M;n_a++){
                        for(int n_b = 0;n_b<M_;n_b++){
                            wf_1 = WaveFunction(n_a,l_a,1.0,basis1);
                            wf_2 = WaveFunction(n_b,l_b,ele.second,basis2);
                            T(n_a,n_b) = wf_1.MatrixElement(n,wf_2,order);
                        }
                    }
                }
            }
        data.push_back(T);
        std::vector<int> parm = {l_max,delta_l,M,M_};
        WriteFile(ele.first,parm,data);
    }

}

int main ()
{
  //RadialMatrixElement(0);
  //XForm();

 WaveFunction wf = WaveFunction(5,4,1.0,Basis::HC);
 for(auto ele:wf.FindRoots()){std::cout<<ele<<std::endl;}
  return 0;
}

























