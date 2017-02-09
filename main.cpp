
//for give paramter, store the magnetization into file
//and store the spin configuration into file2

#include <iostream>
#include <fstream>
#include <stdlib.h>     /* atof atoi */


#include "next.h"
#include "stat.h"

int main(){

  const int Nx=50;
  const int Ntau=200;
  double K=.13;
  double Gamma=0.136;
  double alpha=0.6;
  const int cycle=5000;
  const int TJcycle=20;

  int* spin = new int[Ntau*Nx];
  double* corr = new double[Ntau*Nx];
  double* corrSum = new double[Ntau*Nx];
  int cSum=0;

  int currentCycle=initialize_2d(spin,Nx,Ntau,"Zspin.txt");

  ofstream finfo;
  finfo.open("Zinfo.txt");
  ofstream fmag;
  fmag.open("Zmag.txt",ios::app);
  ofstream fcor;
  fcor.open("Zcor.txt");//,ios::app
  ofstream fspin;
  fspin.open("Zspin.txt",ios::trunc);

  finfo<<"Nx"<<'\t'<<Nx<<'\n';
  finfo<<"Ntau"<<'\t'<<Ntau<<'\n';
  finfo<<"K"<<'\t'<<K<<'\n';
  finfo<<"Gamma"<<'\t'<<Gamma<<'\n';
  finfo<<"alpha"<<'\t'<<alpha<<'\n';
  finfo<<"cycle"<<'\t'<<cycle<<'\n';




  for (size_t tj = 0; tj < TJcycle; tj++) {

    for (size_t i = 0; i < Nx*Ntau; i++) {
      corrSum[i]=0.0;
    }//zero the sum
    cSum=0;

    for(int i=0; i<cycle; i++){
       next_2d_Erdos_L(spin,Nx,Ntau,K,Gamma,alpha);
       fmag<<i+currentCycle+tj*cycle<<'\t'<<magPerSite(spin, Nx*Ntau)<<'\n';

       if(i%20==0){
        corr_2d(spin,Nx,Ntau,corr);
        for (size_t i1 = 0; i1 < Nx; i1++) {
          for (size_t i2 = 0; i2 < Ntau; i2++) {
            corrSum[i2+i1*Ntau]+=corr[i2+i1*Ntau];
          }
        }
        cSum+=1;

       }

       showtime();
       std::cout <<"STD"<<tj<<'\t'  <<i <<" in "<<cycle << '\t'<< i*100.0/cycle<< '%' << '\n'<<'\n';
    }

    fcor<<"#cycle"<<'\t'<<currentCycle+tj*cycle<<'\n';
    for (size_t i1 = 0; i1 < Nx; i1++) {
      for (size_t i2 = 0; i2 < Ntau; i2++) {
        fcor<<corrSum[i2+i1*Ntau]/cSum<<'\t';
      }
      fcor<<'\n';
    }


  }






    fspin<<"#cycle"<<'\t'<<cycle*TJcycle+currentCycle<<'\n';
    for (size_t i1 = 0; i1 < Nx; i1++) {
      for (size_t i2 = 0; i2 < Ntau; i2++) {
        fspin<<spin[i2+i1*Ntau]<<'\t';
      }
      fspin<<'\n';
    }


  finfo.close();
  fmag.close();
  fcor.close();
  fspin.close();
  return 0;
}
