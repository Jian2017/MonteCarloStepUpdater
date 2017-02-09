//#include "next.h"

//WARNING!!
//recursion method is used
//for large system stack overflow may happen
//reset recursion limit, or stack size for your mechine!



#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <math.h>

#include <ctime>

void showtime(){
  time_t now = time(0);
  char* dt = ctime(&now);
  std::cout  << dt  ;
}

#include "rng.h"

#define PI 3.14159265

int mod(int a, int b){
  while (a<0) {
    a=a+b;
  }
  return a%b;
}


int initialize_1d(int arr[], int Ntau){

  for (int i = 0; i < Ntau; i++) {
    arr[i]=1-dis01(gen)*2;
  }

  return 0;
}

int initialize_2d(int arr[], int Nx, int Ntau,const char *filename) {

  std::cout << filename << '\n';

  char tempChar[100];
  int tempInt;

  std::ifstream f1;
  f1.open(filename);
  if(f1.is_open()){
    std::cout << "open existing Zspin.txt" << '\n';
    f1>>tempChar;
    f1>>tempInt;
    std::cout << "last cycle :  "<<tempInt << '\n';

    for(int i=0 ;i <Nx*Ntau; i++){
      f1>>arr[i];
    }
    f1.close();
    return tempInt;
  }
  else {
    std::cout << "failed to open, random initialize with 50% probablity" << '\n';
    for(int i=0 ;i <Nx*Ntau; i++){
      arr[i]=1-dis01(gen)*2;
    }
    return 0;
  }




}








double next_1d_Metropolis_S(int arr[], int Ntau, double Gamma){
  int tau=int(dis(gen)*Ntau);
  tau=mod(tau,Ntau);
  int tauR=mod(tau+1,Ntau);
  int tauL=mod(tau-1,Ntau);
  //change of energy for provisional flip of spin[tau]
  double dS=Gamma*2.0*(arr[tauR]+arr[tauL])*arr[tau];
  if ( dS>0 && dis(gen)>exp(-dS) ) {
      return 0.0;
   } else {
    arr[tau]=-arr[tau];
      return 2.0*arr[tau];
  }
}




double next_1d_Metropolis_L(int arr[], int Ntau, double Gamma,double alpha){
  int tau=int(dis(gen)*Ntau);
  tau=mod(tau,Ntau);
  int tauR=mod(tau+1,Ntau);
  int tauL=mod(tau-1,Ntau);
  //change of energy for provisional flip of spin[tau]
  double dS=Gamma*2.0*(arr[tauR]+arr[tauL])*arr[tau];

  double mf=0.0;
  for (int j = 0; j < Ntau; j++) {
    if(j != tau){
      mf=mf+arr[j]/(sin(PI*(j-tau)/Ntau)*sin(PI*(j-tau)/Ntau));
    }
  }
  mf=mf*(PI/Ntau)*(PI/Ntau)*alpha*arr[tau];
  dS=dS+mf;
  if ( dS>0 && dis(gen)>exp(-dS) ) {
    return 0.0;
   } else {
    arr[tau]=-arr[tau];
    return 2.0*arr[tau];
  }

}




double next_2d_Metropolis_S(int arr[], int Nx, int Ntau, double K, double Gamma ){
  int tau=int(dis(gen)*Ntau);
  tau=mod(tau,Ntau);
  int tauR=mod(tau+1,Ntau);
  int tauL=mod(tau-1,Ntau);

  int x=int(dis(gen)*Nx);
  x=mod(x,Nx);
  int xR=mod(x+1,Nx);
  int xL=mod(x-1,Nx);


  double dS=2.0*(Gamma*(arr[x*Ntau+tauR]+arr[x*Ntau+tauL])+K*(arr[xR*Ntau+tau]+arr[xL*Ntau+tau]))*arr[x*Ntau+tau];
  if ( dS>0 && dis(gen)>exp(-dS) ) {
    return 0.0;

   } else {
    arr[x*Ntau+tau]=-arr[x*Ntau+tau];
    return 2.0*arr[x*Ntau+tau];

  }

}






double next_2d_Metropolis_L(int arr[], int Nx, int Ntau, double K, double Gamma, double alpha ){
  int tau=int(dis(gen)*Ntau);
  tau=mod(tau,Ntau);
  int tauR=mod(tau+1,Ntau);
  int tauL=mod(tau-1,Ntau);

  int x=int(dis(gen)*Nx);
  x=mod(x,Nx);
  int xR=mod(x+1,Nx);
  int xL=mod(x-1,Nx);


  double dS=2.0*(Gamma*(arr[x*Ntau+tauR]+arr[x*Ntau+tauL])+K*(arr[xR*Ntau+tau]+arr[xL*Ntau+tau]))*arr[x*Ntau+tau];

  double mf=0.0;
  for (int j = 0; j < Ntau; j++) {
    if(j != tau){
      mf=mf+arr[x*Ntau+j]/(sin(PI*(j-tau)/Ntau)*sin(PI*(j-tau)/Ntau));
    }
  }
  mf=mf*(PI/Ntau)*(PI/Ntau)*alpha*arr[tau];
  dS=dS+mf;

  if ( dS>0 && dis(gen)>exp(-dS) ) {
    return 0.0;
   } else {
    arr[x*Ntau+tau]=-arr[x*Ntau+tau];
    return 2*arr[x*Ntau+tau];
  }
}










//*************************************cluster algotrim********//

//****** double next_1d_Cluster_S(int arr[], int Ntau, double Gamma)

int condition(int arr[], int cluster[],int clusterSpin,int testSite,double Gamma){
  if((cluster[testSite]==0)&&(arr[testSite]==clusterSpin)&&(dis(gen)> exp(-2*Gamma))){
    return 1;
  }
  else{
    return 0;
  }
}

void growCluster(int arr[],int cluster[],int ic,int clusterSpin,int Ntau,double Gamma){
  if(cluster[ic]==1){
    return ;
  }
  else{
    cluster[ic]=1;
    arr[ic]=-arr[ic];
    if(    condition(arr,cluster,clusterSpin,mod(ic-1,Ntau),Gamma)){
      growCluster(arr,cluster,mod(ic-1,Ntau),clusterSpin,Ntau,Gamma);
    }

    if(    condition(arr,cluster,clusterSpin,mod(ic+1,Ntau),Gamma)){
      growCluster(arr,cluster,mod(ic+1,Ntau),clusterSpin,Ntau,Gamma);
    }

    return ;
  }
}


double next_1d_Cluster_S(int arr[], int Ntau, double Gamma){
    int cluster[Ntau] ;
    for (int i = 0; i < Ntau; i++) {
        cluster[i]=0;
    }

    int ic=int(dis(gen)*Ntau);
    int clusterSpin=arr[ic];

    growCluster(arr,cluster,ic,clusterSpin,Ntau,Gamma);

    return 1.0;


}











//********************************
//**  double next_2d_Cluster_S(int arr[],int Nx, int Ntau,double K, double Gamma)

void growCluster(int arr[],int cluster[],int ix,int itau,int clusterSpin,int Nx,int Ntau,double K,double Gamma){

  if(cluster[ix*Ntau+itau]==1){
    return ;
  }
  else{
    cluster[ix*Ntau+itau]=1;
    arr[ix*Ntau+itau]=-arr[ix*Ntau+itau];

    if(    condition(arr,cluster,clusterSpin,ix*Ntau+mod(itau-1,Ntau),Gamma)){
      growCluster(arr,cluster,ix,mod(itau-1,Ntau),clusterSpin,Nx,Ntau,K,Gamma);
    }

    if(    condition(arr,cluster,clusterSpin,ix*Ntau+mod(itau+1,Ntau),Gamma)){
      growCluster(arr,cluster,ix,mod(itau+1,Ntau),clusterSpin,Nx,Ntau,K,Gamma);
    }

    if(    condition(arr,cluster,clusterSpin,mod(ix+1,Ntau)*Ntau+itau,K)){
      growCluster(arr,cluster,mod(ix+1,Ntau),itau,clusterSpin,Nx,Ntau,K,Gamma);
    }

    if(    condition(arr,cluster,clusterSpin,mod(ix-1,Ntau)*Ntau+itau,K)){
      growCluster(arr,cluster,mod(ix-1,Ntau),itau,clusterSpin,Nx,Ntau,K,Gamma);
    }

    return ;
  }
}


double next_2d_Cluster_S(int arr[],int Nx, int Ntau,double K, double Gamma){
    int cluster[Ntau*Nx] ;
    for (int i = 0; i < Ntau*Nx; i++) {
        cluster[i]=0;
    }

    int itau=int(dis(gen)*Ntau);
    int ix=int(dis(gen)*Nx);
    int clusterSpin=arr[ix*Ntau+itau];

    growCluster(arr,cluster,ix,itau,clusterSpin,Nx,Ntau,K,Gamma);

    return 1.0;
}



















//****** double next_1d_Cluster_L(int arr[], int Ntau, double Gamma)
std::vector<double> v;

void  load_vector(int Ntau, double alpha) {
  //the size of loaded v is from 1,2,3, to Ntau-1
  // with v[index]  from 0,1,2 to Ntau-2
    double temp=0.0;
    for (int i = 1; i < Ntau; i++) {
      temp=temp+alpha*(PI/Ntau)*(PI/Ntau)/(sin(PI*i/Ntau)*sin(PI*i/Ntau));
      v.push_back(1-exp(-temp));
    }
    return ;
}


void growCluster(int arr[],int cluster[],int ic,int clusterSpin,int Ntau,double Gamma,double alpha){

  if(cluster[ic]==1){
    return ;
  }
  else{
    cluster[ic]=1;
    arr[ic]=-arr[ic];

    if(    condition(arr,cluster,clusterSpin,mod(ic-1,Ntau),Gamma)){
      growCluster(arr,cluster,mod(ic-1,Ntau),clusterSpin,Ntau,Gamma, alpha);
    }

    if(    condition(arr,cluster,clusterSpin,mod(ic+1,Ntau),Gamma)){
      growCluster(arr,cluster,mod(ic+1,Ntau),clusterSpin,Ntau,Gamma, alpha);
    }

    double g=dis(gen);
    std::vector<double>::iterator low=std::lower_bound(v.begin(),v.end(),g);
    int it= std::distance(v.begin(), low)+1;
    if(   (it<Ntau)&&(arr[mod(ic+it,Ntau)]==clusterSpin)   ){
      growCluster(arr,cluster,mod(ic+it,Ntau),clusterSpin,Ntau,Gamma, alpha);
     }


    while (it<Ntau) {
      g=dis(gen)*(1-v[it-1])+v[it-1];
      low=std::lower_bound(v.begin(),v.end(),g);
      it=low-v.begin()+1;
      if(   (it<Ntau)&&(arr[mod(ic+it,Ntau)]==clusterSpin)   ){
        growCluster(arr,cluster,mod(ic+it,Ntau),clusterSpin,Ntau,Gamma, alpha);
       }
    }

    return ;
  }
}


double next_1d_Cluster_L(int arr[], int Ntau, double Gamma,double alpha){
    int cluster[Ntau] ;
    for (int i = 0; i < Ntau; i++) {
        cluster[i]=0;
    }
    int ic=int(dis(gen)*Ntau);
    int clusterSpin=arr[ic];

    load_vector(Ntau,alpha);
    growCluster(arr,cluster,ic,clusterSpin,Ntau,Gamma,alpha);
    v.clear();
    return 1.0;
}















//********************************
//**  double next_2d_Cluster_L(int arr[],int Nx, int Ntau,double K, double Gamma,double alpha)

void growCluster(int arr[],int cluster[],int ix,int itau,int clusterSpin,int Nx,int Ntau,double K,double Gamma,double alpha){

  if(cluster[ix*Ntau+itau]==1){
    return ;
  }
  else{
    cluster[ix*Ntau+itau]=1;
    arr[ix*Ntau+itau]=-arr[ix*Ntau+itau];

    if(    condition(arr,cluster,clusterSpin,ix*Ntau+mod(itau-1,Ntau),Gamma)){
      growCluster(arr,cluster,ix,mod(itau-1,Ntau),clusterSpin,Nx,Ntau,K,Gamma,alpha);
    }

    if(    condition(arr,cluster,clusterSpin,ix*Ntau+mod(itau+1,Ntau),Gamma)){
      growCluster(arr,cluster,ix,mod(itau+1,Ntau),clusterSpin,Nx,Ntau,K,Gamma,alpha);
    }

    if(    condition(arr,cluster,clusterSpin,mod(ix+1,Ntau)*Ntau+itau,K)){
      growCluster(arr,cluster,mod(ix+1,Ntau),itau,clusterSpin,Nx,Ntau,K,Gamma,alpha);
    }

    if(    condition(arr,cluster,clusterSpin,mod(ix-1,Ntau)*Ntau+itau,K)){
      growCluster(arr,cluster,mod(ix-1,Ntau),itau,clusterSpin,Nx,Ntau,K,Gamma,alpha);
    }


    double g=dis(gen);
    std::vector<double>::iterator low=std::lower_bound(v.begin(),v.end(),g);
    int it= std::distance(v.begin(), low)+1;
    if(   (it<Ntau)&&(arr[ix*Ntau+mod(itau+it,Ntau)]==clusterSpin)   ){
      growCluster(arr,cluster,ix,mod(itau+it,Ntau),clusterSpin,Nx,Ntau,K,Gamma, alpha);
     }


    while (it<Ntau) {
      g=dis(gen)*(1-v[it-1])+v[it-1];
      low=std::lower_bound(v.begin(),v.end(),g);
      it=low-v.begin()+1;
      if(   (it<Ntau)&&(arr[ix*Ntau+mod(itau+it,Ntau)]==clusterSpin)   ){
        growCluster(arr,cluster,ix,mod(itau+it,Ntau),clusterSpin,Nx,Ntau,K,Gamma, alpha);
       }
    }

    return ;
  }
}


double next_2d_Cluster_L(int arr[], int Nx, int Ntau, double K, double Gamma, double alpha){
    int cluster[Ntau*Nx] ;
    for (int i = 0; i < Ntau*Nx; i++) {
        cluster[i]=0;
    }

    int itau=int(dis(gen)*Ntau);
    int ix=int(dis(gen)*Nx);
    int clusterSpin=arr[ix*Ntau+itau];

    load_vector(Ntau,alpha);
    growCluster(arr,cluster,ix,itau,clusterSpin,Nx,Ntau,K,Gamma,alpha);
    v.clear();

    return 1.0;
}




//*****************************************************************
//****Get Rid of Recursion ****************************************
//*****************************************************************
//****Erdos Number Cluster Labeling method ******** new method !***
//****by Jian Wang, Jan 19th 2017 *********************************
//*****************************************************************

double next_2d_Erdos_L(int arr[], int Nx, int Ntau, double K, double Gamma, double alpha){

//  int erdos[Ntau*Nx]; //************************* function global variable
  int* erdos = new int[Ntau*Nx];

  for (int i = 0; i < Ntau*Nx; i++) {
      erdos[i]=0;
  }

  int itau=int(dis(gen)*Ntau);
  int ix=int(dis(gen)*Nx);

  load_vector(Ntau,alpha);

  int clusterSpin=arr[ix*Ntau+itau];//*********** function gloabl variable

  erdos[ix*Ntau+itau]=1;  //label the core with erdosNumber 1

  int count=1;
  int erdosNumber=1;

  while(count>0){
    count=0;
    for (int i = 0; i < Nx*Ntau; i++) {
      if (erdos[i]==erdosNumber) {

        count=count+1;
        /* code */
        //site i with erdosNumber, find its (erdosNumber+1) neigbour sites
        //label these new sites with erods[newsites]=(erdosNumber+1)
        {
          int iitau=i%Ntau;
          int iix=i/Ntau;

          if(    condition(arr,erdos,clusterSpin,iix*Ntau+mod(iitau-1,Ntau),Gamma)){
            erdos[iix*Ntau+mod(iitau-1,Ntau)]=erdosNumber+1;
          }
          if(    condition(arr,erdos,clusterSpin,iix*Ntau+mod(iitau+1,Ntau),Gamma)){
            erdos[iix*Ntau+mod(iitau+1,Ntau)]=erdosNumber+1;
          }
          if(    condition(arr,erdos,clusterSpin,mod(iix-1,Nx)*Ntau+iitau,K)){
            erdos[mod(iix-1,Nx)*Ntau+iitau]=erdosNumber+1;
          }
          if(    condition(arr,erdos,clusterSpin,mod(iix+1,Nx)*Ntau+iitau,K)){
            erdos[mod(iix+1,Nx)*Ntau+iitau]=erdosNumber+1;
          }


          double g=dis(gen);
          std::vector<double>::iterator low=std::lower_bound(v.begin(),v.end(),g);
          int it= std::distance(v.begin(), low)+1;
          if(   (it<Ntau)&&(arr[ix*Ntau+mod(itau+it,Ntau)]==clusterSpin) &&(erdos[iix*Ntau+mod(iitau+it,Ntau)]==0)  ){
              erdos[iix*Ntau+mod(iitau+it,Ntau)]=erdosNumber+1;
           }
          while (it<Ntau) {
            g=dis(gen)*(1-v[it-1])+v[it-1];
            low=std::lower_bound(v.begin(),v.end(),g);
            it=low-v.begin()+1;
            if(   (it<Ntau)&&(arr[ix*Ntau+mod(itau+it,Ntau)]==clusterSpin) &&(erdos[iix*Ntau+mod(iitau+it,Ntau)]==0)  ){
                erdos[iix*Ntau+mod(iitau+it,Ntau)]=erdosNumber+1;
             }
          }

        }






      }
    }
    erdosNumber=erdosNumber+1;
  }



  v.clear();


  for (size_t i = 0; i < Nx*Ntau; i++) {
    if (erdos[i]>0) {
      arr[i]=-arr[i];
    }
  }

  delete [] erdos;

  return 3.14159;
}
