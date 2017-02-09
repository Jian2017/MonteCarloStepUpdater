#ifndef NEXT_H
#define NEXT_H

//#include "rng.h"

//WARNING!!
//recursion method is used
//for large system stack overflow may happen
//reset recursion limit, or stack size for your mechine!

#include <iostream>
using namespace std;

void showtime();

int initialize_1d(int arr[], int Ntau);
int initialize_2d(int arr[], int Nx, int Ntau,const char *filename);

// void next(int arr[], ......) updates the input arr[]


double next_1d_Metropolis_S(int arr[], int Ntau, double Gamma);
double next_1d_Metropolis_L(int arr[], int Ntau, double Gamma, double alpha);

double next_2d_Metropolis_S(int arr[], int Nx, int Ntau, double K, double Gamma );
double next_2d_Metropolis_L(int arr[], int Nx, int Ntau, double K, double Gamma, double alpha);

double next_1d_Cluster_S(int arr[], int Ntau, double Gamma);
double next_1d_Cluster_L(int arr[], int Ntau, double Gamma, double alpha);

double next_2d_Cluster_S(int arr[], int Nx, int Ntau, double K, double Gamma );
double next_2d_Cluster_L(int arr[], int Nx, int Ntau, double K, double Gamma, double alpha);

//this is the thing!
double next_2d_Erdos_L(int arr[], int Nx, int Ntau, double K, double Gamma, double alpha);


#endif
