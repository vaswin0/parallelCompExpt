#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif




#pragma once
#include<iostream>
#include<cmath>
//#include "eos.h"


using namespace std;

class cell
{
  
  private:

  // density (conserved quantity)
  double rho ;
  double flux ; 

  public:

__host__  cell();
__host__  ~cell();
__device__ __host__ void set_rho(double  );
__device__ __host__  double get_rho( );
__device__ __host__ void add_flux(double ) ; 
__device__ __host__  void clear_flux() ; 
__device__ __host__  void update_rho() ;
friend std::ostream& operator<< (std::ostream & os, const cell & c);

};








