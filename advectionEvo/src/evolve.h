#pragma once
#include<iostream>
#include<cmath>
#include "eos.h"
#include "grid.h"


using namespace std;

class evolve
{
  
  private:
    grid* gr ; 
    double tau ; 
    double dt ; 


  public:

  evolve(grid*, double, double  );
  ~evolve();

__host__ __device__ void calc_flux(int, int, int);
__host__ __device__ void getCellUpdatRho(int, int, int);
__host__ __device__ void getCellClearFlux(int, int, int);
};


