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
  
__host__ __device__ inline double get_tau() {return tau;};
__host__ __device__   evolve(grid*, double, double  );
__host__ __device__   ~evolve();
//__host__ __device__   void evolveIt();
__host__ __device__   void calc_flux(int, int, int);
};


