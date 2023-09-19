#pragma once
#include<iostream>
#include<cmath>
#include "eos.h"


using namespace std;

class cell
{
  
  private:

  // density (conserved quantity)
  double rho ;
  double flux ; 

  public:

__host__ __device__   cell();
__host__ __device__   ~cell();
__host__ __device__   void set_rho(double  );
__host__ __device__   double get_rho( );
__host__ __device__   void add_flux(double ) ; 
__host__ __device__   void clear_flux() ; 
__host__ __device__   void update_rho() ; 

};








