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

  cell();
  ~cell();
  void set_rho(double  );
  double get_rho( );
  void add_flux(double ) ; 
  void clear_flux() ; 
  void update_rho() ; 

};








