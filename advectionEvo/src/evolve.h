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
  void evolveIt();
  void evolveItOmp();
  void calc_flux(int, int, int);
};


