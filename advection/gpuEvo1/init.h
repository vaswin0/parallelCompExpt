#pragma once
#include<iostream>
#include<cmath>
#include "grid.h"



using namespace std;

class init
{
  
  private:
    grid* g ;


  public:

__host__ __device__  init(grid* );
__host__ __device__  ~init();
__host__ __device__  void set_init();

};

