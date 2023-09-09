#pragma once
#include<iostream>
#include<cmath>
#include "grid.h"
#include "idb.h"


using namespace std;

class init
{
  
  private:
    idb *IDB ;


  public:

  init(idb* );
  ~init();
  void set_init(grid* , EoS* );


};

