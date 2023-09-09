#pragma once
#include <iostream>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iomanip>
#include <chrono>
#include <cmath>
#include "idb.h"
#include "grid.h"
#include "eos.h"
#include "cell.h"
#include "init.h"
#include "hydro.h"
#include "eos0.h"



using std::cout;
using std::endl;
using namespace std::chrono;

class master {

public:

  master(idb* );
  ~master();
  void initialize();
  void run_hydro();
  EoS* eos;
  idb* IDB;
  grid* g;
  hydro* h;
  cnvrt* CN;
  init* ic ; 


private:
  void write_grid_info(double t);

} ;






