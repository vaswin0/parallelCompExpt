/////////////////////////////////////////////////////////
//     The viscous part is not                         //
//     following the tilde notation                    //
//     as written in vHLLE.                            //
//     means pi^{\eta \eta} != tau*tau*pi^{\eta \eta}  //
/////////////////////////////////////////////////////////


#pragma once

#include <cmath>
//#include "TMath.h"
#include "cnvrt.h"
#include "grid.h"
#include "cell.h"
//#include "TMath.h"
#include "global.h"
#include <fstream>



using std::cout;
using std::endl;

class hydro
{

public:
  hydro(EoS* , grid* ,idb* ,double , double ,cnvrt*  );
  ~hydro();
  void evolve();
  void hlle_flux(cell* left, cell* right, int direction, int mode);
  void sourcestep(int mode, int ix, int iy,int iz, double _tau);
  inline double get_tau() { return tau; }
  void set_dtau(double deltaTau);  // change the timestep
  inline double get_dtau() { return dt; }  // return current value of timestep


private:
  grid *f;
  EoS *eos;
  double dt;
  double tau;
  idb* IDB;
  cnvrt* CN;



double sign(double x) {
 if (x > 0)
  return 1.;
 else if (x < 0.)
  return -1.;
 else
  return 0.;
}


};

