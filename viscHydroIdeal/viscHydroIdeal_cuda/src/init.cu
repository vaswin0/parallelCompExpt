#include "init.h"
#include <cmath>


init::init(idb *_IDB){
  IDB = _IDB;
}

init::~init(){
}


void init::set_init(grid* g, EoS* eos){

  std::cout << "Initialising ... " << std::endl ; 
  double eps_0 =1.0;
  double n_0 = 0.5;
  double q = 1.0 ;
  cout << "getting tau from idb in set_init fn in init class" <<endl;
  double tau = IDB->tau0 ; 
   cout << "got tau from idb in set_init fn in init class" <<endl;
  cell* c;
  cudaMallocManaged(&c, sizeof(cell));
     cout << "made cell" <<endl;
  for(int ix=0; ix<g->get_nx(); ix++ )
    for(int iy=0; iy<g->get_ny(); iy++ )
      for(int iz=0; iz<g->get_neta(); iz++ ){
        c = g->get_cell(ix,iy,iz);  
             cout << "got cell from grid" <<endl;          
	double x_ = g->get_x(ix);
	double y_ = g->get_y(iy);
        double rt = sqrt(x_*x_+y_*y_);
        double eps = ((eps_0*pow(2.0*q,8.0/3.0))/(pow(tau,4.0/3.0)))
           *(pow((1+2*q*q*(tau*tau+rt*rt)+q*q*q*q*pow(tau*tau-rt*rt,2)),-4.0/3.0));
        double nb = ((n_0*4*q*q)/tau)*(pow((1+2*q*q*(tau*tau+rt*rt)+q*q*q*q*pow(tau*tau-rt*rt,2)),-1.0));
	  
        double k_ = atanh((2.0*q*q*tau*rt)/(1+q*q*tau*tau+q*q*rt*rt));
        double vx =  (x_/rt)*(tanh(k_)); 
        double vy =  (y_/rt)*(tanh(k_)); 
        double vz = 0.0;
        double nq = 0.; double ns=0.0;
                     cout << " about to set_prim_var using  cell pointer " <<endl;
        c->set_prim_var(eos,tau,eps, nb, nq,  ns,  vx,  vy,  vz) ; 
             cout << "set_prim_var using  cell pointer done" <<endl;
  }

  std::cout << "Gubser initialisation done ... " << std::endl ; 

}