#include "init.h"


init::init(grid* g_){
  g = g_ ; 
}

init::~init(){
}


void init::set_init(){

  double ss = 1.0 ; // stad. dev.
  double x0 = 0.0 ; // mean
  cell* c ; 
  for(int ix=0; ix<g->get_nx(); ix++ )
    for(int iy=0; iy<g->get_ny(); iy++ )
      for(int iz=0; iz<g->get_nz(); iz++ ){
        c = g->get_cell(ix,iy,iz); // got the cell
        double xx = g->get_cell_x_cordinate(ix,iy,iz) ; 
        double dens = sqrt( 1. / (2.*3.141*ss*ss) ) * exp( -0.5 * pow( (xx-x0)/ss, 2 ) ) ; 
        if( dens<0.0000001 ){
          dens = 0.0000001 ; 
        }
        c->set_rho(dens);
  }

  std::cout << "Initialisation done ... " << std::endl ; 

}
