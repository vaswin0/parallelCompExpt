#include "evolve.h"


__host__ __device__ evolve::evolve(grid* gr_, double tau_ , double dt_):gr(gr_), tau(tau_), dt(dt_){
}



__host__ __device__ evolve::evolve(){gr = NULL;tau = 1.0; dt = 0.01;};

__host__ __device__ evolve::~evolve(){
}
/*
void evolve::evolveIt(){

  for(int ix=0; ix<gr->get_nx(); ix++ )
    for(int iy=0; iy<gr->get_ny(); iy++ )
      for(int iz=0; iz<gr->get_nz(); iz++ ){
     calc_flux(ix,iy,iz);
  }

 
  for(int ix=0; ix<gr->get_nx(); ix++ )
    for(int iy=0; iy<gr->get_ny(); iy++ )
      for(int iz=0; iz<gr->get_nz(); iz++ ){
        gr->get_cell(ix,iy,iz)->update_rho() ; 
        gr->get_cell(ix,iy,iz)->clear_flux() ; 
   }

  tau += dt ; 


}

*/
__host__ __device__ void evolve::calc_flux(int ix, int iy, int iz){

 double adv_vel = 0.1 ; 
 double flux1 = - adv_vel * dt * gr->get_cell(ix+1,iy,iz)->get_rho() / 2. / gr->get_dx() ;   
 double flux2 =   adv_vel * dt * gr->get_cell(ix-1,iy,iz)->get_rho() / 2. / gr->get_dx() ;   
 gr->get_cell(ix,iy,iz)->add_flux(flux1) ; 
 gr->get_cell(ix,iy,iz)->add_flux(flux2) ; 
printf("%f", flux1 + flux2);
}
