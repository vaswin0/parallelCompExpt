#include "evolve.h"



evolve::evolve(grid* gr_, double tau_ , double dt_):gr(gr_), tau(tau_), dt(dt_){
}

evolve::~evolve(){
}


__host__ __device__ void evolve::calc_flux(int ix, int iy, int iz){

 double adv_vel = 0.1 ; 
 double flux1 = - adv_vel * dt * gr->get_cell(ix+1,iy,iz)->get_rho() / 2. / gr->get_dx() ;   
 double flux2 =   adv_vel * dt * gr->get_cell(ix-1,iy,iz)->get_rho() / 2. / gr->get_dx() ;   
 gr->get_cell(ix,iy,iz)->add_flux(flux1) ; 
 gr->get_cell(ix,iy,iz)->add_flux(flux2) ; 

}

__host__ __device__ void evolve::getCellUpdateRho(int ix, int iy, int iz){


	gr->get_cell(ix, iy, iz)->update_rho();
	
	}


__host__ __device__ void evolve::getCellClearFlux(int ix,     int iy, int iz){ 

	gr->get_cell(ix, iy, iz)->clear_flux();

}


