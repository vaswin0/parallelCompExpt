#include "cell.h"

//using std::cout;
//using std::endl;


// constructor
__host__ __device__ cell::cell(){
  rho = 0. ; 
  flux = 0 ; 
}

// destructor
__host__ __device__ cell::~cell(){
}

__host__ __device__ void cell::set_rho(double xx){
  rho = xx ; 
}

__host__ __device__ double cell::get_rho(){
  return rho ; 
}

__host__ __device__ void cell::add_flux(double xx){
  flux += xx ; 
} 


__host__ __device__ void cell::clear_flux(){
  flux = 0 ; 
}

__host__ __device__ void cell::update_rho(){
  rho += flux ; 
}
















