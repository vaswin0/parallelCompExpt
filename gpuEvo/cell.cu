#include "cell.h"

using std::cout;
using std::endl;


// constructor
cell::cell(){
  rho = 0 ; 
  flux = 0 ; 
}

// destructor
cell::~cell(){
}

void cell::set_rho(double xx){
  rho = xx ; 
}

double cell::get_rho(){
  return rho ; 
}

void cell::add_flux(double xx){
  flux += xx ; 
} 


void cell::clear_flux(){
  flux = 0 ; 
}

void cell::update_rho(){
  rho += flux ; 
}



std::ostream& operator<<(std::ostream& os, const cell &c){

	std::cout << "flux = " << c.flux << " rho = " << c.rho << std::endl;


return os;
}











