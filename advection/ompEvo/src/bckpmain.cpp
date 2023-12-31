#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <chrono>
#include <cmath>
#include <omp.h>
#include "cell.h"
#include "grid.h"
#include "init.h"
#include "eos.h"
#include "evolve.h"

using std::cout;
using std::endl;
using namespace std::chrono;


int main(int argc, char **argv)
{
  double xmin = -12 ; 
  double xmax =  12 ;
  int    nx   =  241000 ;  
  double ymin = 0 ; 
  double ymax = 0 ;
  int    ny   = 1 ;  
  double zmin = 0 ; 
  double zmax = 0 ;
  int    nz   = 1 ;  

  // setting up the eos //
  eos* EoS = new eos();

  // setting up the grid //
  grid* gr = new grid(xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz );

  // setting up the initial condition //
  init* in = new init(gr);
  in->set_init() ; 

  std::ofstream file ;
  file.open("init_dist.dat");
  for(int ix=0; ix<gr->get_nx(); ix++ )
    for(int iy=0; iy<gr->get_ny(); iy++ )
      for(int iz=0; iz<gr->get_nz(); iz++ ){
      file << gr->get_cell_x_cordinate(ix,iy,iz) << "  "
           << gr->get_cell_y_cordinate(ix,iy,iz) << "  "
           << gr->get_cell_z_cordinate(ix,iy,iz) << "  "
           << gr->get_cell(ix,iy,iz)->get_rho() 
           << std::endl ;
  }
  file.close();


  evolve* ev = new evolve(gr, 1.0, 0.01);

auto start = high_resolution_clock::now()    ;   

  for(int it=0; it < 3000; it++){
    ev->evolveIt() ; 
  }


  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> duration = end - start;

  std::cout << "Loop execution time: " << duration.count    () * 1000 << " milliseconds" << std::endl;


  file.open("final_dist.dat");
  for(int ix=0; ix<gr->get_nx(); ix++ )
    for(int iy=0; iy<gr->get_ny(); iy++ )
      for(int iz=0; iz<gr->get_nz(); iz++ ){
      file << gr->get_cell_x_cordinate(ix,iy,iz) << "  "
           << gr->get_cell_y_cordinate(ix,iy,iz) << "  "
           << gr->get_cell_z_cordinate(ix,iy,iz) << "  "
           << gr->get_cell(ix,iy,iz)->get_rho() 
           << std::endl ;
  }
  file.close();


  delete in  ; 
  delete gr  ; 
  delete EoS ; 
  return 0;
}













