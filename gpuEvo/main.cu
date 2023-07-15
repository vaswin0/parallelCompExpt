#include <iostream>
//#include "eos.h"
#include "cell.h"
#include "grid.h"
#include "init.h"
#include "evolve.h"

const int N = 1000000;


__global__ void evolvee(evolve *Ev){
   int i= blockIdx.x*blockDim.x + threadIdx.x;
   if(i<N){
		        
        Ev->calc_flux(i, 0,0);
		__syncthreads();
		Ev->getCellUpdateRho(i, 0,0);
		__syncthreads();
		Ev->getCellClearFlux(i,0,0);

            }
        }			




int main(){

  double xmin = -12 ; 
  double xmax =  12 ;
  int    nx   = N;  
  double ymin = 0 ; 
  double ymax = 0 ;
  int    ny   = 1 ;  
  double zmin = 0 ; 
  double zmax = 0 ;
  int    nz   = 1 ;  

grid  Gr = grid(xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz );
  
 grid* pGr;

 cudaMallocManaged(&pGr, sizeof(grid));
 *pGr = Gr;


init *in = new init(pGr);
in->set_init();

std::cout<<Gr<<std::endl;



evolve Ev = evolve(pGr, 1.0,0.01);
evolve* pEv;
cudaMallocManaged(&pEv, sizeof(evolve));
*pEv = Ev;


evolvee<<<4,26>>>(pEv);
cudaDeviceSynchronize();
            //pR->displayArray();

std::cout<<Gr;
return 0;



}
