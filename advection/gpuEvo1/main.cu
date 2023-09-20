#include <iostream>
//#include "eos.h"
#include "cell.h"
#include "grid.h"
#include "init.h"
#include "evolve.h"
using namespace std;
const int N = 50;
//__device__ grid Gr = grid(-12,12,50,0,0,1,0,0,1);


__global__ void initialize(evolve* Ev){

	grid* gr = new grid (-12,12,50,0,0,1,0,0,1);
	init* in = new init(gr); in->set_init();
	Ev = new evolve(gr, 1.0,0.01);




}

__global__ void evolvee(evolve* Ev){


	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {Ev->calc_flux(i,0,0);printf("%d", i);};
}

int main(){
	
	
	evolve* pEv;
	cudaMalloc(&pEv, sizeof(evolve));


	initialize<<<1,1>>>(pEv);

	cudaDeviceSynchronize();	

	cudaMemcpy(var_d, var, sizeof(int), cudaMemcpyDeviceToDevice);

	
	evolvee<<<5,10>>>(pEv);
	cudaDeviceSynchronize();
	   
	return 0;



}
