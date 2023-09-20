#include <iostream>
//#include "eos.h"
#include "cell.h"
#include "grid.h"
#include "init.h"
#include "evolve.h"
using namespace std;
const int N = 50;
//__device__ grid Gr = grid(-12,12,50,0,0,1,0,0,1);


__global__ void funct(cell* c, grid* g){

g->Cell = c;

printf("%f", c[1].get_rho());
}


int main(){
	
	int nx = N;
	int ny = 1;
	int nz = 1;

	cell *pCell;
	
	size_t cellSize = nx*ny*nz*sizeof(cell);

	cudaMalloc(&pCell, cellSize); // cout<<"cudaMalloc(&pCell, cellSize) worked"<<endl;

	grid Gr = grid(-12,12,N,0,0,1,0,0,1);
//	cout << Gr.Cell<<"from host" <<endl;

	grid* pGr_d;

	
	init *in = new init(&Gr);// cout<<"made init obj"<<endl;
	in->set_init(); //cout<< "initialized" <<endl;
//	cout<<Gr.Cell[1].get_rho()<<endl;

	cudaMemcpy(pCell, &Gr.Cell, cellSize, cudaMemcpyHostToDevice);
	
	size_t gridSize = sizeof(grid);
	cudaMalloc(&pGr_d, gridSize);
	cudaMemcpy(pGr_d, &Gr, gridSize, cudaMemcpyHostToDevice);

	evolve Ev = evolve(pGr_d, 1,0.01);
	evolve* pEv_d;
	size_t evosz = sizeof(Ev);
	cudaMalloc(&pEv_d, evosz);
	cudaMemcpy(pEv_d, &Ev, evosz, cudaMemcpyHostToDevice);

	funct<<<1,1>>>(pCell, pGr_d);
	cudaDeviceSynchronize();	
	   
	return 0;



}
