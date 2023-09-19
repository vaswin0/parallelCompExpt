#include <iostream>
//#include "eos.h"
#include "cell.h"
#include "grid.h"
#include "init.h"
#include "evolve.h"
using namespace std;
const int N = 50;
//__device__ grid Gr = grid(-12,12,50,0,0,1,0,0,1);

int main(){
	
	int nx = N;
	int ny = 1;
	int nz = 1;

	cell *pCell;
	


	
	size_t cellSize = nx*ny*nz*sizeof(cell);

	cudaMalloc(&pCell, cellSize);  cout<<"cudaMalloc(&pCell, cellSize) worked"<<endl;

	grid Gr = grid(-12,12,N,0,0,1,0,0,1);
	cout << Gr.Cell <<endl;
	Gr.Cell = pCell; cout<<Gr.Cell<<endl;
	
	grid *pGr_d;

	
	std::cout<<4;

	init *in = new init(&Gr); cout<<"made init obj"<<endl;
	in->set_init(); cout<< "initialized" <<endl;
//	std::cout<<Gr<<std::endl;
 	std::cout<<5;

//	evolve Ev = evolve(pGr_d, 1.0,0.01);


	size_t gridSize = sizeof(grid);
	cudaMalloc(&pGr_d, gridSize);
	cudaMemcpy(pGr_d, &Gr,gridSize, cudaMemcpyHostToDevice);
	
	evolve Ev = evolve(pGr_d, 1.0,0.01);
	evolve* pEv_d;
	size_t evoSz = sizeof(Ev);
	cudaMalloc(&pEv_d, evoSz);
	cudaMemcpy(pEv_d, &Ev, evoSz, cudaMemcpyHostToDevice);


//	evolvee<<<1,1>>>(pEv_d);
	//std::cout<<122;
//	cudaDeviceSynchronize();
	
		   
	return 0;



}
