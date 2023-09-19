#pragma once
#include<iostream>
#include<cmath>
#include "cell.h" 



using namespace std;

class grid
{
  
  private:

    double xmin ; 
    double xmax ; 
    double ymin ; 
    double ymax ; 
    double zmin ; 
    double zmax ; 

    int nx ; 
    int ny ; 
    int nz ; 
    
//    cell* Cell ; 

  public:
  cell* Cell;

__host__ __device__    grid(double xmin_, double xmax_, int nx_, double ymin_, double ymax_, int ny_ , double zmin_, double zmax_, double nz_);
__host__ __device__    ~grid();

__host__ __device__    inline int get_nx()   { return nx; }
__host__ __device__     inline int get_ny()   { return ny; }
__host__ __device__     inline int get_nz()   { return nz; }

__host__ __device__     inline double get_dx()   { if(nx <= 1 ){return 0;}else{return (xmax - xmin)/(nx-1);} }
__host__ __device__     inline double get_dy()   { if(ny <= 1 ){return 0;}else{return (ymax - ymin)/(ny-1);} }
__host__ __device__     inline double get_dz()   { if(nz <= 1 ){return 0;}else{return (zmax - zmin)/(nz-1);} }

__host__ __device__     cell* get_cell(int ix, int iy, int iz);
__host__ __device__     inline double get_cell_x_cordinate(int ix, int iy, int iz){ return xmin + ix * get_dx() ;  }
__host__ __device__     inline double get_cell_y_cordinate(int ix, int iy, int iz){ return ymin + iy * get_dy() ;  }
__host__ __device__     inline double get_cell_z_cordinate(int ix, int iy, int iz){ return zmin + iz * get_dz() ;  }

};





