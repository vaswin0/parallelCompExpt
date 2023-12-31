#include "master.h"
#include "eos0.h"
#include "eos.h"

master::master(idb* _IDB)
{
  // Creating the "hydro_output" directory if it doesn't exist
  if (mkdir("outputs", 0777) == -1)
    cerr << "[Info] " << strerror(errno) << endl;
  else
    cout << "[Info] hydro_output directory is created" << endl;

  IDB = _IDB;
//  eos = new EoS0();
 EoS0 eoss = EoS0();
 cudaMallocManaged(&eos, sizeof(EoS0));
 memcpy(eos,&eoss, sizeof(EoS0));
  //CN = new cnvrt();
  cnvrt cc = cnvrt();
  cudaMallocManaged(&CN, sizeof(cnvrt));
   memcpy(CN,&cc, sizeof(cnvrt));

  
}

master::~master()
{
 // delete IDB;
 cudaFree(IDB);
 // delete eos;
 cudaFree(eos);
  //delete h;
  cudaFree(h);
  //delete g;
 // delete CN;
 cudaFree(g);
 cudaFree(CN);
}

void master::initialize()
{
cout<< "in class master, initializing" <<endl;
  // make the grid
  //g = new grid(IDB,CN,eos);
  grid gr = grid(IDB,CN,eos);
  cudaMallocManaged(&g, sizeof(grid));
  memcpy(g,&gr, sizeof(grid));

  g->make_grid();

  //ic = new init(IDB);

  init in = init(IDB);
  
  cudaMallocManaged(&ic, sizeof(init));
  memcpy(ic,&in, sizeof(init));




  ic->set_init(g,eos);
    
}



void master::run_hydro(){

  double current_tau ; 
  
  //h = new hydro(eos, g , IDB , IDB->tau0 , IDB->dtau, CN);
  
  hydro hy =  hydro(eos, g , IDB , IDB->tau0 , IDB->dtau, CN);
  
  
  cudaMallocManaged(&h, sizeof(hydro));
  memcpy(h,&hy, sizeof(hydro));
  
  
  int nstep = (IDB->tauMax - IDB->tau0)/IDB->dtau + 1;
  std::cout << "Hydro evolution ..." << std::endl ;
  std::cout << "evolution till max. tau : " << IDB->tauMax << std::endl ;   
  for(int istep=0; istep<nstep ; istep++){
    current_tau =  IDB->tau0 + istep * IDB->dtau ;
  
    if(istep%10==0){ 
    std::cout << "tau : " << current_tau << std::endl ;
    } 

    if(istep%100 == 0 ){
      write_grid_info(current_tau) ; 
    }

    // entire hydro evolution here
    h->evolve();      
  }// step loop
}




void master::write_grid_info(double t){

  std::cout << "storing fluid information at tau :"
       << t << " (fm)" << std::endl ; 

  std::ofstream mFile;
  std::stringstream output_filename;

  output_filename.str("");
  output_filename << "outputs/fluid_info_at_tau_" << t << ".txt";

  mFile.open(output_filename.str().c_str(), std::ios::out );
  double epsilon, pressure, vx, vy, vz, nb, nq, ns ; 
  for(int ix=0; ix<IDB->nx; ix++){
    for(int iy=0; iy<IDB->ny; iy++){
       double xx = g->get_x(ix);;
       double yy = g->get_y(iy);;
       double rr = sqrt(xx*xx+yy*yy);
       g->get_cell(ix,iy,0)->get_physical_var(eos, t, epsilon, pressure, nb, nq, ns, vx, vy, vz); 
       mFile << xx << "  " << yy << "  " << rr << "  " << epsilon << "  " << nb << "  " << vx << "  " << vy << std::endl ; 
    }
  }

  mFile.close() ; 
}


















