#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <chrono>
#include <cmath>
#include "idb.h"
#include "read_id.h"
#include "master.h"
#include <fstream>


using std::cout;
using std::endl;
using namespace std::chrono;


int main(int argc, char **argv)
{
  
  // reading the input data from the file
  string input_file_name;char* event_no_s ;

  if(argc == 3){
      event_no_s = *(argv+1); input_file_name = *(argv+2);}
  else
     {
      cout<<"[Info] Please give 2 arguments\n"
            "1st argument - event no\n    "
            "2nd argument - input filename"<<endl;
      exit(1);
     }

read_id* reader;
cudaMallocManaged(&reader, sizeof(read_id));


  //read_id* reader = new read_id(); 
  //idb *IDB = new idb;

  idb *IDB;
  cudaMallocManaged(&IDB, sizeof(idb));
  reader->read_id_from_file(IDB, input_file_name); // read input data base and store it.
  int event_no = atof(event_no_s) ;

  cout<<"event no : "<<event_no<<endl;


  master head =  master(IDB);

  //head.initialize();

  master *pHd;
  cudaMallocManaged(&pHd, sizeof(master));
  *pHd = head;
  //pHd->initialize();
    //  cout <<"****initialized from head***"<<endl;

 //EoS0 eoss = EoS0();
 EoS0* eos;
  cudaMallocManaged(&eos, sizeof(EoS0));
  *eos = eoss;
  cout<<eos->pressure(1.0,1.0,1.0,1.0)<<endl;

  //head->initialize(); // initialize
 // head.run_hydro(); // run hydro 
 

 

cudaFree(pHd);
 return 0;
}


