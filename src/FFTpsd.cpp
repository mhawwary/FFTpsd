
#include "fft_driver.h"
#include <time.h>
#include "../include/error.h"
#include <random>
#include "general_tools.h"
#include <cmath>

using namespace std;

fftspc::FFT_Driver *fft_driver_=nullptr;

void logo();

int main(int argc, char** argv){

  if (argc < 2)
    FatalError_exit("ERROR: no input file is specified ............ ");

  string in_fname = argv[argc-1];

  logo();
  fft_driver_=new fftspc::FFT_Driver(in_fname);
  fft_driver_->Init();
  fft_driver_->Compute();
  fft_driver_->DumpOutputs();

  fft_driver_->Reset();

  if(fft_driver_!=nullptr)
    delete fft_driver_;

  return 0;
}

void logo(){

    cout<<"_________________________________________________________________________________________"<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"              "<<"  Welcome to the FFT, PSD, and SPL computation toolbox   "<<"          "<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"                     Mohammad Alhawwary, PhD. Student, mhawwary@ku.edu                   "<<endl;
    cout<<"                Aerospace Engineering Department, University of Kansas, USA              "<<endl;
    cout<<"_______________________________________01/03/2019________________________________________\n"<<endl;

    return;
}


