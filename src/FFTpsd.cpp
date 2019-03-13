
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

  logo();
  if (argc ==1){
    std::cout<<"help with implemented options"<<std::endl;
    std::cout<<"--------------------------------"<<std::endl;
    std::cout<<"       minimum arguments for fft computations"<<std::endl;
    std::cout<<"       ---------------------------------------"<<std::endl;
    std::cout<<"           [-i inputfile] "<<std::endl;
    std::cout<<"           [-n N] "<<std::endl;

    std::cout<<"\n       minimum arguments for psd, power & spl computations"<<std::endl;
    std::cout<<"       ------------------------------------------------------"<<std::endl;
    std::cout<<"           [-i inputfile]"<<std::endl;
    std::cout<<"           [-n N] "<<std::endl;
    std::cout<<"           [-psd] or [-pow]"<<std::endl;

    std::cout<<"\n\nfor a detailed input list and discription use [--help]\n";

    std::cout<<"\nAlternatively, you can just parse an input parsing file as in the examples\n\n";

    _print_log("Exiting the program");
    exit(0);

  }else if(argc==2){ // use an input file
    string in_fname = argv[argc-1];
    if(in_fname==string("-h") || in_fname==string("--help"))
      fft_driver_=new fftspc::FFT_Driver(argc,argv);
    else
      fft_driver_=new fftspc::FFT_Driver(in_fname);
  }else{ // use command line argument
    fft_driver_=new fftspc::FFT_Driver(argc,argv);
  }

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


