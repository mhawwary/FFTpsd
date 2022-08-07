
#include "fft_driver.h"
#include <time.h>
#include "../include/error.h"
#include <random>
#include "general_tools.h"
#include <cmath>

using namespace std;

fftspc::FFT_Driver *fft_driver_=nullptr;

void logo();
void help_detailed();
void help_short();

int main(int argc, char** argv){

  logo();
  if (argc ==1){
    help_short();
  }else if(argc==2){
    string in_fname = argv[argc-1];
    if(in_fname==string("-h"))
      help_short();
    else if(in_fname==string("--help"))
      help_detailed();
    else{
      std::ifstream input;
      input.open(in_fname);
      if(!input.is_open()){
        std::string o_error_message_="Failed to open "+ in_fname +", or file does not exist!";
        FatalError(o_error_message_.c_str());
        help_short();
      }
      input.close(); // this was only to check if the file exist
      fft_driver_=new fftspc::FFT_Driver(in_fname);
    }
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

void help_detailed(){
  std::cout<<"help with implemented options"<<std::endl;
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"       minimum arguments for fft computations"<<std::endl;
  std::cout<<"       --------------------------------------"<<std::endl;
  std::cout<<"           [-i datafile] "<<std::endl;
  std::cout<<"           [-l time length] "<<std::endl;

  std::cout<<"\n       minimum arguments for psd, power & spl computations"<<std::endl;
  std::cout<<"       -----------------------------------------------------"<<std::endl;
  std::cout<<"           [-i datafile]"<<std::endl;
  std::cout<<"           [-l time length] "<<std::endl;
  std::cout<<"           [-psd] or [-pow] or [-spl]"<<std::endl;

  std::cout<<"\n       If DFT is preferred just add to your arguments the flag [-dft]"<<std::endl;

  std::cout<<"\nDetailed option list for more control over how the spectral analysis is performed"<<std::endl;
  std::cout<<"---------------------------------------"<<std::endl;
  std::cout<<"  [-i datafile]  name of the data file with its address"<<std::endl;
  std::cout<<"  [-o directory] output directory, default is the same as input"<<std::endl;
  std::cout<<"  [-spec] for output files with headers and spec suffix"<<std::endl;
  std::cout<<"  [-r first_row]  of the signal data in the data file, default 1"<<std::endl;
  std::cout<<"  [-q last_row ]  of the signal data in the data file, default is the end of the file"<<std::endl;
  std::cout<<"  [-c column] of the signal data in the data file, default 2, assuming time is at the first column"<<std::endl;
  std::cout<<"\n  [-n N] number of data in one window subset, or in the whole sample if no shifting/averaging."<<std::endl;
  std::cout<<"         If using fft it must be 2^{k}, k is an integer"<<std::endl;
  std::cout<<"  [-l time length] window subset length in sec or the whole sample length if no shifting/averaging"<<std::endl;
  std::cout<<"  [-s shift] a value in [0.,1.] to indicate the ratio of data to be shifted"<<std::endl;
  std::cout<<"             shift defaults: 0.5 for psd/pow/spl, 0. for fft/dft"<<std::endl;
  std::cout<<"  [-dt time_step] of the window subset or the whole sample in sec"<<std::endl;
  std::cout<<"  [-w window] name of the window functions, for psd default is hann, for fft/dft default is rectangular"<<std::endl;
  std::cout<<"              windows: hann,hamm,bartlett,welch,blackman,triangular,rectangular"<<std::endl;
  std::cout<<"  [-m number] either 1 or 0 for mean substract, default is 1 to substract the mean"<<std::endl;
  std::cout<<"\n  [-psd filename]      compute the power spectral density (psd), filename is optional"<<std::endl;
  std::cout<<"  [-pow filename]      compute the power spectrum, filename is optional"<<std::endl;
  std::cout<<"  [-spl filename]      compute the sound pressure level (SPL), filename is optional"<<std::endl;
  std::cout<<"  [-peak]     peak preserving averaging psd mode, useful for ensemble averaging with shift=1."<<std::endl;
  std::cout<<"  [-variance] variance preserving averaging psd mode for conserving the signal energy, by default it is forced for psd, power and spl computations"<<std::endl;
  std::cout<<"  [-pref pressure] reference pressure for SPL computation, default is 2e-5 Pa"<<std::endl;
  // std::cout<<"  [-bin] bin-averaging for psd and spl"<<std::endl;
  // std::cout<<"  [-oct number] the number nth octave averaging that is required, default is 1"<<std::endl;
  //std::cout<<"\nAlternatively, you can just parse an input parsing file as follows: fftpsd inputfile\n";
  exit(0);
}

void help_short(){
  std::cout<<"\nhelp with implemented options"<<std::endl;
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"       minimum arguments for fft computations"<<std::endl;
  std::cout<<"       --------------------------------------"<<std::endl;
  std::cout<<"           [-i data file] "<<std::endl;
  std::cout<<"           [-l time length] "<<std::endl;

  std::cout<<"\n       minimum arguments for psd, power & spl computations"<<std::endl;
  std::cout<<"       -----------------------------------------------------"<<std::endl;
  std::cout<<"           [-i data file]"<<std::endl;
  std::cout<<"           [-l time lengthh] "<<std::endl;
  std::cout<<"           [-psd] or [-pow] or [-spl]"<<std::endl;

  std::cout<<"\n       If DFT is preferred just add to your arguments the flag [-dft]"<<std::endl;

  std::cout<<"\nfor a detailed input list and discription use [--help]\n";

  //std::cout<<"\nAlternatively, you can just parse an input parsing file as follows: fftpsd inputfile\n";
  exit(0);
}


