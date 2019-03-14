#ifndef FFT_PARAM_H
#define FFT_PARAM_H

#include"../include/getpot.h"
#include"general_tools.h"
#include"fft_type.h"
#include"string_to_type.h"
#include"fft.h"

namespace fftspc{

struct PSD_Param{
  int bin_aver_flag=0;
  double octave_n=1.;
  double p_ref=2e-5; //Pa
  int mode=1;  //1:PSD only, 2:PSD/SPL/OASPL
  PSD_Type type_=DENSITY;
  double scaling_factor=1.;

  void Parse(const std::string &in_fname){
    GetPot input(in_fname.c_str());
    bin_aver_flag=input("psd/bin_aver_flag",0,false);
    octave_n=input("psd/octave_n",1.,false);
    p_ref=input("psd/p_ref",2.e-5,false);
    type_=string_to_enum<PSD_Type>(input("psd/type","density",false));
    return;
  }

  void Parse(GetPot& cmdline){
    if(cmdline.search("-bin"))
      bin_aver_flag=1;
    octave_n=cmdline.follow(1.,"-oct");
    if(cmdline.search("-psd"))
      type_=DENSITY;
    if(cmdline.search("-pow"))
      type_=POWER;
    p_ref=cmdline.follow(2.e-5,"-pref");
  }
};

struct FFT_Param{
  uint Nt;             // no. of points for the whole simulation
  uint Nt_sub;         // no. of points in each subset of Nt, for averaging
  int  Navg=0;         // total number of averages
  double dt_sub=0.;       // physical dt of the subset
  double dt=0.;       // physical dt of the whole simulation time
  double Lt=0.;           // width/period of the whole simulation time
  double Lt_sub=0.;    // width/period of the subset time signal, or T_period
  int Nt_shifted=0;    // the shift position
  double shift=0.5;    // the shift distance between different subsets
  int aver_count=-1;   // counter for the number of averages
  bool aver_flag=false;// a flag for averaging
  bool mean_substract=true; // a flag for substracting the mean
  FFT_WINDOW_Type window_type=RECTANGULAR;   // default "Hann" windowing, no windowing is "NO"
  AVGFFT_Mode_Type avgfft_mode=VARIANCE; // default: "variance" to preserve the variance and energy of the signal
  std::vector<double> Wwind;   // array of the windowing values
  double wind_scaling=1.;      // for windowing scaling in consistence with avgfft_mode

  virtual FFT_Param& operator =(const FFT_Param& Rdata_in);
  void SetupFFTData(const double sample_dt);

  void Parse(const std::string &in_fname){
    GetPot input(in_fname.c_str());
    Nt_sub=input("fft/N",0,false);
    Lt_sub=input("fft/window_length",0.,false);
    double dt_temp=0.;
    if(Nt_sub>0)
        dt_temp=Lt_sub/(Nt_sub-1);  // a default value
    dt_sub=input("fft/dt",dt_temp,false);
    if(Lt_sub<=1.e-5)
      Lt_sub=dt_sub*Nt_sub;
    shift=input("fft/shift",0.,false); //default no shifting or windowing
    window_type=string_to_enum<FFT_WINDOW_Type>(input("fft/window","RECTANGULAR",false));
    avgfft_mode=string_to_enum<AVGFFT_Mode_Type>(input("fft/magnitude_scaling_mode","peak",false));
    int mean_subs_flag;
    mean_subs_flag=input("fft/mean_substract",1,false); // default is mean substract
    if(mean_subs_flag==1)
      mean_substract=true;
    else
      mean_substract=false;

    return;
  }

  void Parse_psd_defaults(const std::string &in_fname){
    // defaults if the psd flag is turned on
    // hann window with 50% shift and variance preserving
    GetPot input(in_fname.c_str());
    Nt_sub=input("fft/N",0,false);
    Lt_sub=input("fft/window_length",0.,false);
    double dt_temp=0.;
    if(Nt_sub>0)
        dt_temp=Lt_sub/(Nt_sub-1);  // a default value
    dt_sub=input("fft/dt",dt_temp,false);
    if(Lt_sub<=1.e-5)
      Lt_sub=dt_sub*Nt_sub;
    shift=input("fft/shift",0.5,false); //default 0.5 shifting
    window_type=string_to_enum<FFT_WINDOW_Type>(input("fft/window","HANN",false));
    avgfft_mode=VARIANCE; // force to this value for psd and spl computation
    int mean_subs_flag;
    mean_subs_flag=input("fft/mean_substract",1,false); // default is mean substract
    if(mean_subs_flag==1)
      mean_substract=true;
    else
      mean_substract=false;

    return;
  }

  void Parse(GetPot& cmdline){
    Nt_sub=cmdline.follow(0,"-n");
    Lt_sub=cmdline.follow(0.,"-l");
    double dt_temp=0.;
    if(Nt_sub>0)
        dt_temp=Lt_sub/(Nt_sub-1);  // a default value
    dt_sub=cmdline.follow(dt_temp,"-dt");
    if(Lt_sub<=1.e-5)
      Lt_sub=dt_sub*Nt_sub;
    shift=cmdline.follow(0.,"-s"); //default no shifting or windowing
    window_type=string_to_enum<FFT_WINDOW_Type>(cmdline.follow("RECTANGULAR","-w"));
    avgfft_mode=PEAK; // default for fft averaging
    if(cmdline.search("-variance"))
      avgfft_mode=VARIANCE;

    int mean_subs_flag=cmdline.follow(1,"-m"); // default is mean substract
    if(mean_subs_flag==1)
      mean_substract=true;
    else
      mean_substract=false;
    return;
  }

  void Parse_psd_defaults(GetPot& cmdline){
    Nt_sub=cmdline.follow(0,"-n");
    Lt_sub=cmdline.follow(0.,"-l");
    double dt_temp=0.;
    if(Nt_sub>0)
        dt_temp=Lt_sub/(Nt_sub-1);  // a default value
    dt_sub=cmdline.follow(dt_temp,"-dt");
    if(Lt_sub<=1.e-5)
      Lt_sub=dt_sub*Nt_sub;
    shift=cmdline.follow(0.5,"-s"); //default 50% overlap
    window_type=string_to_enum<FFT_WINDOW_Type>(cmdline.follow("HANN","-w"));
    avgfft_mode=VARIANCE; // force to this value for psd and spl computation

    int mean_subs_flag=cmdline.follow(1,"-m"); // default is mean substract
    if(mean_subs_flag==1)
      mean_substract=true;
    else
      mean_substract=false;
    return;
  }
};

struct Wave{
  const double PI=3.1415926535897932384626433832795029L;

  Wave_Form_Type type=SINE;
  double t_start=0.;
  double t_final=-1e-5;
  double T=0.;
  double N_periods=-1e-5;
  int Nnodes=0;
  double dt=0.;
  double sample_freq=0.;

  // Trigonometric wave data
  double twave_freq=0.;
  double twave_amp=0.;
  double twave_const=0.;
  double twave_phase_shift=0.;

  //Gaussian wave data
  double gwave_factor=0.;
  double gwave_amp=0.;

  void Parse(const std::string &in_fname);
  void GetWaveSignal(std::vector<double> &time_vec, std::vector<double> &data_vec);
  void DumpWaveSignal(std::string fname_in,const std::vector<double> time_vec,std::vector<double> data_vec);
};

struct Case_Param{
  // Case parameters:
  //-----------------------------
  std::string data_fname="NO_INPUT_FILE";    // parent folder for the case results
  std::string output_data_name="NO_INPUT_FILE";
  std::string input_dir;
  std::string output_dir;
  unsigned int data_row=0;  // row in the file to extract data
  unsigned int data_col=1;  // col in the file to extract data
  unsigned int psd_flag=0;
  unsigned int spl_flag=0;

  std::string psd_output_file="DEFAULT";
  std::string spl_output_file="DEFAULT";
  unsigned int noheader=1;

  FFT_Param fft_param;
  PSD_Param psd_param;
  Wave wave_param;
  // Function members:
  //--------------------------
  void Parse(const std::string &fname);
  void Parse(int argc, char** argv);

//  void dump_python_inputfile();

//  void Reset();
};

void help_short();


}

#endif
