#include"fft_param.h"

namespace fftspc {

void Case_Param::Parse(const std::string &in_fname){
  GetPot gp_input(in_fname.c_str());

  data_fname=gp_input("input_file","NO_INPUT_FILE",false);
  input_dir=GetFileDirectory(data_fname);

  std::string temp=remove_extension(data_fname);
  std::size_t last_slash=temp.find_last_of("/");
  if (last_slash == std::string::npos) // meaning the name is absolute without any directory, i.e., working in current directory
    output_data_name=temp;
  else
    output_data_name=temp.substr(last_slash+1,temp.size());

  output_dir=gp_input("output_dir",input_dir.c_str(),false);
  if(output_dir.find_last_of("/") != output_dir.size()-1) // correct if it does not include "/" at the end
    output_dir+="/";

  data_row=gp_input("data_row",0,false);
  data_col=gp_input("data_col",1,false);

  if(data_row>0) data_row--;
  if(data_col>1) data_col--;

  if(data_fname=="NO_INPUT_FILE"){
    _print_log("Reading the data of an analytical wave");
    wave_param.type=string_to_enum<Wave_Form_Type>(gp_input("wave/wave_form","sine",false));
    wave_param.Parse(in_fname);
  }
  psd_flag=gp_input("psd_flag",0,false);
  spl_flag=gp_input("spl_flag",0,false);

  if(psd_flag || spl_flag){
    fft_param.Parse_psd_defaults(in_fname);
    fft_param.avgfft_mode==VARIANCE; // force it to use this mode for both types of PSD/POWER and SPL
    if(psd_flag){
      psd_param.Parse(in_fname);
      psd_param.scaling_factor=FFT<double>::GetWindowScalingFactor(fft_param.window_type,psd_param.type_);
    }
  }else{
    fft_param.Parse(in_fname);
  }

  return;
}

void Case_Param::Parse(int argc, char **argv){
  // right now this works only with an input data file and not for analytical waves
  GetPot cmdline(argc,argv);

  data_fname=cmdline.follow("NO_INPUT_FILE",2,"-i","--file");
  if(data_fname=="NO_INPUT_FILE"){
    FatalError("There is no data file");
    help_short();
    exit(0);
  }

  input_dir=GetFileDirectory(data_fname);
  std::string temp=remove_extension(data_fname);
  std::size_t last_slash=temp.find_last_of("/");
  if (last_slash == std::string::npos) // meaning the name is absolute without any directory, i.e., working in current directory
    output_data_name=temp;
  else
    output_data_name=temp.substr(last_slash+1,temp.size());

  output_dir=cmdline.follow(input_dir.c_str(),2,"-o","--outdir");
  if(output_dir.find_last_of("/") != output_dir.size()-1) // correct if -o does not include "/" at the end
    output_dir+="/";

  data_row=cmdline.follow(0,"-r");
  data_col=cmdline.follow(1,"-c");

  if(data_row>0) data_row--;
  if(data_col>1) data_col--;

  psd_flag=0;
  if(cmdline.search("-psd")){
    psd_flag=1;
    psd_param.type_=DENSITY;
    std::vector<std::string> nominus_vec=cmdline.nominus_followers("-psd");
    if(nominus_vec.size()>0)
      psd_output_file=nominus_vec[0];
    else
      psd_output_file="DEFAULT";
  }
  if(cmdline.search("-pow")){
    psd_flag=1;
    psd_param.type_=POWER;
    std::vector<std::string> nominus_vec=cmdline.nominus_followers("-pow");
    if(nominus_vec.size()>0)
      psd_output_file=nominus_vec[0];
    else
      psd_output_file="DEFAULT";
  }
  if(cmdline.search("-spl")){
    spl_flag=1;
    std::vector<std::string> nominus_vec=cmdline.nominus_followers("-spl");
    if(nominus_vec.size()>0)
      spl_output_file=nominus_vec[0];
    else
      spl_output_file="DEFAULT";
  }

  if(psd_flag || spl_flag){
    fft_param.Parse_psd_defaults(cmdline);
    fft_param.avgfft_mode==VARIANCE; // force it to use this mode for both types of PSD/POWER and SPL
    if(psd_flag){
      psd_param.Parse(cmdline);
      psd_param.scaling_factor=FFT<double>::GetWindowScalingFactor(fft_param.window_type,psd_param.type_);
    }
  }else{
    fft_param.Parse(cmdline);
  }

  if(cmdline.search("-spec"))
    noheader=0;

  return;

}

void Wave::GetWaveSignal(std::vector<double> &time_vec, std::vector<double> &data_vec){
  time_vec.resize(Nnodes,0.);
  for(int i=0; i<Nnodes; i++)
    time_vec[i]=t_start+i*dt;
  data_vec.resize(Nnodes,0.);
  if(type==SINE)
    for(int i=0; i<Nnodes; i++)
      data_vec[i]=twave_amp*sin(2*PI*twave_freq*time_vec[i]+twave_phase_shift)+twave_const;
  else if(type==COSINE)
    for(int i=0; i<Nnodes; i++)
      data_vec[i]=twave_amp*cos(2*PI*twave_freq*time_vec[i]+twave_phase_shift)+twave_const;
  else if(type==GAUSSIAN)
    for(int i=0; i<Nnodes; i++)
      data_vec[i]=gwave_amp*exp(-gwave_factor*time_vec[i]*time_vec[i]);

  return;
}

void Wave::DumpWaveSignal(std::string fname_in,const std::vector<double> time_vec,std::vector<double> data_vec){
  std::ofstream fout;
  fout.open(fname_in);
  if( !fout.is_open() ) // ...else, create new file...
    fout.open(fname_in, std::ios_base::out | std::ios_base::trunc);
  fout<<"#Nnodes="<<Nnodes;
  fout<<",\tL="<<*(time_vec.end()-1)<<",\tT="<<T;
  fout<<",\tsample_freq="<<sample_freq<<std::endl;
  fout<<"#time, u"<<std::endl;
  fout<<std::setprecision(20);
  for(int itt=0; itt<Nnodes-1; itt++)
    fout<<time_vec[itt]<<" "<<data_vec[itt]<<std::endl;
  return;
}

void Wave::Parse(const std::string &in_fname){
  GetPot input(in_fname.c_str());
  t_start=input("wave/t_start",0.,false);
  t_final=input("wave/t_end",-1e-5,false);
  if(t_final<1e-5){ // to make sure it is not specified
    N_periods=input("wave/N_periods",1.,false);
  }else{
    N_periods=-1e-5;
  }
  Nnodes=input("wave/N_elem",64,false)+1;

  if(type==SINE || type==COSINE){
    twave_freq=input("wave/trigonometric/freq",1.,false);
    twave_amp=input("wave/trigonometric/amp",1.,false);
    twave_const=input("wave/trigonometric/const",1.,false);
    twave_phase_shift=input("wave/trigonometric/phase_shift",1.,false);
    if(t_final<1e-5){
      T=1./twave_freq;
      t_final=t_start+N_periods*T;
    }else{
      N_periods=(t_final-t_start)/T;
    }
  }else if(type==GAUSSIAN){
    gwave_amp=input("wave/Gaussian/amp",1.,false);
    gwave_factor=input("wave/Gaussian/factor",1.,false);
  }

  dt=(t_final-t_start)/(Nnodes-1);
  sample_freq=1./dt;

  std::cout<<"Wave: type="<<enum_to_string<Wave_Form_Type>(type)<<std::endl;
  if(type==SINE ||  type==COSINE){
    std::cout<<"Wave: amp="<<twave_amp<<std::endl;
    std::cout<<"Wave: freq="<<twave_freq<<std::endl;
    std::cout<<"Wave: phase_shift="<<twave_phase_shift<<std::endl;
    std::cout<<"Wave: const="<<twave_const<<std::endl;
  }else if(type==GAUSSIAN){
    std::cout<<"Wave: amp="<<gwave_amp<<std::endl;
    std::cout<<"Wave: factor="<<gwave_factor<<std::endl;
  }
  std::cout<<"Wave: T_period="<<T<<std::endl;
  std::cout<<"Wave: N_periods="<<N_periods<<std::endl;
  std::cout<<"Wave: sample_freq="<<sample_freq<<std::endl;
  std::cout<<"Wave: dt="<<dt<<std::endl;
  std::cout<<"Wave: t_start="<<t_start<<std::endl;
  std::cout<<"Wave: t_final="<<t_final<<std::endl;
  std::cout<<"Wave: N_nodes="<<Nnodes<<std::endl;

  return;
}


FFT_Param& FFT_Param::operator =(const FFT_Param& Rdata_in){
  Nt=Rdata_in.Nt;
  Lt=Rdata_in.Lt;
  dt=Rdata_in.dt;
  aver_flag=Rdata_in.aver_flag;
  Navg=Rdata_in.Navg;
  Lt_sub=Rdata_in.Lt_sub;
  Nt_sub=Rdata_in.Nt_sub;
  dt_sub=Rdata_in.dt_sub;
  shift=Rdata_in.shift;
  Nt_shifted=Rdata_in.Nt_shifted;
  aver_count=Rdata_in.aver_count;
  window_type=Rdata_in.window_type;
  Wwind=Rdata_in.Wwind;
  avgfft_mode=Rdata_in.avgfft_mode;
  wind_scaling=Rdata_in.wind_scaling;
  mean_substract=Rdata_in.mean_substract;
  return *this;
}

void FFT_Param::SetupFFTData(const double sample_dt){
  // note that sample_dt cannot be zero, otherwise this function cannot work
  if(Lt_sub<=sample_dt || Lt_sub<=1.e-10){
    FatalErrorST("Fatal error, fft_param window length is too short < 1e-10 or less than dt");
  }else if(sample_dt<=1.e-11){
    FatalErrorST("Fatal error, global dt is < 1e-11");
  }

  if(dt_sub>sample_dt && fabs(dt_sub-sample_dt)>=1e-2*sample_dt){ // dt is specified as an input and not equal to the simulation dt
    // first check if dt_sub is an integer multiple of sample_dt
    int temp_int=std::round(dt_sub/sample_dt);
    dt_sub = temp_int*sample_dt;
    // second check if Lt_sub is an integer number of dt_sub
    int N_temp=std::round(Lt_sub/dt_sub)+1; // this ensures we are not less than Lt_sub
    Lt_sub=(N_temp-1)*dt_sub; // at this point N_temp is integer and Lt_sub is adjusted

    // determine the nearest power of 2 using the adjusted dt_sub and window length
    int exponent=std::log2(N_temp);
    double diff0 = abs(N_temp-std::pow(2,exponent));
    if(diff0!=0){ // N_temp is not a power of two
      double diff1=abs(N_temp-std::pow(2,exponent+1));
      if(diff0<diff1)
        N_temp=std::pow(2,exponent);
      else
        N_temp=std::pow(2,exponent+1);
    }
    Nt_sub=N_temp;
    Lt_sub=(Nt_sub-1)*dt_sub;

  }else{ // dt is either unspecified or very close to the simulation dt
    dt_sub=sample_dt;  // let it be the dt_sim at the beginning
    // first check if Lt_sub is an integer number of sample_dt
    int N_temp=std::round(Lt_sub/sample_dt)+1; // this ensures we are not less than Lt_sub
    Lt_sub=(N_temp-1)*sample_dt; // adjust to the new Lt_sub

    if(Nt_sub>=2) // making sure it is a valid specified input
      N_temp=Nt_sub;
    // determine the nearest power of 2 using the simulation dt and the specified window length
    int exponent=std::log2(N_temp);
    double diff0 = abs(N_temp-std::pow(2,exponent));
    if(diff0!=0){ // N_temp is not a power of two
      double diff1=abs(N_temp-std::pow(2,exponent+1));
      if(diff0<diff1)
        N_temp=std::pow(2,exponent);
      else
        N_temp=std::pow(2,exponent+1);
    }
    double Lt0=(N_temp-1)*sample_dt;
    int ratio=std::round(Lt_sub/Lt0);
    if(ratio%2==0){
      Lt_sub=ratio*Lt0;
      if(ratio>=1)
        dt_sub=Lt_sub/(N_temp-1);
      else
        dt_sub=sample_dt;
      Nt_sub=N_temp;
      Lt_sub=(Nt_sub-1)*dt_sub;
    }else{ // they are not a power of two apart
      if(ratio<=1){
        Lt_sub=Lt0;
        Nt_sub=N_temp;
        dt_sub=Lt0/(N_temp-1);
      }else{
        Lt_sub=2*Lt0;
        dt_sub=Lt_sub/(N_temp-1); // this settings keeps the dt the same as in Lt0
        Nt_sub=N_temp;
      }
    }
  }

  if(shift>5.e-2){ // anything less than 5% is no averaging or overlapping
    aver_count=0;
    aver_flag=true;
    if(shift+1.e-5<1.)
      Nt_shifted=std::round(Nt_sub*shift);
    else // shift==1 meaning complete overlap and ensemble average
      Nt_shifted=Nt_sub;

    Navg=(Nt-Nt_sub)/Nt_shifted;
    aver_count=Navg+1;

  }else{ //shift<=5.e-2;
    aver_flag=false;
    Nt_shifted=0;
    Lt=Lt_sub;
    Nt=Nt_sub;
    dt=dt_sub;
    Navg=0;
  }
  //note: shift should be thought of as the overlap distance of two Lt_sub

  if(window_type==HANN){
    FFT<double>::Hanning(Nt_sub,"Periodic",avgfft_mode,Wwind,wind_scaling);
  }else if(window_type==HAMM){
    FFT<double>::Hamming(Nt_sub,"Periodic",avgfft_mode,Wwind,wind_scaling);
  }else if(window_type==TRIANGULAR){
    FFT<double>::Triangular(Nt_sub,Nt_sub,avgfft_mode,Wwind,wind_scaling);
  }else if(window_type==BARTLETT){
    FFT<double>::Bartlett(Nt_sub,"Periodic",avgfft_mode,Wwind,wind_scaling);
  }else if(window_type==WELCH){
    FFT<double>::Welch(Nt_sub,"Periodic",avgfft_mode,Wwind,wind_scaling);
  }else if(window_type==BLACKMAN){
    FFT<double>::Blackman(Nt_sub,"Periodic",avgfft_mode,Wwind,wind_scaling);
  }else if(window_type==RECTANGULAR){  // rectangular or no windowing
    Wwind.resize(Nt_sub);
    std::fill_n(Wwind.begin(),Nt_sub,1.);
    wind_scaling=1.;
  }

  return;
}

void help_short(){
  std::cout<<"\nhelp with implemented options"<<std::endl;
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<"       minimum arguments for fft computations"<<std::endl;
  std::cout<<"       --------------------------------------"<<std::endl;
  std::cout<<"           [-i datafile] "<<std::endl;
  std::cout<<"           [-n N] "<<std::endl;

  std::cout<<"\n       minimum arguments for psd, power & spl computations"<<std::endl;
  std::cout<<"       -----------------------------------------------------"<<std::endl;
  std::cout<<"           [-i datafile]"<<std::endl;
  std::cout<<"           [-n N] "<<std::endl;
  std::cout<<"           [-psd] or [-pow] or [-spl]"<<std::endl;

  std::cout<<"\nfor a detailed input list and discription use [--help]\n";

  std::cout<<"\nAlternatively, you can just parse an input parsing file as follows: fftpsd inputfile\n";
  exit(0);
}


}
