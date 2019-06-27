#include"fft_driver.h"

namespace fftspc {

void FFT_Driver::Init(){

  fft_param_=case_param_.fft_param;
  if(case_param_.psd_flag)
    psd_param_=case_param_.psd_param;

  if(case_param_.data_fname=="NO_INPUT_FILE"){
    std::cout<<"output_dir=" <<case_param_.output_dir<<std::endl;
    fft_param_.Nt=case_param_.wave_param.Nnodes-1;
    fft_param_.dt=case_param_.wave_param.dt;
    fft_param_.Lt=case_param_.wave_param.t_final-case_param_.wave_param.t_start-fft_param_.dt;
    case_param_.wave_param.GetWaveSignal(time_vec_,udata_vec_);
    std::string ofname=case_param_.output_dir+std::string("u_data.dat");
    case_param_.wave_param.DumpWaveSignal(ofname,time_vec_,udata_vec_);
    // removing the last periodic element for fft to work
    time_vec_.erase(time_vec_.end()-1);
    udata_vec_.erase(udata_vec_.end()-1);

  }else{ // use input file data to setup the fft Lt,Nt,dt data
    std::cout<<"input_dir ="<<case_param_.input_dir<<std::endl;
    std::cout<<"output_dir=" <<case_param_.output_dir<<std::endl;
    ReadTimeData(case_param_.data_fname);
    fft_param_.Nt=time_vec_.size();
    fft_param_.Lt=time_vec_[fft_param_.Nt-1]-time_vec_[0];
    fft_param_.dt=time_vec_[1]-time_vec_[0];
    if(fft_param_.Nt_sub>0)
      fft_param_.Lt_sub=fft_param_.dt*(fft_param_.Nt_sub-1);
  }

  if(fft_param_.DFT_mode){
    fft_param_.SetupDFTData(fft_param_.dt);
    if(fft_param_.Nt_sub>10000){
      std::cout<<"WARNING: the number of data points for DFT exceeded 10000, this will take some time to compute";
      std::cout<<"\n         please consider using an FFT algorithm instead, i.e., N=2^{k}, k is integer"<<std::endl;
    }
  }else{
    fft_param_.SetupFFTData(fft_param_.dt);
  }

  if(fft_param_.Nt_sub>(case_param_.data_lastrow-case_param_.data_row)){
    std::string str_message="The input data range (r:q,c) is less than the required N="
                             +std::to_string(fft_param_.Nt_sub)+" number of data points ";
    FatalErrorST(str_message.c_str());
  }

  if(fft_param_.DFT_mode){
    fft_= new FFT<double>(fft_param_.Nt_sub,"DFT");
  }else{
    fft_= new FFT<double>(fft_param_.Nt_sub,"real");
  }

  std::cout<<"FFT: N data actually used     ="<<fft_param_.Nt_sub+fft_param_.Navg*fft_param_.Nt_shifted<<std::endl;
  std::cout<<"FFT: N data in a window subset="<<fft_param_.Nt_sub<<std::endl;
  std::cout<<"FFT: N data shifted           ="<<fft_param_.Nt_shifted<<std::endl;
  std::cout<<"FFT: N of window subsets      ="<<fft_param_.Navg+1<<std::endl;
  std::cout<<"FFT: L,     sample length ="<<fft_param_.Lt<<std::endl;
  std::cout<<"FFT: L_sub, window length ="<<fft_param_.Lt_sub<<std::endl;
  std::cout<<"FFT: dt,    in the sample ="<<fft_param_.dt<<std::endl;
  std::cout<<"FFT: dt_sub, actually used="<<fft_param_.dt_sub<<std::endl;
  std::cout<<"FFT: window function      =\""<<enum_to_string<FFT_WINDOW_Type>(fft_param_.window_type)<<"\""<<std::endl;
  std::cout<<"FFT: window_scaling_factor="<<fft_param_.wind_scaling<<std::endl;
  if(case_param_.psd_flag>0)
    std::cout<<"PSD: psd type = \""<<enum_to_string<PSD_Type>(psd_param_.type_)<<"\""<<std::endl;


  //Setup the frequency vector
  fft_->rfftfreq(fft_param_.Nt_sub,fft_param_.dt_sub,freq_);
  Nfreq_=freq_.size();  // frequency size
  omega_.resize(Nfreq_,0.); // angular frequency
  for(int iff=0; iff<Nfreq_; iff++)
    omega_[iff]=2.*PI*freq_[iff];

  return;
}

void FFT_Driver::ComputeFFT(){
  vector<double> u_sub(fft_param_.Nt_sub,0.); // subset data vector
  std::copy_n(udata_vec_.begin(),fft_param_.Nt_sub,u_sub.begin());
  if(fft_param_.window_type!=RECTANGULAR && fft_param_.mean_substract) // if windowing substract the mean first
    SubstractMean(u_sub);

  for(int itt=0; itt<fft_param_.Nt_sub; itt++)
    u_sub[itt]*=fft_param_.Wwind[itt];

  if(fft_param_.DFT_mode)
    fft_->dft(u_sub,fft_vec_);
  else
    fft_->rfft(u_sub,fft_vec_);

  Convert2MagPhaseScaled();

  return;
}

void FFT_Driver::ComputeFFTavg(){

  vector<double> u_sub(fft_param_.Nt_sub,0.); // subset data vector
  fft_vec_.resize(Nfreq_,I0);
  vector<double> mag_sum(Nfreq_,0.),fft_sum_real(Nfreq_,0.),fft_sum_imag(Nfreq_,0.);

  int local_fft_iter=fft_param_.dt_sub/fft_param_.dt;

  double u_mean = accumulate( udata_vec_.begin(), udata_vec_.end(), 0.0)/udata_vec_.size();
  std::cout<<"data mean="<<u_mean<<std::endl;

  if(fft_param_.window_type!=RECTANGULAR && fft_param_.mean_substract) // if windowing substract the mean first
    SubstractMean(udata_vec_);

  for(int i_avg=0; i_avg<fft_param_.aver_count; i_avg++){ // averaging loop
    for(int itt=0; itt<fft_param_.Nt_sub; itt++){
      int global_iter=local_fft_iter*itt+i_avg*fft_param_.Nt_shifted;
      u_sub[itt]=udata_vec_[global_iter]*fft_param_.Wwind[itt];
    }
//    if(fft_param_.window_type!=RECTANGULAR && fft_param_.mean_substract) // if windowing substract the mean first
//      SubstractMean(u_sub);

    if(fft_param_.DFT_mode)
      fft_->dft(u_sub,fft_vec_);
    else
      fft_->rfft(u_sub,fft_vec_);
    AccumulateFFT(fft_param_.avgfft_mode,fft_param_.wind_scaling,fft_vec_
                  ,mag_sum,fft_sum_real,fft_sum_imag);
  }

  for(int iff=0; iff<Nfreq_; iff++)
    AverageFFTOneFreq(fft_param_.avgfft_mode,fft_param_.aver_count
                      ,mag_sum[iff],fft_sum_real[iff],fft_sum_imag[iff],fft_vec_[iff]);

  Convert2MagPhaseScaled();

  return;
}

void FFT_Driver::Convert2MagPhaseScaled(){

  fft_mag_.resize(Nfreq_,0.);
  fft_phase_.resize(Nfreq_,0.);
  for(int iff=1; iff<Nfreq_-1; iff++){ // fftreal scaling
    fft_vec_[iff]*= (fft_param_.wind_scaling/((double)fft_param_.Nt_sub));
    fft_mag_[iff]=std::abs(2.*fft_vec_[iff]);
    fft_phase_[iff]=std::arg(2.*fft_vec_[iff]);
  }
  int iff=0;
  fft_vec_[iff]*= (fft_param_.wind_scaling/((double)fft_param_.Nt_sub));
  fft_mag_[iff]=std::abs(fft_vec_[iff]);
  fft_phase_[iff]=std::arg(fft_vec_[iff]);
  iff=Nfreq_-1;
  fft_vec_[iff]*= (fft_param_.wind_scaling/((double)fft_param_.Nt_sub));
  fft_mag_[iff]=std::abs(fft_vec_[iff]);
  fft_phase_[iff]=std::arg(fft_vec_[iff]);

  return;
}

void FFT_Driver::ComputePSD(){

  double psd_scaling=1.0;  // since real fft was normalized by Nt_sub/2
  if(psd_param_.type_==DENSITY)
    psd_scaling=1./(freq_[1]);
  else if(psd_param_.type_==POWER){ // needs to make it general for any window functions
    psd_scaling=psd_param_.scaling_factor*psd_param_.scaling_factor
                    /(fft_param_.wind_scaling*fft_param_.wind_scaling);
  }

  if(!psd_param_.bin_aver_flag){
    psd_vec_.resize(Nfreq_,0.);
    psd_vec_[0]=std::pow(std::abs(fft_vec_[0]),2)*psd_scaling;
    for(int iff=1; iff<Nfreq_-1; iff++)
      psd_vec_[iff]=2.*std::pow(std::abs(fft_vec_[iff]),2)*psd_scaling;
    psd_vec_[Nfreq_-1]=std::pow(std::abs(fft_vec_[Nfreq_-1]),2)*psd_scaling;

    Nbins_=Nfreq_;
  }else{

  }

  return;
}

void FFT_Driver::ComputeSPL(){

  double p_abs_sum=0.,p_abs;
  double sqroot2=std::sqrt(2.);
  double spl_scaling=1.0;

  //if(!psd_param_.bin_aver_flag){
    spl_vec_.resize(Nfreq_,0.);

    p_abs=std::abs(fft_vec_[0])*spl_scaling;
    spl_vec_[0]=20.*std::log10(p_abs/psd_param_.p_ref);
    p_abs_sum=p_abs;

    for(int iff=1; iff<Nfreq_-1; iff++){
      p_abs=sqroot2*std::abs(fft_vec_[iff])*spl_scaling;
      spl_vec_[iff]=20.*std::log10(p_abs/psd_param_.p_ref);
      p_abs_sum+=p_abs;
    }

    p_abs=std::fabs(fft_vec_[Nfreq_-1])*spl_scaling;
    spl_vec_[Nfreq_-1]=20.*std::log10(p_abs/psd_param_.p_ref);
    p_abs_sum+=p_abs;
    oaspl_=20.*std::log10(p_abs_sum/psd_param_.p_ref);

  //}else{

  //}

  return;
}

void FFT_Driver::ReadTimeData(const std::string in_fname){
  std::ifstream input;
  std::string local_fname=in_fname;
  open_inputfile_forreading(local_fname,input);

  if(case_param_.data_lastrow<0.99e7)
    std::cout<<"Reading data in the range (row,col) = ( "<<case_param_.data_row+1<<":"<<case_param_.data_lastrow+1<<","
             <<case_param_.data_col+1<<" )"<<std::endl;
  else
    std::cout<<"Reading data in the range (row,col) = ( "<<case_param_.data_row+1<<":end,"
             <<case_param_.data_col+1<<" )"<<std::endl;


  unsigned long int i=0;
  std::string line;
  int n_numbers;

  while(std::getline(input,line) && i<=case_param_.data_lastrow ){
    if(i>=case_param_.data_row ){
      if(!is_a_comment_line(line) && !is_a_text_line(line)){ // make sure you skip text lines
        std::vector<double> line_data;
        line2doubledata(line,line_data,n_numbers);
        time_vec_.push_back(line_data[0]);
        udata_vec_.push_back(line_data[case_param_.data_col]);
      }
    }
    i++;
  }

  std::cout<<std::setprecision(16);
  std::cout<<"N data read="<<time_vec_.size()<<std::endl;
  std::cout<<"dt="<<time_vec_[1]-time_vec_[0]<<std::endl;
  std::cout<<"Fs="<<1./(time_vec_[1]-time_vec_[0])<<"\n\n";

  return;
}

void FFT_Driver::DumpOutputs(){

  if(case_param_.data_fname=="NO_INPUT_FILE"){
    std::cout<<"\nfreq vs fft_data:\n______________________"<<std::endl;
    for(int i=0; i<Nfreq_; i++)
      std::cout<<freq_[i]<<"     "<<fft_mag_[i]<<std::endl;
  }

  std::string fname_fft;
  if(fft_param_.DFT_mode)
    fname_fft=case_param_.output_dir+string("dft_")+case_param_.output_data_name+string(".dat");
  else
    fname_fft=case_param_.output_dir+string("fft_")+case_param_.output_data_name+string(".dat");
  std::cout<<"\nfft_fname="<<fname_fft<<std::endl;
  dump_fft_results(fname_fft);

  if(case_param_.psd_flag==1){
    std::string fname_psd;
    if(psd_param_.type_==DENSITY){
      if(case_param_.psd_output_file!="DEFAULT")
        fname_psd=case_param_.psd_output_file;
      else
        fname_psd=case_param_.output_dir+string("psd_")+case_param_.output_data_name+string(".dat");
    }else if(psd_param_.type_==POWER){
      if(case_param_.psd_output_file!="DEFAULT")
        fname_psd=case_param_.psd_output_file;
      else
        fname_psd=case_param_.output_dir+string("power_")+case_param_.output_data_name+string(".dat");
    }
    std::cout<<"psd_fname="<<fname_psd<<std::endl;
    dump_psd_results_nohead(fname_psd);
    if(!case_param_.noheader){
      fname_psd=remove_extension(fname_psd);
      fname_psd+=string("_spec.dat");
      dump_psd_results(fname_psd);
    }
  }

  if(case_param_.spl_flag==1){
    std::string fname_spl;
    if(case_param_.spl_output_file!="DEFAULT")
      fname_spl=case_param_.spl_output_file;
    else
      fname_spl=case_param_.output_dir+string("SPL_")+case_param_.output_data_name+string(".dat");
    std::cout<<"spl_fname="<<fname_spl<<std::endl;
    dump_spl_results_nohead(fname_spl);
    if(!case_param_.noheader){
      fname_spl=remove_extension(fname_spl);
      fname_spl=fname_spl+string("_spec.dat");
      dump_spl_results(fname_spl);
    }
  }

}

void FFT_Driver::AccumulateFFT(const AVGFFT_Mode_Type mag_scaling_mode,const double in_scaling
                               ,const vector<fftspc::Complex> q_fft_in,vector<double> &mag_sum
                               ,vector<double> &fft_sum_real,vector<double> &fft_sum_imag){

  fftspc::Complex temp_fft=I0;
  if(mag_scaling_mode==PEAK)        // for peak preserving and power spectrum
    for(int iff=0; iff<Nfreq_; iff++){
      temp_fft=q_fft_in[iff];
      mag_sum[iff]+=std::abs(temp_fft);
      fft_sum_real[iff]+=temp_fft.real();
      fft_sum_imag[iff]+=temp_fft.imag();
    }
  else if(mag_scaling_mode==VARIANCE)  // for variance preserving and power spectral density (psd)
    for(int iff=0; iff<Nfreq_; iff++){
      temp_fft=q_fft_in[iff];
      mag_sum[iff]+=std::pow(std::abs(temp_fft),2);
      fft_sum_real[iff]+=temp_fft.real();
      fft_sum_imag[iff]+=temp_fft.imag();
    }

  return;
}

void FFT_Driver::WindowScalingFFT(const double in_scaling
                                  ,vector<fftspc::Complex> &fft_in){
  int Ntemp=fft_in.size();
  for(int iff=0; iff<Ntemp; iff++)
    fft_in[iff]*=in_scaling;
  return;
}

void FFT_Driver::WindowScalingPSD(const double in_scaling
                                  ,vector<double> &psd_in){
  int Ntemp=psd_in.size();
  for(int iff=0; iff<Ntemp; iff++)
    psd_in[iff]*=in_scaling;
  return;
}

void FFT_Driver::AverageFFTOneFreq(const AVGFFT_Mode_Type mag_scaling_mode,const double in_aver_count
                                   ,const double mag_sum_in,const double fft_sum_real
                                   ,const double fft_sum_imag,fftspc::Complex& avgfft_out){
  double mag,phase;
  if(mag_scaling_mode==PEAK)   // for peak preserving
    mag= mag_sum_in/in_aver_count;
  else if(mag_scaling_mode==VARIANCE) // for variance preserving
    mag=std::sqrt(mag_sum_in/in_aver_count);

  phase=std::arg(fftspc::Complex(fft_sum_real,fft_sum_imag)); // fix me!: check the angle if it like atan2 or atan in MATLAB
  avgfft_out=mag*std::exp(I*phase);

  return;
}

void FFT_Driver::SubstractMean(vector<double> &u_data){
  double u_mean = accumulate( u_data.begin(), u_data.end(), 0.0)/u_data.size();
  for(int itt=0; itt<u_data.size(); itt++) // substracting the mean of the data
    u_data[itt]=u_data[itt]-u_mean;
  return;
}

void FFT_Driver::dump_fft_results(const string fname_in){

  std::ofstream fout;
  fout.open(fname_in.c_str());
  if( !fout.is_open() ) // ...else, create new file...
    fout.open(fname_in.c_str(), std::ios_base::out | std::ios_base::trunc);

  fout<<"# INFO: mohammadalhawwary@gmail.com, 01/03/2019"<<std::endl;
  fout<<"# Nt="<<fft_param_.Nt;
  fout<<",\t\tLt="<<fft_param_.Lt<<",\t\t\tdt="<<fft_param_.dt<<std::endl;
  fout<<"# Nt_sub="<<fft_param_.Nt_sub;
  fout<<",\tLt_sub="<<fft_param_.Lt_sub<<",\tdt_sub="<<fft_param_.dt_sub<<std::endl;
  fout<<"# Nt data actually used="<<fft_param_.Nt_sub+fft_param_.Navg*fft_param_.Nt_shifted<<std::endl;
  fout<<"# shift="<<fft_param_.shift*100<<"\%";
  fout<<",\twindow="<<enum_to_string<FFT_WINDOW_Type>(fft_param_.window_type);
  fout<<",\tN of windows="<<fft_param_.aver_count<<std::endl;
  fout<<"# Nfreq="<<freq_.size()<<",\tdf="<<freq_[1]<<",\tfs="<<1./fft_param_.dt_sub<<std::endl;
  fout<<"# freq(Hz), mag, phase"<<std::endl;

  fout<<std::setprecision(16);
  for(int iff=0; iff<Nfreq_; iff++)
    fout<<freq_[iff]<<"   "<<fft_mag_[iff]<<" "<<fft_phase_[iff]<<std::endl;


  return;
}

void FFT_Driver::dump_psd_results(const string in_fname){

  std::ofstream fout;
  fout.open(in_fname.c_str());
  if( !fout.is_open() ) // ...else, create new file...
    fout.open(in_fname.c_str(), std::ios_base::out | std::ios_base::trunc);

  fout<<"# INFO: mohammadalhawwary@gmail.com, 01/03/2019"<<std::endl;
  fout<<"# Ntotal="<<fft_param_.Nt<<",\tNt data used="<<fft_param_.Nt_sub+fft_param_.Navg*fft_param_.Nt_shifted<<std::endl;
  fout<<"# Nt_window_subset="<<fft_param_.Nt_sub<<",\tLt_sub="<<fft_param_.Lt_sub<<" sec,\tdt_sub="<<fft_param_.dt_sub<<std::endl;
  fout<<"# shift="<<fft_param_.shift*100<<"\%";
  fout<<",\twindow="<<enum_to_string<FFT_WINDOW_Type>(fft_param_.window_type);
  fout<<",\tN of windows="<<fft_param_.aver_count<<std::endl;
  if(!psd_param_.bin_aver_flag){
    fout<<"# Nfreq="<<Nfreq_<<",\tdf="<<freq_[1]<<",\tfs="<<1./fft_param_.dt_sub<<std::endl;
    if(psd_param_.type_==DENSITY)
      fout<<"# freq(Hz), PSD"<<std::endl;
    else if(psd_param_.type_==POWER)
      fout<<"# freq(Hz), Power"<<std::endl;

    fout<<std::setprecision(16);
    for(int iff=0; iff<Nfreq_; iff++)
      fout<<freq_[iff]<<"   "<<psd_vec_[iff]<<std::endl;
  }

  return;
}

void FFT_Driver::dump_spl_results(const string in_fname){

  std::ofstream fout;
  fout.open(in_fname.c_str());
  if( !fout.is_open() ) // ...else, create new file...
    fout.open(in_fname.c_str(), std::ios_base::out | std::ios_base::trunc);

  fout<<"# INFO: mohammadalhawwary@gmail.com, 01/03/2019"<<std::endl;
  fout<<"# Ntotal="<<fft_param_.Nt<<",\tNt data used="<<fft_param_.Nt_sub+fft_param_.Navg*fft_param_.Nt_shifted<<std::endl;
  fout<<"# Nt_window_subset="<<fft_param_.Nt_sub<<",\tLt_sub="<<fft_param_.Lt_sub<<" sec,\tdt_sub="<<fft_param_.dt_sub<<std::endl;
  fout<<"# shift="<<fft_param_.shift*100<<"\%";
  fout<<",\twindow="<<enum_to_string<FFT_WINDOW_Type>(fft_param_.window_type);
  fout<<",\tN of windows="<<fft_param_.aver_count<<std::endl;
  //if(!psd_param_.bin_aver_flag){
    fout<<"# Nfreq="<<Nfreq_<<",\tdf="<<freq_[1]<<",\tfs="<<1./fft_param_.dt_sub<<std::endl;
    fout<<"# OASPL(dB)="<<oaspl_<<std::endl;
    fout<<"# freq(Hz), SPL(dB)"<<std::endl;
    fout<<std::setprecision(16);
    for(int iff=0; iff<Nfreq_; iff++)
      fout<<freq_[iff]<<"   "<<spl_vec_[iff]<<std::endl;
  //}else{

  //}

  return;
}

void FFT_Driver::dump_fft_results_nohead(const string fname_in){

  std::ofstream fout;
  fout.open(fname_in.c_str());
  if( !fout.is_open() ) // ...else, create new file...
    fout.open(fname_in.c_str(), std::ios_base::out | std::ios_base::trunc);

  fout<<std::setprecision(16);
  for(int iff=0; iff<Nfreq_; iff++)
    fout<<freq_[iff]<<"   "<<fft_mag_[iff]<<" "<<fft_phase_[iff]<<std::endl;

  return;
}

void FFT_Driver::dump_psd_results_nohead(const string in_fname){

  std::ofstream fout;
  fout.open(in_fname.c_str());
  if( !fout.is_open() ) // ...else, create new file...
    fout.open(in_fname.c_str(), std::ios_base::out | std::ios_base::trunc);

  fout<<std::setprecision(16);
  for(int iff=0; iff<Nfreq_; iff++)
    fout<<freq_[iff]<<"   "<<psd_vec_[iff]<<std::endl;

  return;
}

void FFT_Driver::dump_spl_results_nohead(const string in_fname){

  std::ofstream fout;
  fout.open(in_fname.c_str());
  if( !fout.is_open() ) // ...else, create new file...
    fout.open(in_fname.c_str(), std::ios_base::out | std::ios_base::trunc);


  fout<<std::setprecision(16);
  for(int iff=0; iff<Nfreq_; iff++)
    fout<<freq_[iff]<<"   "<<spl_vec_[iff]<<std::endl;

  return;
}





}



