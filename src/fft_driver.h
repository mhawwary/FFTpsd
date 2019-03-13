#ifndef FFT_DRIVER_H
#define FFT_DRIVER_H

#include<complex>
#include<vector>
#include<iostream>
#include<string>
#include<cmath>
#include"../include/error.h"
#include<numeric>
#include<fstream>
#include"fft_param.h"
#include"fft.h"
#include"general_tools.h"

using namespace std;

namespace fftspc {

/*!
 * \brief The FFT_Driver class
 * \author Mohammad Alhawwary
 * \date 01/03/2019
 */
class FFT_Driver{

public:
  //  Construction Functions :
  FFT_Driver(std::string in_fname){
    case_param_.Parse(in_fname);

    return;
  }
  ~FFT_Driver(){}

  void Reset(){}

  void SetupDriver(){}
  void Init();

  void ReadTimeData(const std::string in_fname);

  void Compute(){

    if(fft_param_.aver_flag)
      ComputeFFTavg();
    else
      ComputeFFT();

    if(case_param_.psd_flag==1){
      ComputePSD();
    }else if(case_param_.psd_flag==2){
      ComputePSD();
      ComputeSPL();
    }

    return;
  }
  void DumpOutputs();
  void ComputeDFT(){}
  void ComputeDFTavg(){}
  void ComputeFFT();
  void ComputeFFTavg();
  void ComputePSD();
  void ComputeSPL();

protected:


  void AccumulateFFT(const AVGFFT_Mode_Type mag_scaling_mode, const double in_scaling
                     , const vector<fftspc::Complex> q_fft_in
                     , vector<double> &mag_sum
                     , vector<double> &fft_sum_real
                     , vector<double> &fft_sum_imag);

  void WindowScalingFFT(const double in_scaling,vector<fftspc::Complex> &fft_in);
  void WindowScalingPSD(const double in_scaling,vector<double> &psd_in);

  void AverageFFTOneFreq(const AVGFFT_Mode_Type mag_scaling_mode, const double in_aver_count
                         , const double mag_sum_in, const double fft_sum_real
                         , const double fft_sum_imag, fftspc::Complex& avgfft_out);

  void SubstractMean(vector<double> &u_data);

  void Convert2MagPhaseScaled();

  void dump_fft_results(const string fname);
  void dump_psd_results(const string fname);
  void dump_spl_results(const string in_fname);

protected:
  FFT<double> *fft_=nullptr;
  Case_Param case_param_;
  FFT_Param fft_param_;
  PSD_Param psd_param_;
  std::vector<double> freq_,omega_,fft_mag_,fft_phase_;
  std::vector<double> time_vec_,udata_vec_;
  std::vector<fftspc::Complex> fft_vec_;
  std::vector<double> psd_vec_,spl_vec_;
  double oaspl_=0.;
  int Nfreq_;
  int Nbins_;
};

}



#endif
