#ifndef FFT_H
#define FFT_H

//using namespace std;

#include<complex>
#include<vector>
#include<algorithm>  // for swap function
#include<iostream>
#include<string>
#include<cmath>
#include"fft_type.h"
#include"string_to_type.h"

namespace  fftspc {

  typedef unsigned int uint;
  typedef std::complex<double>  Complex;
  const double PI=3.1415926535897932384626433832795029L;

  const Complex I=Complex(0.,1.);
  const Complex I0=Complex(0.,0.);

  /**
   *\brief The FFT class
   *\author Mohammad Alhawwary
   *\date 01/12/2018
   */
  template<typename T>
  class FFT{  // FFT on complex data, outputs both positive&negative frequency data

  protected:
    typedef std::complex<T> complex;
    typedef std::vector<T>  vector;
    typedef std::vector<vector>  matrix;
    typedef std::vector<complex> cvector;
    typedef std::vector<cvector> cmatrix;
    typedef std::vector<int>     ivector;
    const T PI=3.1415926535897932384626433832795029L;
    const complex I=complex(0.,1.);   //1.i
    const complex I0=complex(0.,0.);  //0.i

  public:
    //  Construction Functions :
    FFT(void){}
    FFT(const int n_data_size, const std::string in_mode){
      if(in_mode=="real"){
        real_mode=true;
        Init_real(n_data_size);
      }else{
        Init(n_data_size);
      }
    }

    void set_data_size(const int n_data_size){
      if(real_mode)
        Init_real(n_data_size);
      else
        Init(n_data_size);
      return;
    }

    void set_real_mode(const std::string in_mode){
      if(in_mode=="real")
        real_mode=true;
      else
        real_mode=false;
      return;
    }

    virtual ~FFT(void);
    template<typename T1, typename T2>
    void   fft(const T1 u_data,T2 &u_fft);  //in:real/complex data, out:complex fft
    void  ifft(const cvector u_fft, cvector &u_data); //in:complex fft,  out:real data
    void  rfft(const vector u_data, cvector &u_fft);  //in:real data,    out:complex fft
    void irfft(const cvector u_fft, vector &u_data); //in:complex fft,  out:real data
    void fftfreq(const int local_n, const T sample_spacing, vector &freq_); // freq positive and negative
    void fftintfreq(const int local_n, ivector &freq_); // for integer frequencies
    void rfftfreq(const int local_n, const T sample_spacing, vector &freq_); // positive freq
    void rfftintfreq(const int local_n, ivector &freq_); // for positive integer freq
    template<typename TT>
    void fftshift(TT &u_fft);

    static void Hanning(const int n_size_in,const std::string in_wind_mode,const AVGFFT_Mode_Type in_scaling_mode
                        ,vector& Wwind,T& scaling_factor);
    static void Hamming(const int n_size_in,const std::string in_wind_mode,const AVGFFT_Mode_Type in_scaling_mode
                        ,vector& Wwind,T& scaling_factor);
    static void Blackman(const int n_size_in,const std::string in_wind_mode,const AVGFFT_Mode_Type in_scaling_mode
                         ,vector& Wwind,T& scaling_factor);
    static void Triangular(const int n_size_in,const int L,const AVGFFT_Mode_Type in_scaling_mode
                           ,vector& Wwind,T& scaling_factor);
    static void Bartlett(const int n_size_in,const std::string in_wind_mode,const AVGFFT_Mode_Type in_scaling_mode
                         ,vector& Wwind,T& scaling_factor);
    static void Welch(const int n_size_in,const std::string in_wind_mode,const AVGFFT_Mode_Type in_scaling_mode
                      ,vector& Wwind,T& scaling_factor);

    static T GetWindowScalingFactor(FFT_WINDOW_Type window_type_in,PSD_Type psd_type_in);

  protected:
    void Init(const int n_data_size);
    void Init_real(const int n_data_size);
    template<typename TT>
    void bit_reverse(TT &u_data, const int in_data_size); // bit reverse of 1D array
    void fft_inplace(cvector &u_fft, const int local_n, const int local_n_levels);
    void ifft_inplace(cvector &u_fft, const int local_n, const int local_n_levels);
    void packing_oddeven_real2complexdata(const vector u_data, cvector& h_packed, const int n_packed);
    void folding_realimag_complex2realdata(const cvector h_data, vector& u_data, const int n_packed);

  protected:
    bool real_mode=false;  // if true then we only compute the positive frequency part for a real input
    cmatrix W_arr;
    cmatrix Wi_arr;
    ivector m_arr;
    ivector m_2_arr;
    int n_data;          // n
    int n_data_levels;   // log2(n)

    // for real fft:
    int n_2_data;        // n/2
    int n_4_data;        // n/4
    int n_2_data_levels; // log2(n/2)
    cvector Wn,Wn2;

  };

  // Windowing functions:
  template<typename T>
  void FFT<T>::Hanning(const int n_size_in,const std::string in_wind_mode,const AVGFFT_Mode_Type in_scaling_mode
                       ,vector& Wwind,T& scaling_factor){
    //W.H. Press, B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling, “Numerical Recipes”
    //Oppenheim, Alan V., Ronald W. Schafer, and John R. Buck. Discrete-Time Signal Processing. Upper Saddle River, NJ: Prentice Hall, 1999, pp. 468–471
    const T PI=3.1415926535897932384626433832795029L;
    int n_window;
    Wwind.resize(n_size_in);
    if(in_wind_mode=="SYMMETRIC" || in_wind_mode=="Symmetric" || in_wind_mode=="symmetric")
      n_window=n_size_in-1;
    else if(in_wind_mode=="PERIODIC" ||in_wind_mode=="Periodic" || in_wind_mode=="periodic")
      n_window=n_size_in;
    else{
      std::cout<<"Fatal Error: windowing mode is wrong, use either \"periodic\" or \"symmetric\""<<std::endl;
      return;
    }

    for(int i=0; i<n_size_in; i++)
      Wwind[i]=0.5*(1. - cos(2.*PI*i/((double)n_window)));

    if(in_scaling_mode==PEAK) // peak magnitude preserving
      scaling_factor=2.0;
    else if(in_scaling_mode==VARIANCE) // variance preserving
      scaling_factor=sqrt(8./3.);
    else{
      std::cout<<"Warning: wrong scaling, use either \"v\" or \"p\""<<std::endl;
      scaling_factor=1000.;
    }

    return;
  }

  template<typename T>
  void FFT<T>::Hamming(const int n_size_in,const std::string in_wind_mode,const AVGFFT_Mode_Type in_scaling_mode
                       ,vector& Wwind,T& scaling_factor){
    //Blackman, R.B. and Tukey, J.W., (1958) The measurement of power spectra, Dover Publications, New York.
    const T PI=3.1415926535897932384626433832795029L;
    int n_window;
    Wwind.resize(n_size_in);
    if(in_wind_mode=="SYMMETRIC" || in_wind_mode=="Symmetric" || in_wind_mode=="symmetric")
      n_window=n_size_in-1;
    else if(in_wind_mode=="PERIODIC" ||in_wind_mode=="Periodic" || in_wind_mode=="periodic")
      n_window=n_size_in;
    else{
      std::cout<<"Fatal Error: windowing mode is wrong, use either \"periodic\" or \"symmetric\""<<std::endl;
      return;
    }

    for(int i=0; i<n_size_in; i++)
      Wwind[i]=0.54 - 0.46*cos(2.*PI*i/n_window);

    if(in_scaling_mode==PEAK) // peak magnitude preserving
      scaling_factor=50./27.;
    else if(in_scaling_mode==VARIANCE) // variance preserving
      scaling_factor=(50.*sqrt(3974.))/1987.;
    else{
      std::cout<<"Warning: wrong scaling, use either \"v\" or \"p\""<<std::endl;
      scaling_factor=1000.;
    }

    return;
  }

  template<typename T>
  void FFT<T>::Blackman(const int n_size_in,const std::string in_wind_mode,const AVGFFT_Mode_Type in_scaling_mode
                        ,vector& Wwind,T& scaling_factor){
    //Blackman, R.B. and Tukey, J.W., (1958) The measurement of power spectra, Dover Publications, New York.
    //Oppenheim, Alan V., Ronald W. Schafer, and John R. Buck. Discrete-Time Signal Processing. Upper Saddle River, NJ: Prentice Hall, 1999, pp. 468–471
    const T PI=3.1415926535897932384626433832795029L;
    int n_window;
    Wwind.resize(n_size_in);
    if(in_wind_mode=="SYMMETRIC" || in_wind_mode=="Symmetric" || in_wind_mode=="symmetric")
      n_window=n_size_in-1;
    else if(in_wind_mode=="PERIODIC" ||in_wind_mode=="Periodic" || in_wind_mode=="periodic")
      n_window=n_size_in;
    else{
      std::cout<<"Fatal Error: windowing mode is wrong, use either \"periodic\" or \"symmetric\""<<std::endl;
      return;
    }

    for(int i=0; i<n_size_in; i++)
      Wwind[i]= 0.42 - 0.5*cos(2.*PI*i/n_window) + 0.08*cos(4.*PI*i/n_window);

    if(in_scaling_mode==PEAK) // peak magnitude preserving
      scaling_factor=50./21.;
    else if(in_scaling_mode==VARIANCE) // variance preserving
      scaling_factor=(50.*sqrt(3046.))/1523.;
    else{
      std::cout<<"Warning: wrong scaling, use either \"v\" or \"p\""<<std::endl;
      scaling_factor=1000.;
    }

    return;
  }

  template<typename T>
  void FFT<T>::Triangular(const int n_size_in, const int L,const AVGFFT_Mode_Type in_scaling_mode
                          ,vector& Wwind,T& scaling_factor){
    //Oppenheim, Alan V., Ronald W. Schafer, and John R. Buck. Discrete-Time Signal Processing. Upper Saddle River, NJ: Prentice Hall, 1999, pp. 468–471
    // MATLAB and WIKIPEDIA

    if(L==n_size_in){ //MATLAB version
      if(L%2==1){ // odd
        Wwind.resize(n_size_in);
        for(int i=0; i<((L+1)/2); i++)
          Wwind[i]=2.*(i+1)/(L+1.);
        for(int i=((L+1)/2); i<L; i++)
          Wwind[i]=2.-2.*(i+1)/(L+1.);

      }else if(L%2==0){ // even
        Wwind.resize(n_size_in);
        for(int i=0; i<(L/2); i++)
          Wwind[i]=(2.*i+1.)/L;
        for(int i=(L/2); i<L; i++)
          Wwind[i]=2.-(2.*i+1)/L;
      }

    }else{// wikipedia version
      int n_window=n_size_in-1;
      int N_2_=n_window/2,L_2_=L/2;
      Wwind.resize(n_size_in);
      for(int i=0; i<n_size_in; i++)
        Wwind[i]= 1. - fabs((i-N_2_)/L_2_);
    }

    if(in_scaling_mode==PEAK) // peak magnitude preserving
      scaling_factor=2.;
    else if(in_scaling_mode==VARIANCE) // variance preserving
      scaling_factor=sqrt(3.);
    else{
      std::cout<<"Warning: wrong scaling, use either \"v\" or \"p\""<<std::endl;
      scaling_factor=1000.;
    }

    return;
  }

  template<typename T>
  void FFT<T>::Bartlett(const int n_size_in,const std::string in_wind_mode,const AVGFFT_Mode_Type in_scaling_mode
                        ,vector& Wwind,T& scaling_factor){
    // or Bartlett, from--> W.H. Press, B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling, “Numerical Recipes”
    // Oppenheim, Alan V., Ronald W. Schafer, and John R. Buck. Discrete-Time Signal Processing. Upper Saddle River, NJ: Prentice Hall, 1999, pp. 468–471
    int n_window;
    Wwind.resize(n_size_in);

    if(in_wind_mode=="SYMMETRIC" || in_wind_mode=="Symmetric" || in_wind_mode=="symmetric")
      n_window=n_size_in-1; //n_size should be odd
    else if(in_wind_mode=="PERIODIC" ||in_wind_mode=="Periodic" || in_wind_mode=="periodic")
      n_window=n_size_in;   // n_size should be even
    else{
      std::cout<<"Fatal Error: windowing mode is wrong, use either \"periodic\" or \"symmetric\""<<std::endl;
      return;
    }

    int N_2_=n_window/2;

    for(int i=0; i<n_size_in; i++)
      Wwind[i]= 1. - fabs((i-N_2_)/N_2_);

    if(in_scaling_mode==PEAK) // peak magnitude preserving
      scaling_factor=2.;
    else if(in_scaling_mode==VARIANCE) // variance preserving
      scaling_factor=sqrt(3.);
    else{
      std::cout<<"Warning: wrong scaling, use either \"v\" or \"p\""<<std::endl;
      scaling_factor=1000.;
    }

    return;
  }

  template<typename T>
  void FFT<T>::Welch(const int n_size_in,const std::string in_wind_mode,const AVGFFT_Mode_Type in_scaling_mode
                     ,vector& Wwind,T& scaling_factor){
    // W.H. Press, B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling, “Numerical Recipes”
    int n_window;
    Wwind.resize(n_size_in);

    if(in_wind_mode=="SYMMETRIC" || in_wind_mode=="Symmetric" || in_wind_mode=="symmetric")
      n_window=n_size_in-1; //n_size should be odd
    else if(in_wind_mode=="PERIODIC" ||in_wind_mode=="Periodic" || in_wind_mode=="periodic")
      n_window=n_size_in;   // n_size should be even
    else{
      std::cout<<"Fatal Error: windowing mode is wrong, use either \"periodic\" or \"symmetric\""<<std::endl;
      return;
    }

    int N_2_=n_window/2;

    for(int i=0; i<n_size_in; i++)
      Wwind[i]= 1. - pow(((i-N_2_)/N_2_),2);

    if(in_scaling_mode==PEAK) // peak magnitude preserving
      scaling_factor=3./2.;
    else if(in_scaling_mode==VARIANCE) // variance preserving
      scaling_factor=sqrt(30.)/4.;
    else{
      std::cout<<"Warning: wrong scaling, use either \"v\" or \"p\""<<std::endl;
      scaling_factor=1000.;
    }

    return;
  }

  template<typename T>
  T FFT<T>::GetWindowScalingFactor(FFT_WINDOW_Type window_type_in, PSD_Type psd_type_in){
    T scaling=0.;
    if(window_type_in==HANN){
      if(psd_type_in==DENSITY)
        scaling=sqrt(8./3.);
      else if(psd_type_in==POWER)
        scaling=2.;
    }else if(window_type_in==HAMM){
      if(psd_type_in==DENSITY)
        scaling=(50.*sqrt(3974.))/1987.;
      else if(psd_type_in==POWER)
        scaling=50./27.;
    }else if(window_type_in==BLACKMAN){
      if(psd_type_in==DENSITY)
        scaling=(50.*sqrt(3046.))/1523.;
      else if(psd_type_in==POWER)
        scaling=50./21.;
    }else if(window_type_in==TRIANGULAR){
      if(psd_type_in==DENSITY)
        scaling=sqrt(3.);
      else if(psd_type_in==POWER)
        scaling=2.;
    }else if(window_type_in==BARTLETT){
      if(psd_type_in==DENSITY)
        scaling=sqrt(3.);
      else if(psd_type_in==POWER)
        scaling=2.;
    }else if(window_type_in==WELCH){
      if(psd_type_in==DENSITY)
        scaling=sqrt(30.)/4.;
      else if(psd_type_in==POWER)
        scaling=3./2.;
    }else if(window_type_in==RECTANGULAR){
      scaling=1.;
    }else{
      std::cout<<"Warning: Window type is not implemented, assume it is rectangular\n";
      scaling=1.;
    }

    return scaling;
  }

}

#endif
