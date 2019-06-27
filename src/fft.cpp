#include"fft.h"

namespace fftspc{

template<typename T>
FFT<T>::~FFT(){}

//------------------------------------------------
// Private functions:
//------------------------------------------------
template<typename T>
void FFT<T>::Init(const int n_data_size){
  n_data=n_data_size;
  n_data_levels = std::log2(n_data);
  n_2_data=n_data/2;
  W_arr.clear();
  Wi_arr.clear();
  m_arr.clear();
  m_2_arr.clear();
  W_arr.resize(n_data_levels);
  Wi_arr.resize(n_data_levels);
  m_arr.resize(n_data_levels);
  m_2_arr.resize(n_data_levels);

  for(register int s=0; s<n_data_levels; s++){
    complex ws = std::exp(-2.*I*PI/std::pow(2,s+1));
    complex ws_inv = std::exp(2.*I*PI/std::pow(2,s+1));
    m_arr[s]=pow(2,s+1);
    m_2_arr[s]=m_arr[s]/2;
    for(int j=0; j<m_2_arr[s]; j++){
      W_arr[s].push_back(pow(ws,j));
      Wi_arr[s].push_back(pow(ws_inv,j));
    }
  }

  return;
}

template<typename T>
void FFT<T>::Init_real(const int n_data_size){
  n_data=n_data_size;    // n
  n_data_levels = std::log2(n_data);
  n_2_data=n_data/2;   // n/2
  n_2_data_levels = std::log2(n_2_data);
  n_4_data=n_2_data/2;  // n/4
  W_arr.clear();
  Wi_arr.clear();
  m_arr.clear();
  m_2_arr.clear();
  W_arr.resize(n_2_data_levels);
  Wi_arr.resize(n_2_data_levels);
  m_arr.resize(n_2_data_levels);
  m_2_arr.resize(n_2_data_levels);

  for(register int s=0; s<n_2_data_levels; s++){
    complex ws = std::exp(-2.*I*PI/std::pow(2,s+1));
    complex ws_inv = std::exp(2.*I*PI/std::pow(2,s+1));
    m_arr[s]=pow(2,s+1);
    m_2_arr[s]=m_arr[s]/2;
    for(int j=0; j<m_2_arr[s]; j++){
      W_arr[s].push_back(pow(ws,j));
      Wi_arr[s].push_back(pow(ws_inv,j));
    }
  }

  T nn_n;
  complex argument,argument2;
  Wn.clear();
  Wn2.clear();
  Wn.resize(n_4_data+1);
  Wn2.resize(n_4_data+1);
  for(register int s=0; s<n_4_data+1; s++){
    nn_n=(T) s/n_data;
    argument  = (-2*nn_n+0.5)*PI*I;
    argument2 = ( 2*nn_n-0.5)*PI*I;
    Wn[s]  = std::exp(argument);
    Wn2[s] = std::exp(argument2);
  }

  return;
}

template<typename T>
void FFT<T>::Init_dft(const int n_data_size){
  n_data=n_data_size;
  int fftsize=0;
  if (n_data%2==0)
    fftsize = int(n_data/2)+1;
  else
    fftsize = int((n_data-1)/2)+1;
  W_arr.clear();
  Wi_arr.clear();
  W_arr.resize(fftsize);
  Wi_arr.resize(fftsize);
  complex factor0=I*(2.*PI/n_data),factor;
  for(int s=0; s<fftsize; s++){
    W_arr[s].resize(n_data,I0);
    Wi_arr[s].resize(n_data,I0);
    for(int j=0; j<n_data; j++){
      factor = (T) s*j*factor0;
      W_arr[s][j]=exp(-factor);
      Wi_arr[s][j]=exp(factor);
    }
  }

  return;
}

// An inefficient bit reverse
//-Reference: Numerical Recipes Press et. al 2007 3rd edition
template<typename T>
template<typename TT>
void FFT<T>::bit_reverse(TT &u_data, const int in_data_size){
  int local_n=in_data_size;
  int local_n_2 = local_n/2;  // so it works only if n is a power of 2
  int j=1,m;
  for(int i=1; i<local_n; i++){
    if(j>i)
      std::swap(u_data[i-1],u_data[j-1]);
    m = local_n_2;
    while( m>=2 && j>m ){
      j-=m;
      m=m/2;
    }
    j+=m;
  }

  return;
}

template<typename T>
void FFT<T>::fft_inplace(cvector &u_fft, const int local_n, const int local_n_levels){
  bit_reverse(u_fft,local_n);
  for(register int s=0; s<local_n_levels; s++){
    for(int k=0; k<local_n; k+=m_arr[s]){
      complex y0,temp;
      for(int j=0; j<m_2_arr[s]; j++){
        temp = W_arr[s][j]*u_fft[k+j+m_2_arr[s]];
        y0 = u_fft[k+j];
        u_fft[k+j]= y0 + temp;
        u_fft[k+j+m_2_arr[s]]=y0-temp;
      }
    }
  }
  return;
}

template<typename T>
void FFT<T>::ifft_inplace(cvector &u_data, const int local_n, const int local_n_levels){
  bit_reverse(u_data,local_n);
  for(register int s=0; s<local_n_levels; s++){
    for(int k=0; k<local_n; k+=m_arr[s]){
      complex y0,temp;
      for(int j=0; j<m_2_arr[s]; j++){
        temp = Wi_arr[s][j]*u_data[k+j+m_2_arr[s]];
        y0 = u_data[k+j];
        u_data[k+j]= y0 + temp;
        u_data[k+j+m_2_arr[s]]=y0-temp;
      }
    }
  }

  for(register int j=0;j<local_n; j++)
    u_data[j]/=local_n;
  return;
}

template<typename T>
void FFT<T>::packing_oddeven_real2complexdata(const vector u_data, cvector &h_packed, const int n_packed){
  int k=0;
  for(register int j=0; j<n_packed; j++){
    //h_packed.push_back(u_data[k] + u_data[k+1]*I);
    h_packed[j]=u_data[k] + u_data[k+1]*I;
    k+=2;
  }
  return;
}

template<typename T>
void FFT<T>::folding_realimag_complex2realdata(const cvector h_packed, vector &u_data, const int n_packed){
  int k=0;
  for(register int j=0; j<n_packed; j++){
    u_data[k]=h_packed[j].real();
    u_data[k+1]=h_packed[j].imag();
    k+=2;
  }
  return;
}

//------------------------------------------------
// Public functions:
//------------------------------------------------
template<typename T>
void FFT<T>::dft(const vector u_data, cvector &u_fft){
  int fftsize=0;
  if (u_data.size()%2==0)
    fftsize = int(u_data.size()/2)+1;
  else
    fftsize = int((u_data.size()-1)/2)+1;
  if(u_fft.size()!=fftsize)
    u_fft.resize(fftsize); // for safety
  std::fill_n(u_fft.begin(),fftsize,I0);

  if(u_data.size()!=n_data || dft_mode==false){
    n_data = u_data.size();
    dft_mode=true;
    Init_dft(n_data);
    std::cout<<"\nReset the infrastructure of complex fft..........\n\n";
  }

  for(int s=0; s<fftsize; s++)
    for(int j=0; j<n_data; j++)
      u_fft[s]+=u_data[j]*W_arr[s][j];

  return;
}

template<typename T>
template<typename T1, typename T2>
void FFT<T>::fft(const T1 u_data, T2 &u_fft){
  u_fft.resize(u_data.size()); // for safety
  std::copy(u_data.begin(),u_data.end(),u_fft.begin());
  if(u_data.size()!=n_data || real_mode==true){
    n_data = u_data.size();
    real_mode=false;
    Init(n_data);
    std::cout<<"\nReset the infrastructure of complex fft..........\n\n";
  }
  fft_inplace(u_fft,n_data,n_data_levels);
  return;
}

template<typename T>
void FFT<T>::ifft(const cvector u_fft, cvector &u_data){
  u_data.resize(u_fft.size()); // for safety
  std::copy(u_fft.begin(),u_fft.end(),u_data.begin());
  if(u_fft.size()!=n_data || real_mode==true){
    n_data = u_data.size();
    real_mode=false;
    Init(n_data);
    std::cout<<"\nReset the infrastructure of complex fft..........\n\n";
  }
  ifft_inplace(u_data,n_data,n_data_levels);

  return;
}

template<typename T>
void FFT<T>::fftfreq(const int local_n, const T sample_spacing, vector &freq_){
  freq_.clear();
  freq_.resize(local_n);
  int local_n_2=local_n/2;
  T factor = sample_spacing*local_n;
  freq_[0]=0;  // zero/DC or mean value frequency
  for(register int i=1; i<local_n_2; i++){
    freq_[i]=i/factor;
    freq_[local_n-i]=-freq_[i];
  }
  freq_[local_n_2]=-local_n_2/factor;  // Nyquist cuttoff
  return;
}

template<typename T>
void FFT<T>::fftintfreq(const int local_n, ivector &freq_){
  freq_.clear();
  freq_.resize(local_n);
  int local_n_2=local_n/2;
  freq_[0]=0;  // zero/DC or mean value frequency
  for(register int i=1; i<local_n_2; i++){
    freq_[i]=i;
    freq_[local_n-i]=-i;
  }
  freq_[local_n_2]=-local_n_2;  // Nyquist cuttoff

  return;
}

template<typename T>
template<typename TT>
void FFT<T>::fftshift(TT &u_fft){
  int n_2_local = u_fft.size()/2;
  // Swapping the 0 and Nyquist frequencies
  std::swap(u_fft[n_2_local],u_fft[0]);
  // Swapping the negative and positive frequencies
  std::swap_ranges(u_fft.begin()+1,u_fft.begin()+n_2_local,u_fft.begin()+n_2_local+1);
  return;
}

// Real fft functions
template<typename T>
void FFT<T>::rfft(const vector u_data, cvector &u_fft){
  // First: split the data into even and odd and pack them together to form
  // a complex array of size n/2
  if(u_data.size()!=n_data || real_mode==false){
    n_data=u_data.size();
    real_mode=true;
    Init_real(n_data);
    std::cout<<"\nReset the infrastructure of real fft..........\n\n";
  }

  /*  cvector h_fft(n_2_data);
  //  packing_oddeven_real2complexdata(u_data,h_fft,n_2_data);
  //  fft_inplace(h_fft,n_2_data,n_2_data_levels);  // compute the FFT of packed data in place
  //  h_fft.push_back(h_fft[0]);  // length is now N/2+1, adding u_fft_[N/2]=u_fft_[0]
  //  u_fft.resize(n_2_data+1); // the output FFT of size N/2+1
  //  register int i;
  //  #pragma omp parallel for //num_threads(1)
  //  for(i=0; i<n_4_data+1; i++){ // does not necessarily need to be done in place, so some parallization is possible
  //    complex h1,h2,conj_h1,conj_h2,conj_h_fft,conjn2_h_fft;
  //    conj_h_fft=std::conj(h_fft[i]);
  //    conjn2_h_fft=std::conj(h_fft[n_2_data-i]);
  //    h1 = 0.5*(h_fft[i]+conjn2_h_fft);
  //    h2 = 0.5*(h_fft[i]-conjn2_h_fft);
  //    conj_h1 = 0.5*(conj_h_fft+h_fft[n_2_data-i]) ;
  //    conj_h2 = 0.5*(conj_h_fft-h_fft[n_2_data-i]) ;
  //    u_fft[i]=          h1 - Wn[i] *h2       ;  // u_{i}
  //    u_fft[n_2_data-i]= conj_h1 + Wn2[i]*conj_h2  ;  // u_{N/2-i}
  //  }*/

  u_fft.resize(n_2_data); // the output FFT of size N/2 and for saftey
  packing_oddeven_real2complexdata(u_data,u_fft,n_2_data);
  fft_inplace(u_fft,n_2_data,n_2_data_levels);  // compute the FFT of packed data in place
  u_fft.push_back(u_fft[0]);  // length is now N/2+1, adding u_fft_[N/2]=u_fft_[0]
  // Folding loop to compute the final output
  complex h1,h2,conj_h1,conj_h2,conj_h_fft,conjn2_h_fft;
  for(register int i=0; i<n_4_data+1; i++){ // does not necessarily need to be done in place, so some parallization is possible
    conj_h_fft=std::conj(u_fft[i]);
    conjn2_h_fft=std::conj(u_fft[n_2_data-i]);
    h1 = 0.5*(u_fft[i]+conjn2_h_fft);
    h2 = 0.5*(u_fft[i]-conjn2_h_fft);
    conj_h1 = 0.5*(conj_h_fft+u_fft[n_2_data-i]) ;
    conj_h2 = 0.5*(conj_h_fft-u_fft[n_2_data-i]) ;
    u_fft[i]=          h1 - Wn[i] *h2       ;  // u_{i}
    u_fft[n_2_data-i]= conj_h1 + Wn2[i]*conj_h2  ;  // u_{N/2-i}
  }
  return;
}

template<typename T>
void FFT<T>::irfft(const cvector h_fft, vector &u_data){
  // First: split the data into even and odd and pack them together to form
  // a complex array of size n/2
  if(2*(h_fft.size()-1)!=n_data || real_mode==false){
    n_data = 2*(h_fft.size()-1);
    real_mode=true;
    Init_real(n_data);
    std::cout<<"\nReset the infrastructure of real fft..........\n\n";
  }
  cvector h_fft_(n_2_data+1);
  std::copy(h_fft.begin(),h_fft.end(),h_fft_.begin());
  // Folding loop to compute the final output
  complex h1,h2,conj_h1,conj_h2,conj_h_fft,conjn2_h_fft;
  for(register int i=0; i<n_4_data+1; i++){
    conj_h_fft=std::conj(h_fft_[i]);
    conjn2_h_fft=std::conj(h_fft_[n_2_data-i]);
    h1 = 0.5*(h_fft_[i]+conjn2_h_fft);
    h2 = 0.5*(h_fft_[i]-conjn2_h_fft);
    conj_h1 = 0.5*(conj_h_fft+h_fft_[n_2_data-i]) ;
    conj_h2 = 0.5*(conj_h_fft-h_fft_[n_2_data-i]) ;
    h_fft_[i]=          h1 - Wn2[i] *h2       ;  // u_{i}
    h_fft_[n_2_data-i]= conj_h1 + Wn[i]*conj_h2  ;  // u_{N/2-i}
  }
  u_data.resize(n_data); // for safety
  h_fft_.pop_back(); // need to erase the last element of h_fft_, h_fft_[0]=h_fft_[end]
  ifft_inplace(h_fft_,n_2_data,n_2_data_levels);
  folding_realimag_complex2realdata(h_fft_,u_data,n_2_data);

  return;
}

template<typename T>
void FFT<T>::rfftfreq(const int local_n, const T sample_spacing, vector &freq_){
  int local_n_2=local_n/2;
  freq_.clear();
  freq_.resize(local_n_2+1);
  T factor = sample_spacing*local_n;
  for(register int i=0; i<local_n_2+1; i++)
    freq_[i]=i/factor;
  return;
}

template<typename T>
void FFT<T>::rfftintfreq(const int local_n, ivector &freq_){
  int local_n_2=local_n/2;
  freq_.clear();
  freq_.resize(local_n_2+1);
  for(register int i=0; i<local_n_2+1; i++)
    freq_[i]=i;

  return;
}

template class FFT<double>;
template void FFT<double>::bit_reverse<>(std::vector<double>&, const int);
template void FFT<double>::bit_reverse<>(std::vector<std::complex<double>>&, const int);
template void FFT<double>::fft<>(const std::vector<double>, std::vector<std::complex<double>>&);
template void FFT<double>::fft<>(const std::vector<std::complex<double>>, std::vector<std::complex<double>>&);
template void FFT<double>::fftshift<>(std::vector<std::complex<double>>&);
template void FFT<double>::fftshift<>(std::vector<double>&);
template void FFT<double>::fftshift<>(std::vector<int>&);

}
//template class FFT<long double>;
//template void FFT<long double>::bit_reverse<>(std::vector<long double>&, const int);
//template void FFT<long double>::bit_reverse<>(std::vector<std::complex<long double>>&, const int);
//template void FFT<long double>::fft<>(const std::vector<long double>, std::vector<std::complex<long double>>&);
//template void FFT<long double>::fft<>(const std::vector<std::complex<long double>>, std::vector<std::complex<long double>>&);
//template void FFT<long double>::fftshift<>(std::vector<std::complex<long double>>&);
//template void FFT<long double>::fftshift<>(std::vector<long double>&);
//template void FFT<long double>::fftshift<>(std::vector<int>&);
