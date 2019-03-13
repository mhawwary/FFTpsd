// C++ includes
#include <algorithm>
#include <map>
#include <sstream>

// Local includes
#include "string_to_type.h"
#include "../include/error.h"
#include "fft_type.h"


// ------------------------------------------------------------
// Anonymous namespace to hold local data & methods
namespace {
// Reverse a map
template <typename MapIter, class MapType>
inline
void build_reverse_map (MapIter it, MapIter end, MapType& reverse)
{
  reverse.clear();

  for (; it != end; ++it)
  {
    // If the forward map is not invertible, we might already have
    // found a preimage of it->second.  Choose the "largest"
    // preimage according to operator<; for std::string this will
    // give us the longest, hopefully most specific name
    // corresponding to an enum.
    typename MapType::iterator preimage = reverse.find(it->second);
    if (preimage == reverse.end())
      reverse.insert (std::make_pair(it->second, it->first));
    else if (preimage->second < it->first)
      preimage->second = it->first;
  }
}


inline void string_to_enum_error( const std::string& _s )
{
  std::stringstream ss;
  ss << "\nCannot convert string\"" << _s << "\" to enumeration." << std::endl;
  FatalError_exit( ss.str().c_str() );
}

// Undefined at end of file
#define enum_to_string_error( _eType, _eValue )                                                                                     \
  do{                                                                                                                             \
  std::stringstream ss;                                                                                                           \
  ss<<"\nCannot convert the enumeration value " << #_eType << " = " << static_cast<int>( _eValue ) << " to a string." << std::endl;\
  FatalError_exit(ss.str().c_str());                                                                                                          \
}while(0)

} // end anonymous namespace


/**
 * FFT Window type
 */
std::map<std::string, fftspc::FFT_WINDOW_Type> fftwindow_type_to_enum;
void init_fftwindow_type_to_enum ()
{
  if (fftwindow_type_to_enum.empty())
  {
    fftwindow_type_to_enum["HANN"]=fftspc::HANN;
    fftwindow_type_to_enum["HAMM"]=fftspc::HAMM;
    fftwindow_type_to_enum["BLACKMAN"]=fftspc::BLACKMAN;
    fftwindow_type_to_enum["WELCH"]=fftspc::WELCH;
    fftwindow_type_to_enum["BARTLETT"]=fftspc::BARTLETT;
    fftwindow_type_to_enum["TRIANGULAR"]=fftspc::TRIANGULAR;
    fftwindow_type_to_enum["RECTANGULAR"]=fftspc::RECTANGULAR;
  }
}

std::map<fftspc::FFT_WINDOW_Type, std::string> enum_to_fftwindow_type;

void init_enum_to_fftwindow_type ()
{
  // Build reverse map
  if (enum_to_fftwindow_type.empty())
  {
    init_fftwindow_type_to_enum();

    build_reverse_map (fftwindow_type_to_enum.begin(),
                       fftwindow_type_to_enum.end(),
                       enum_to_fftwindow_type);
  }
}

// specialization
template <>
fftspc::FFT_WINDOW_Type string_to_enum<fftspc::FFT_WINDOW_Type> (const std::string& s)
{
  init_fftwindow_type_to_enum();

  std::string upper(s);
  std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

  if (!fftwindow_type_to_enum.count(upper))
    string_to_enum_error(s);

  return fftwindow_type_to_enum[upper];
}

template <>
std::string enum_to_string<fftspc::FFT_WINDOW_Type> (const fftspc::FFT_WINDOW_Type i)
{
  init_enum_to_fftwindow_type();

  if (!enum_to_fftwindow_type.count(i))
    enum_to_string_error( fftspc::FFT_WINDOW_Type, i );

  return enum_to_fftwindow_type[i];
}

//--------------------------------------------------------
/**
 * Wave Form Type
 */
std::map<std::string, fftspc::Wave_Form_Type> waveform_type_to_enum;
void init_waveform_type_to_enum ()
{
  if (waveform_type_to_enum.empty())
  {
    waveform_type_to_enum["SINE"]=fftspc::SINE;
    waveform_type_to_enum["COSINE"]=fftspc::COSINE;
    waveform_type_to_enum["GAUSSIAN"]=fftspc::GAUSSIAN;
  }
}

std::map<fftspc::Wave_Form_Type, std::string> enum_to_waveform_type;

void init_enum_to_waveform_type ()
{
  // Build reverse map
  if (enum_to_waveform_type.empty())
  {
    init_waveform_type_to_enum();

    build_reverse_map (waveform_type_to_enum.begin(),
                       waveform_type_to_enum.end(),
                       enum_to_waveform_type);
  }
}

// specialization
template <>
fftspc::Wave_Form_Type string_to_enum<fftspc::Wave_Form_Type> (const std::string& s)
{
  init_waveform_type_to_enum();

  std::string upper(s);
  std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

  if (!waveform_type_to_enum.count(upper))
    string_to_enum_error(s);

  return waveform_type_to_enum[upper];
}

template <>
std::string enum_to_string<fftspc::Wave_Form_Type> (const fftspc::Wave_Form_Type i)
{
  init_enum_to_waveform_type();

  if (!enum_to_waveform_type.count(i))
    enum_to_string_error( fftspc::Wave_Form_Type, i );

  return enum_to_waveform_type[i];
}

//--------------------------------------------------------
/**
 * PSD Type
 */
std::map<std::string, fftspc::PSD_Type> psd_type_to_enum;
void init_psd_type_to_enum ()
{
  if (psd_type_to_enum.empty())
  {
    psd_type_to_enum["DENSITY"]=fftspc::DENSITY;
    psd_type_to_enum["POWER"]=fftspc::POWER;
  }
}

std::map<fftspc::PSD_Type, std::string> enum_to_psd_type;

void init_enum_to_psd_type ()
{
  // Build reverse map
  if (enum_to_psd_type.empty())
  {
    init_psd_type_to_enum();

    build_reverse_map (psd_type_to_enum.begin(),
                       psd_type_to_enum.end(),
                       enum_to_psd_type);
  }
}

// specialization
template <>
fftspc::PSD_Type string_to_enum<fftspc::PSD_Type> (const std::string& s)
{
  init_psd_type_to_enum();

  std::string upper(s);
  std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

  if (!psd_type_to_enum.count(upper))
    string_to_enum_error(s);

  return psd_type_to_enum[upper];
}

template <>
std::string enum_to_string<fftspc::PSD_Type> (const fftspc::PSD_Type i)
{
  init_enum_to_psd_type();

  if (!enum_to_psd_type.count(i))
    enum_to_string_error( fftspc::PSD_Type, i );

  return enum_to_psd_type[i];
}

//--------------------------------------------------------
/**
 * fftspc::AVGFFT_Mode_Type
 */
std::map<std::string, fftspc::AVGFFT_Mode_Type> avgfft_mode_type_to_enum;
void init_avgfft_mode_type_to_enum ()
{
  if (avgfft_mode_type_to_enum.empty())
  {
    avgfft_mode_type_to_enum["PEAK"]=fftspc::PEAK;
    avgfft_mode_type_to_enum["VARIANCE"]=fftspc::VARIANCE;
    avgfft_mode_type_to_enum["ENSEMBLE"]=fftspc::ENSEMBLE;
  }
}

std::map<fftspc::AVGFFT_Mode_Type, std::string> enum_to_avgfft_mode_type;

void init_enum_to_avgfft_mode_type ()
{
  // Build reverse map
  if (enum_to_avgfft_mode_type.empty())
  {
    init_avgfft_mode_type_to_enum();

    build_reverse_map (avgfft_mode_type_to_enum.begin(),
                       avgfft_mode_type_to_enum.end(),
                       enum_to_avgfft_mode_type);
  }
}

// specialization
template <>
fftspc::AVGFFT_Mode_Type string_to_enum<fftspc::AVGFFT_Mode_Type> (const std::string& s)
{
  init_avgfft_mode_type_to_enum();

  std::string upper(s);
  std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

  if (!avgfft_mode_type_to_enum.count(upper))
    string_to_enum_error(s);

  return avgfft_mode_type_to_enum[upper];
}

template <>
std::string enum_to_string<fftspc::AVGFFT_Mode_Type> (const fftspc::AVGFFT_Mode_Type i)
{
  init_enum_to_avgfft_mode_type();

  if (!enum_to_avgfft_mode_type.count(i))
    enum_to_string_error( fftspc::AVGFFT_Mode_Type, i );

  return enum_to_avgfft_mode_type[i];
}

#undef enum_to_string_error


