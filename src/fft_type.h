#ifndef FFT_TYPE_H
#define FFT_TYPE_H

namespace fftspc{

  enum FFT_WINDOW_Type{
    HANN,
    HAMM,
    BLACKMAN,
    WELCH,
    BARTLETT,
    TRIANGULAR,
    RECTANGULAR      // is equivalent to No windowing
  };

  enum Wave_Form_Type{
    SINE,
    COSINE,
    GAUSSIAN
  };

  enum PSD_Type{
    DENSITY,        //!< refers to power spectral density
    POWER           //!< refers to power spectrum
  };

  enum AVGFFT_Mode_Type{
    PEAK,
    VARIANCE,
    ENSEMBLE
  };

}

#endif // FFT_TYPE_H

