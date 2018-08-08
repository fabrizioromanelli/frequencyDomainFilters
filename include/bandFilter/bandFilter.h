#ifndef BAND_FILTER_H
#define BAND_FILTER_H

#include <iostream>
#include <math.h>

//#define DBGMSG 1

namespace filter
{
  /**
    * @brief The bandPassParams struct
    * Data structure for the band-pass filter.
    */
  struct bandPassParams
  {
    double a0, a1, a2, b1, b2;

    bandPassParams(double _BW, double cf, double sampleFreq)
    {
      // Compute parameters from BW and f
      double f  = cf / sampleFreq;  // Adjust it from Hz to percentage of sample frequency
      double BW = _BW / sampleFreq; // Adjust it from Hz to percentage of sample frequency
      double R  = 1 - 3*BW;
      double K  = (1 - 2*R*cos(2*M_PI*f)+R*R) / (2 - 2*cos(2*M_PI*f));
#ifdef DBGMSG
      std::cout << "**********************" << std::endl << "Band-pass filter paramers:" << std::endl;
      std::cout << "R: " << R << " K: " << K << std::endl;
#endif

      // Fill in the FIR coefficients
      a0 = 1 - K;
      a1 = 2*(K-R)*cos(2*M_PI*f);
      a2 = R*R - K;
      b1 = 2*R*cos(2*M_PI*f);
      b2 = -R*R;
#ifdef DBGMSG
      std::cout << "a0: " << a0 << " a1: " << a1 << " a2: " << a2 << std::endl;
      std::cout << "b1: " << b1 << " b2: " << b2 << std::endl;
      std::cout << "**********************" << std::endl;
#endif
    }
  };

  /**
    * @brief The bandRejectParams struct
    * Data structure for the band-reject (notch) filter.
    */
  struct bandRejectParams
  {
    double a0, a1, a2, b1, b2;

    bandRejectParams(double _BW, double cf, double sampleFreq)
    {
      // Compute parameters from BW and f
      double f  = cf / sampleFreq;  // Adjust it from Hz to percentage of sample frequency
      double BW = _BW / sampleFreq; // Adjust it from Hz to percentage of sample frequency
      double R  = 1 - 3*BW;
      double K  = (1 - 2*R*cos(2*M_PI*f)+R*R) / (2 - 2*cos(2*M_PI*f));
#ifdef DBGMSG
      std::cout << "**********************" << std::endl << "Notch filter paramers:" << std::endl;
      std::cout << "R: " << R << " K: " << K << std::endl;
#endif

      // Fill in the FIR coefficients
      a0 = K;
      a1 = -2*K*cos(2*M_PI*f);
      a2 = K;
      b1 = 2*R*cos(2*M_PI*f);
      b2 = -R*R;
#ifdef DBGMSG
      std::cout << "a0: " << a0 << " a1: " << a1 << " a2: " << a2 << std::endl;
      std::cout << "b1: " << b1 << " b2: " << b2 << std::endl;
      std::cout << "**********************" << std::endl;
#endif
    }
  };

  template <typename T, typename U>
  class band
  {
    public:
      band(U params) : _a0(params.a0), _a1(params.a1), _a2(params.a2), _b1(params.b1), _b2(params.b2),
                       _in_1(0), _in_2(0), _filtered_1(0), _filtered_2(0) {}

      /**
       *  @brief This function visualizes the coefficients from the FIR.
       *  @param none
       */
      void visualizeCoefficients(void) { std::cout << "A0: " << _a0 << " A1: " << _a1 << " A2: " << _a2 << " B1: " << _b1 << " B2: " << _b2 << std::endl; }

      /**
       *  @brief This function visualizes the coefficients from the FIR.
       *  @param output a0, a1, a2, b1, b2: coefficients from the FIR
       */
      void getCoefficients(double &a0, double &a1, double &a2, double &b1, double &b2) { a0 = _a0; a1 = _a1; a2 = _a2; b1 = _b1; b2 = _b2; }

      /**
       *  @brief This function computes the band-pass or band-reject (notch) filter offline given the input parameters.
       *  @param input in: input signal
       *  @param input n_samples: number of signal samples
       *  @param input params: parameters for the band filter
       *  @param input filtered: out filtered signal
       */
      void offlineUpdate(const T *in, T *filtered, int n_samples)
      {
        for (int i = 0; i < n_samples; ++i)
        {
          filtered[i] = _a0 * in[i] + _a1 * _in_1 + _a2 * _in_2 + _b1 * _filtered_1 + _b2 * _filtered_2;
          _in_2 = _in_1;
          _in_1 = in[i];
          _filtered_2 = _filtered_1;
          _filtered_1 = filtered[i];
        }
      }

      /**
       *  @brief This function computes the band-pass or band-reject (notch) filter online given the input parameters.
       *  @param input in: input signal
       *  @param input params: parameters for the band filter
       *  @param input filtered: out filtered signal
       */
      void onlineUpdate(const T in, T &filtered)
      {
        filtered = _a0 * in + _a1 * _in_1 + _a2 * _in_2 + _b1 * _filtered_1 + _b2 * _filtered_2;
        _in_2 = _in_1;
        _in_1 = in;
        _filtered_2 = _filtered_1;
        _filtered_1 = filtered;
      }

    private:
      double _a0, _a1, _a2, _b1, _b2;
      T      _in_1, _in_2, _filtered_1, _filtered_2;
  };
} // end namespace filter

#endif // BAND_FILTER_H
