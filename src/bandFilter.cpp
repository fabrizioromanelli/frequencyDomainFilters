#include <bandFilter/bandFilter.h>

/**
 *  @brief This function computes the band-pass or band-reject (notch) filter offline given the input parameters.
 *  @param input in: input signal
 *  @param input n_samples: number of signal samples
 *  @param input params: parameters for the band filter
 *  @param input filtered: out filtered signal
 */
//template <class T, class U>
//void bandOffline(const T *in, T *filtered, int n_samples, U params)
//{
//  double in_2 = 0.0;
//  double in_1 = 0.0;
//  double filtered_2 = 0.0;
//  double filtered_1 = 0.0;

//  for (int i = 0; i < n_samples; ++i)
//  {
//    filtered[i] = params.a0 * in[i] + params.a1 * in_1 + params.a2 * in_2 + params.b1 * filtered_1 + params.b2 * filtered_2;
//    in_2 = in_1;
//    in_1 = in[i];
//    filtered_2 = filtered_1;
//    filtered_1 = filtered[i];
//  }
//}

/**
 *  @brief This function computes the band-pass or band-reject (notch) filter online given the input parameters.
 *  @param input in: input signal
 *  @param input params: parameters for the band filter
 *  @param input filtered: out filtered signal
 */
//template <class T, class U>
//void bandOnline(const T in, T &filtered, U params)
//{
//  static double in_2 = 0.0;
//  static double in_1 = 0.0;
//  static double filtered_2 = 0.0;
//  static double filtered_1 = 0.0;

//  filtered = params.a0 * in + params.a1 * in_1 + params.a2 * in_2 + params.b1 * filtered_1 + params.b2 * filtered_2;
//  in_2 = in_1;
//  in_1 = in;
//  filtered_2 = filtered_1;
//  filtered_1 = filtered;
//}
