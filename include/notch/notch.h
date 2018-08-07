#ifndef NOTCH_H
#define NOTCH_H

#include <bitset>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdint>
#include <math.h>

using std::cin;
using std::cout;
using std::endl;
using std::fstream;
using std::string;

typedef struct WAV_HEADER
{
  /* RIFF Chunk Descriptor */
  uint8_t         RIFF[4];        // RIFF Header Magic header
  uint32_t        ChunkSize;      // RIFF Chunk Size
  uint8_t         WAVE[4];        // WAVE Header
  /* "fmt" sub-chunk */
  uint8_t         fmt[4];         // FMT header
  uint32_t        Subchunk1Size;  // Size of the fmt chunk
  uint16_t        AudioFormat;    // Audio format 1=PCM,6=mulaw,7=alaw,     257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM
  uint16_t        NumOfChan;      // Number of channels 1=Mono 2=Stereo
  uint32_t        SamplesPerSec;  // Sampling Frequency in Hz
  uint32_t        bytesPerSec;    // bytes per second
  uint16_t        blockAlign;     // 2=16-bit mono, 4=16-bit stereo
  uint16_t        bitsPerSample;  // Number of bits per sample
  /* "data" sub-chunk */
  uint8_t         Subchunk2ID[4]; // "data"  string
  uint32_t        Subchunk2Size;  // Sampled data length
} wav_hdr;

struct notch_params
{
  double a0, a1, a2, b1, b2;

  //  const double BW = 0.0066;          // bandwidth
  //  const double cf = 440;             // cutoff frequency
  //  const double sampleFreq = 44100;   // sampling frequency
  notch_params(double BW, double cf, double sampleFreq)
  {
    // Compute parameters from BW and f
    double f = cf / sampleFreq;
    double R = 1 - 3*BW;
    double K = (1 - 2*R*cos(2*M_PI*f)+R*R) / (2 - 2*cos(2*M_PI*f));
    cout << "**********************" << endl << "Notch filter paramers:" << endl;
    cout << "R: " << R << " K: " << K << endl;

    // Fill in the FIR coefficients
    a0 = K;
    a1 = -2*K*cos(2*M_PI*f);
    a2 = K;
    b1 = 2*R*cos(2*M_PI*f);
    b2 = -R*R;
    cout << "a0: " << a0 << " a1: " << a1 << " a2: " << a2 << endl;
    cout << "b1: " << b1 << " b2: " << b2 << endl;
    cout << "**********************" << endl;
  }
};

// Function prototypes
int getFileSize(FILE* inFile);
template <class T> void notchFilterIntegerOffline(const T *in, T *filtered, int n_samples, notch_params notchParams);
template <class T> void notchFilterIntegerOnline(const T in, T &filtered, notch_params notchParams);

#endif // NOTCH_H
