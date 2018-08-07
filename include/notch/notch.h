#ifndef NOTCH_H
#define NOTCH_H

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

typedef struct NOTCH_PARAMETERS
{
  double BW = 0.0;         // bandwidth
  double cf = 0.0;         // cutoff frequency
  double sampleFreq = 0.0; // frequency
  //  const double BW = 0.0066;          // bandwidth
  //  const double cf = 440;             // cutoff frequency
  //  const double sampleFreq = 44100;   // sampling frequency
} notch_params;

// Function prototypes
int getFileSize(FILE* inFile);
void notchFilter(const short int *in, short int *filtered, int n_samples, notch_params notchParams);

#endif // NOTCH_H
