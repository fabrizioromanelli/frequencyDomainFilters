#include <bitset>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdint>
#include <math.h>

#define PI 3.14159265

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
  uint16_t        NumOfChan;      // Number of channels 1=Mono 2=Sterio
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

int main(int argc, char* argv[])
{
  wav_hdr wavHeader;
  notch_params notchParameters;
  int headerSize = sizeof(wav_hdr), filelength = 0;

  const char* filePath;
  string input;
  if (argc <= 1)
  {
    cout << "Input wave file name: ";
    cin >> input;
    cin.get();
    filePath = input.c_str();
  }
  else
  {
    filePath = argv[1];
    cout << "Input wave file name: " << filePath << endl;
  }

  FILE* wavFile = fopen(filePath, "r");
  if (wavFile == nullptr)
  {
    fprintf(stderr, "Unable to open wave file: %s\n", filePath);
    return 1;
  }

  //Read the header
  size_t bytesRead = fread(&wavHeader, sizeof(char), headerSize, wavFile);
  cout << "Header Read " << bytesRead << " bytes." << endl;

  if (bytesRead > 0)
  {
    //Read the data
    uint16_t bytesPerSample = wavHeader.bitsPerSample / 8;      //Number     of bytes per sample
    uint64_t numSamples = wavHeader.ChunkSize / bytesPerSample; //How many samples are in the wav file?
    short int *value_i = new short int[numSamples];
    double *value_d = new double[numSamples];

    //Reading data
    for (int i = 0; i < numSamples; i++)
      {
        fread(&value_i[i], bytesPerSample, 1, wavFile);
        value_d[i] = value_i[i] / 32768.0;
        //cout << value_i[i] << "|" << value_d[i] << endl;
      }

    filelength = getFileSize(wavFile);

    cout << "File is                    :" << filelength << " bytes." << endl;
    cout << "RIFF header                :" << wavHeader.RIFF[0] << wavHeader.RIFF[1] << wavHeader.RIFF[2] << wavHeader.RIFF[3] << endl;
    cout << "WAVE header                :" << wavHeader.WAVE[0] << wavHeader.WAVE[1] << wavHeader.WAVE[2] << wavHeader.WAVE[3] << endl;
    cout << "FMT                        :" << wavHeader.fmt[0] << wavHeader.fmt[1] << wavHeader.fmt[2] << wavHeader.fmt[3] << endl;
    cout << "Data size                  :" << wavHeader.ChunkSize << endl;

    // Display the sampling Rate from the header
    cout << "Sampling Rate              :" << wavHeader.SamplesPerSec << endl;
    cout << "Number of bits used        :" << wavHeader.bitsPerSample << endl;
    cout << "Number of channels         :" << wavHeader.NumOfChan << endl;
    cout << "Number of bytes per second :" << wavHeader.bytesPerSec << endl;
    cout << "Data length                :" << wavHeader.Subchunk2Size << endl;
    cout << "Audio Format               :" << wavHeader.AudioFormat << endl;
    // Audio format 1=PCM,6=mulaw,7=alaw, 257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM

    cout << "Block align                :" << wavHeader.blockAlign << endl;
    cout << "Data string                :" << wavHeader.Subchunk2ID[0] << wavHeader.Subchunk2ID[1] << wavHeader.Subchunk2ID[2] << wavHeader.Subchunk2ID[3] << endl;

    // Ask the user for notch filter parameters
    cout << "Input bandwidth (0-0.5): ";
    cin >> input;
    cin.get();
    notchParameters.BW = atof(input.c_str());
    cout << "Input sampling frequency: ";
    cin >> input;
    cin.get();
    notchParameters.sampleFreq = atof(input.c_str());
    cout << "Input cut-off frequency ( < " << notchParameters.sampleFreq / 2.0 << "): ";
    cin >> input;
    cin.get();
    notchParameters.cf = atof(input.c_str());

    short int *filtered_i = new short int[numSamples];
    notchFilter(value_i, filtered_i, numSamples, notchParameters);

    const char* filePath2 = "filtered.wav";
    FILE* wavFile2 = fopen(filePath2, "wb");

    fwrite (&wavHeader , sizeof(char), sizeof(wavHeader), wavFile2);
    for (int i = 0; i < numSamples; i++)
      fwrite(&filtered_i[i], bytesPerSample, 1, wavFile2);

    fclose(wavFile2);
  }
  fclose(wavFile);

  return 0;
}

void notchFilter(const short int *in, short int *filtered, int n_samples, notch_params notchParams)
{
  double in_2 = 0.0;
  double in_1 = 0.0;
  double filtered_2 = 0.0;
  double filtered_1 = 0.0;

  // Parameters for the notch filter
  const double BW = notchParams.BW;
  const double cf = notchParams.cf;
  const double  f = cf / notchParams.sampleFreq;
  double K = 0.0;
  double R = 0.0;
  double a0 = 0.0;
  double a1 = 0.0;
  double a2 = 0.0;
  double b1 = 0.0;
  double b2 = 0.0;

  // Compute parameters from BW and f
  R = 1 - 3*BW;
  K = (1 - 2*R*cos(2*PI*f)+R*R) / (2 - 2*cos(2*PI*f));
  cout << "**********************" << endl << "Notch filter paramers:" << endl;
  cout << "R: " << R << " K: " << K << endl;

  // Fill in the FIR coefficients
  a0 = K;
  a1 = -2*K*cos(2*PI*f);
  a2 = K;
  b1 = 2*R*cos(2*PI*f);
  b2 = -R*R;
  cout << "a0: " << a0 << " a1: " << a1 << " a2: " << a2 << endl;
  cout << "b1: " << b1 << " b2: " << b2 << endl;
  cout << "**********************" << endl;

  for (int i = 0; i < n_samples; ++i)
  {
    filtered[i] = a0 * in[i] + a1 * in_1 + a2 * in_2 + b1 * filtered_1 + b2 * filtered_2;
    in_2 = in_1;
    in_1 = in[i];
    filtered_2 = filtered_1;
    filtered_1 = filtered[i];
  }
}

// find the file size
int getFileSize(FILE* inFile)
{
  int fileSize = 0;
  fseek(inFile, 0, SEEK_END);

  fileSize = ftell(inFile);

  fseek(inFile, 0, SEEK_SET);
  return fileSize;
}
