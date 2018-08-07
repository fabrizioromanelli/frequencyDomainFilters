#include <notch/notch.h>

int main(int argc, char* argv[])
{
  wav_hdr wavHeader;
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
      }

    filelength = getFileSize(wavFile);

    cout << "File is                    : " << filelength << " bytes." << endl;
    cout << "RIFF header                : " << wavHeader.RIFF[0] << wavHeader.RIFF[1] << wavHeader.RIFF[2] << wavHeader.RIFF[3] << endl;
    cout << "WAVE header                : " << wavHeader.WAVE[0] << wavHeader.WAVE[1] << wavHeader.WAVE[2] << wavHeader.WAVE[3] << endl;
    cout << "FMT                        : " << wavHeader.fmt[0] << wavHeader.fmt[1] << wavHeader.fmt[2] << wavHeader.fmt[3] << endl;
    cout << "Data size                  : " << wavHeader.ChunkSize << endl;

    // Display the sampling Rate from the header
    cout << "Sampling Rate              : " << wavHeader.SamplesPerSec << endl;
    cout << "Number of bits used        : " << wavHeader.bitsPerSample << endl;
    cout << "Number of channels         : " << wavHeader.NumOfChan << endl;
    cout << "Number of bytes per second : " << wavHeader.bytesPerSec << endl;
    cout << "Data length                : " << wavHeader.Subchunk2Size << endl;
    cout << "Audio Format               : " << wavHeader.AudioFormat << endl;
    // Audio format 1=PCM,6=mulaw,7=alaw, 257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM

    cout << "Block align                : " << wavHeader.blockAlign << endl;
    cout << "Data string                : " << wavHeader.Subchunk2ID[0] << wavHeader.Subchunk2ID[1] << wavHeader.Subchunk2ID[2] << wavHeader.Subchunk2ID[3] << endl;

    // Ask the user for notch filter parameters
    cout << "Input bandwidth (0-0.5)    : ";
    cin >> input;
    cin.get();
    double BW = atof(input.c_str());

    double sampleFreq = wavHeader.SamplesPerSec;

    cout << "Input cut-off frequency ( < " << sampleFreq / 2.0 << "): ";
    cin >> input;
    cin.get();
    double cf = atof(input.c_str());

    notch_params notchParameters(BW, cf, sampleFreq);

    short int *filtered_i = new short int[numSamples];
    bandFilterOffline<short int, notch_params>(value_i, filtered_i, numSamples, notchParameters);

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

/**
 *  @brief This function computes the band-pass or band-reject filter offline given the input parameters.
 */
template <class T, class U>
void bandFilterOffline(const T *in, T *filtered, int n_samples, U params)
{
  double in_2 = 0.0;
  double in_1 = 0.0;
  double filtered_2 = 0.0;
  double filtered_1 = 0.0;

  for (int i = 0; i < n_samples; ++i)
  {
    filtered[i] = params.a0 * in[i] + params.a1 * in_1 + params.a2 * in_2 + params.b1 * filtered_1 + params.b2 * filtered_2;
    in_2 = in_1;
    in_1 = in[i];
    filtered_2 = filtered_1;
    filtered_1 = filtered[i];
  }
}

/**
 *  @brief This function computes the band-pass or band-reject filter online given the input parameters.
 */
template <class T, class U>
void bandFilterOnline(const T in, T &filtered, U params)
{
  static double in_2 = 0.0;
  static double in_1 = 0.0;
  static double filtered_2 = 0.0;
  static double filtered_1 = 0.0;

  filtered = params.a0 * in + params.a1 * in_1 + params.a2 * in_2 + params.b1 * filtered_1 + params.b2 * filtered_2;
  in_2 = in_1;
  in_1 = in;
  filtered_2 = filtered_1;
  filtered_1 = filtered;
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
