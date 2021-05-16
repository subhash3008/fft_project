#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <complex>
#include <valarray>


const double PI = 3.141592653589793238460;


std::vector<double> real;
std::vector<double> img;


// Read the data from file
void readData() {
  std::ifstream rawData { "alice.txt" };
  if (!rawData) {
    std::cout << "Error reading alice.txt\n";
    return;
  }

  long rawDataIndex = 0;
  while(!rawData.eof()) {
    std::string val;
    std::getline(rawData, val);
    if (rawDataIndex % 2 == 0) {
      real.push_back(std::stod(val));
    } else {
      img.push_back(std::stod(val));
    }
    ++rawDataIndex;
  }
}

// Write the plot file
void generatePlot(const std::string& fileName) {
  std::ofstream plotData { fileName };
  if (!plotData) {
    std::cout << "Cannot open plotData file\n";
  }
  for (int x = 0; x < real.size(); ++x) {
    plotData << x << " " << real.at(x) << std::endl;
  }
}

// Print the Data
void printData() {
  for (int x = 0; x < real.size(); ++x) {
    std::cout << x << " " << real.at(x) << "+ i" << img.at(x) << std::endl;
  }
}

void reverseBits() {
  long datasetSize = real.size();
  long i2 = datasetSize >> 1;
  long j = 0;
  for (long i = 0; i < datasetSize - 1; ++i) {
    if (i < j) {
      double tempReal = real.at(i);
      double tempImg = img.at(i);
      real[i] = real[j];
      img[i] = img[j];
      real[j] = tempReal;
      img[j] = tempImg;
    }
    long k = i2;
    while (k <= j) {
      j = j - k;
      k = k >> 1;
    }
    j += k;
  }
}


// Danielson-Lanzcos routine
void fft(long nearestPowerOf2Num, int isInverseFFT) {
  long mmax = 2;

  while (nearestPowerOf2Num > mmax) {
    long istep = mmax << 1;
    double theta = isInverseFFT * (2 * PI) / mmax;
    double wtemp = std::sin(.5 * theta);
    double wpr = -2.0 * wtemp * wtemp;
    double wpi = std::sin(theta);
    double wr = 1.0;
    double wi = 0.0;
    for (long m = 1; m < mmax; m+=2) {
      for (long i = m; i < real.size(); i += istep) {
        long j = i + mmax;
        double tempr = wr * real.at(j) - wi * img.at(j);
        double tempi = wr * img.at(j) + wi * real.at(j);
        real.at(j) = real.at(i) - tempr;
        img.at(j) = img.at(i) - tempi;
        real.at(i) += tempr;
        img.at(i) += tempi;
      }
      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = (wi * wpr) + wtemp * wpi + wi;
    }
    mmax = istep;
  }
}

int main() {
  readData();

  generatePlot("inputSignal.dat");

  // for (int gy = 0; gy < 8; ++gy) {
  //   if (gy < 4) {
  //     real.push_back(1.0);
  //   } else {
  //     real.push_back(0.0);
  //   }
  //   img.push_back(0.0);
  // }

  int nearestPowerOf2 = 1;
  long nearestPowerOf2Num = 1;
  while (nearestPowerOf2Num < real.size()) {
    nearestPowerOf2Num = nearestPowerOf2Num << 1;
    ++nearestPowerOf2;
  }


  // Padding 0
  long loopCount = nearestPowerOf2Num - real.size();
  for (long i = 0; i < loopCount; ++i) {
    real.push_back(0.0);
    img.push_back(0.0);
  }
  reverseBits();

  fft(nearestPowerOf2Num, 1);

  generatePlot("plotFFT.dat");

  fft(nearestPowerOf2Num, -1);
  generatePlot("inverseFFT.dat");

  return 0;
}