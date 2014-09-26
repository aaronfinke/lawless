#ifndef SCALA_UTIL
#define SCALA_UTIL

#include <string>
#include "hkl_datatypes.hh"

namespace scala
{
  //--------------------------------------------------------------
  // Extract base file name, ie strip off leading directory
  // paths and trailing extension (if NoExt true)
  std::string BaseFileName(const std::string& Name, const bool NoExt=true);
  //--------------------------------------------------------------
  std::string FileNameNoExtension(const std::string& Name);
  //--------------------------------------------------------------
  // Return filename extension if present, excluding ".", else ""
  std::string FileNameExtension(const std::string& name);
  //--------------------------------------------------------------
  void AddFileExtension(std::string& name, const std::string& ext);
  //--------------------------------------------------------------
  // Return directory component of filename Name (if any)
  std::string Directory(const std::string& Name);
  //--------------------------------------------------------------
  // Return formatted unit cell
  std::string FormatCell(const std::vector<double>& cell,
			 const int w=7, const int p=2);
  //--------------------------------------------------------------
  // Return formatted unit cell
  std::string FormatCell(const Scell& cell, const int w=7, const int p=2);
  //--------------------------------------------------------------
  RPair MnSd(const std::vector<double>& val);
  // Return mean & SD of vector elements
  //--------------------------------------------------------------
  int IRandom(const int& MaxVal);
    // Return random integer between 0 and MaxVal
  //--------------------------------------------------------------
  float FRandom(const float& MaxVal);
  // Return random float between 0 and MaxVal
  //--------------------------------------------------------------
  double FRandom(const double& MaxVal);
  // Return random double between 0 and MaxVal
  // //--------------------------------------------------------------
  // returns sqrt(a) if a >=0, else - sqrt(-a)
  double SafeSqrt(const double& a);
  //--------------------------------------------------------------
  // Average unit cells over all datasets & store average
  // On entry:
  //  datasets     list of datasets
  // Returns:   average cell
  Scell AverageDsetCell(const std::vector<Dataset>& datasets);
  //--------------------------------------------------------------
  // Average unit cells over all batches for each dataset
  // On entry:
  //  batches     list of batches
  // Returns:   average cell for each dataset
  //            averageMosaicity   average mosaicity  for each dataset
  //            averageWavelength  average wavelength for each dataset
  std::vector<Scell> AverageBatchCell(const std::vector<Batch>& batches,
				      const int& Ndatasets,
				      std::vector<float>& averageMosaicity,
				      std::vector<float>& averageWavelength);
  //--------------------------------------------------------------
  // Average wavelengths over all datasets & store average
  // On entry:
  //  datasets     list of datasets
  // Returns:   average wavelength
  float AverageDsetWavelength(const std::vector<Dataset>& datasets);
  //--------------------------------------------------------------
  // Average list of wavelengths
  // if idxexclude >= 0, exclude entry with this index
  float AverageWavelength(const std::vector<float>& allwavelengths,
			  const int& idxexclude=-1);
  //======================================================================
  //======================================================================
  class MeanSD
  // Mean & SD of numbers (SD of distribution, not of mean), unit weights
  {
  public:
    MeanSD() : sum_sc(0.0), sum_sc2(0.0), count(0) {}
    MeanSD(const std::vector<float>& list);
    MeanSD(const std::vector<double>& list);

    // if idxexclude >= 0, exclude entry with this index
    void initExclude(const std::vector<float>& list, const int& idxexclude=-1);
    void initExclude(const std::vector<double>& list, const int& idxexclude=-1);

    void Add(const float& v);
    void Add(const double& v);
    void clear() {sum_sc=0.0; sum_sc2=0.0; count=0;}

    double Mean() const;
    double SD() const;
    double Var() const;  // variance
    int Count() const {return count;}

    std::string format() const;

    MeanSD& operator+=(const MeanSD& other);
    friend MeanSD& operator+ (const MeanSD& a, const MeanSD& b);
    static bool MeanSDsmallerSD(const MeanSD& a, const MeanSD& b);

  private:
    double sum_sc;
    double sum_sc2;
    int count;
  };
  //======================================================================
  class MeanValue
  // Mean  of numbers
  {
  public:
    MeanValue() : sum_sc(0.0), sum_w(0.0), count(0) {}
    MeanValue(const std::vector<float>& list);
    MeanValue(const std::vector<double>& list);
  
    void Add(const float& v, const float& w = 1.0f);
    void Add(const double& v, const double& w = 1.0);
    void clear() {sum_sc=0.0; sum_w=0.0; count=0;}

    double Mean() const;
    int Count() const {return count;}

    std::string format() const;

    MeanValue& operator+=(const MeanValue& other);
    friend MeanValue& operator+ (const MeanValue& a, const MeanValue& b);

  private:
    double sum_sc, sum_w;
    int count;
  };
  //======================================================================
  class MeanVariance
  // Weighted mean & variance of mean
  {
  public:
    MeanVariance() : sum_sc(0.0), sum_w(0.0), sum_w2(0.0),
		     sum_sc2(0.0), count(0) {}
    MeanVariance(const std::vector<float>& list);
    MeanVariance(const std::vector<double>& list);
  
    void Add(const float& v, const float& w = 1.0f);
    void Add(const double& v, const double& w = 1.0);
    void clear();

    double Mean() const;

    // variance of mean from the weights ie 1/Sum(weights)
    double VarianceFromWeights() const;
    // SD of mean from the weights ie sqrt(1/Sum(weights))
    double SDfromWeights() const;

    // Variance of mean
    double Variance() const;
    double SD() const;



    // Variance of distribution
    double SampleVariance() const;
    double SampleSD() const;

    int Count() const {return count;}

    std::string format() const;

    MeanVariance& operator+=(const MeanVariance& other);
    friend MeanVariance& operator+ (const MeanVariance& a, const MeanVariance& b);

  private:
    double sum_sc, sum_w, sum_w2;
    double sum_sc2;
    int count;
  };

}

#endif

