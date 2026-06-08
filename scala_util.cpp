// scala_util.cpp
#include <cstdlib>
#include "scala_util.hh"
#include "string_util.hh"
#include "dataset.hh"

#include <assert.h>
#define ASSERT assert

namespace scala
{
//--------------------------------------------------------------
  std::string BaseFileName(const std::string& Name, const bool NoExt)
  // Extract base file name, ie strip off leading directory
  // paths and trailing extension (if NoExt true)
  {
    std::string bfn = Name.substr(Name.rfind('/')+1);
    if (NoExt) {
      return  bfn.substr(0,bfn.find('.'));
    } else {
      return bfn;
    }
  }
//--------------------------------------------------------------
  std::string FileNameNoExtension(const std::string& Name)
  // Strip off trailing extension from filename
  {
    // allow for "." preceding a "/"
    size_t jdot   = Name.rfind('.');   // find last '.' if any
    size_t jslash = Name.rfind('/');   // find last '/' if any
    if ((jslash != std::string::npos && jdot < jslash)
        || jdot == std::string::npos)
      // no extension
      {return Name;}
    else
      // strip extension
      {return Name.substr(0,jdot);}
  }
  //--------------------------------------------------------------
  std::string FileNameExtension(const std::string& name)
  // Return filename extension if present, excluding ".", else ""
  {
    if (name.size() > 0) {
      bool dot = false;
      int i = name.size();
      for (i=name.size();i>=0;i--) {
        if (name[i] == '/') break; // No dot found
        if (name[i] == '.') {
          dot = true;
          break;
        }
      }
      if (dot && i < int(name.size())-1) {
        return name.substr(i+1, name.size()-1-i);
      }
    }
    return "";
  }
  //--------------------------------------------------------------
  void AddFileExtension(std::string& name, const std::string& ext)
  // Add filename extension if not present
  // also remove leading spaces
  {
    // remove leading & trailing spaces
    name = StringUtil::Trim(name);
    if (name.size() > 0) {
      bool dot = false;
      for (int i=name.size();i>=0;i--) {
        if (name[i] == '/') break; // No dot found
        if (name[i] == '.') {
          dot = true;
          break;
        }
      }
      if (ext.size() > 0 && !dot) name += "."+ext;
    }
    return;
  }
//--------------------------------------------------------------
  std::string Directory(const std::string& Name)
// Return directory component of filename Name (if any)
  {
    if (Name.rfind('/') == std::string::npos) return "";
    return Name.substr(0, Name.rfind('/'));
  }
  //--------------------------------------------------------------
  std::string FormatCell(const std::vector<double>& cell,
                         const int w, const int p)
  {
    clipper::String line;
    for (int i=0;i<3;i++) line += StringUtil::ftos(cell[i],w,p);
    line += "   ";
    for (int i=0;i<3;i++) line += StringUtil::ftos(cell[i+3],w,p);
    return line;
  }
  //--------------------------------------------------------------
  std::string FormatCell(const Scell& cell, const int w, const int p)
  {
    clipper::String line;
    for (int i=0;i<3;i++) line += StringUtil::ftos(cell[i],w,p);
    line += "   ";
    for (int i=0;i<3;i++) line += StringUtil::ftos(cell[i+3],w,p);
    return line;
  }
  //--------------------------------------------------------------
  RPair MnSd(const std::vector<double>& val)
    // Return mean & SD of vector elements
  {
    int n = val.size();
    if (n < 2) return RPair(0.0,0.0);
    double sx = 0.0;
    double sxx = 0.0;
    double an = n;
    for (int i=0;i<n;i++)
      {
        sx += val[i];
        sxx += val[i]*val[i];
      }
    return RPair(sx/an, sqrt((an*sxx - sx*sx)/(an*(an-1.0))));
  }
  //--------------------------------------------------------------
  int IRandom(const int& MaxVal)
    // Return random integer between 0 and MaxVal
  {
    return Min(int((double(std::rand())/RAND_MAX) * MaxVal), MaxVal-1);
  }
  //--------------------------------------------------------------
  float FRandom(const float& MaxVal)
  // Return random float between 0 and MaxVal
  {
    return float(std::rand())/RAND_MAX * MaxVal;
  }
  //--------------------------------------------------------------
  double FRandom(const double& MaxVal)
  // Return random double between 0 and MaxVal
  {
    return double(std::rand())/RAND_MAX * MaxVal;
  }
  //--------------------------------------------------------------
  double SafeSqrt(const double& a)
    // returns sqrt(a) if a >=0, else - sqrt(-a)
  {
    return (a>=0.0) ? sqrt(a) : -sqrt(-a);
  }
  //--------------------------------------------------------------
  Scell AverageDsetCell(const std::vector<Dataset>& datasets)
  // Average unit cells over all datasets & store average
  // On entry:
  //  datasets     list of datasets
  // Returns:   average cell
  {
    UnitCellSet cellset;
    for (size_t k=0; k<datasets.size(); k++) {
      cellset.AddCellSet(datasets[k].AllCellSet());
    }
    return cellset.AverageCell();
  }
  //--------------------------------------------------------------
  std::vector<Scell> AverageBatchCell(const std::vector<Batch>& batches,
                                      const int& ndatasets,
                                      std::vector<double>& averageMosaicity,
                                      std::vector<double>& averageWavelength)
  // Average unit cells over all batches for each dataset
  // On entry:
  //  batches     list of batches
  // Returns:   average cell for each dataset
  //            averageMosaicity   average mosaicity for each dataset
  //            averageWavelength  average wavelength for each dataset
  {
    std::vector<Scell> averagecell(ndatasets);
    averageMosaicity.assign(ndatasets,0.0);
    averageWavelength.assign(ndatasets,0.0);
    int nbatches = batches.size();
    std::vector<int> n(ndatasets, 0);

    std::vector<std::vector<Scell> > allcells(ndatasets); // batch cell for each dataset

    for (int k=0; k<nbatches; k++)  {
      int idx = batches[k].datasetindex(); // dataset index
      allcells[idx].push_back(batches[k].cell()); // add batch cell
      n[idx]++;
      averageMosaicity[idx] += batches[k].Mosaicity();
      averageWavelength[idx] += batches[k].Wavelength();
    }
    for (int j=0;j<ndatasets;j++) {
      if (n[j] > 0) {
        averagecell[j] = UnitCellSet(allcells[j]).AverageCell();
        averageMosaicity[j] /= double(n[j]);
        averageWavelength[j]  /= double(n[j]);
        //^
        //std::cout << "AverageBatchCell: dataset "<<j<<"  "<<nbatches<<" batches "
        //                <<"\nAverage cell: "<< averagecell[j].format() <<"\n";
        //^-
      }
    }
    return averagecell;
  }
  //--------------------------------------------------------------
  Scell AverageBatchCellforDataset(const std::vector<Batch>& batches,
                                   const int& idts)
  // Average unit cells over all batches for specified dataset
  // On entry:
  //  batches     list of batches
  //  idts        dataset index
  // Returns:   average cell for dataset
  {
    int nbatches = batches.size();
    UnitCellSet cellset;

    for (int k=0; k<nbatches; k++)  {
      int idx = batches[k].datasetindex(); // dataset index
      if (idx == idts) {
        cellset.AddCell(batches[k].cell()); // add batch cell
      }
    }
    return cellset.AverageCell();
  }

  //--------------------------------------------------------------
  void AverageCells::init(const std::vector<Batch>& batches)
  // store batch cells, indexed by run
  {
    cellsets.clear();
    volumerange.clear();
    int nbatches = batches.size();
    for (int k=0; k<nbatches; k++)  {
      int idx = batches[k].RunIndex(); // run index
      cellsets.AddCell(batches[k].cell()); // add batch cell
      volumerange.update(batches[k].cell().Volume());
    }
  }
  //--------------------------------------------------------------
  //! return true if worst deviation is unacceptable
  bool AverageCells::badcell(const double& celldeviationlimit) const
  {
    return (WorstDeviation() > celldeviationlimit);
  }
  //--------------------------------------------------------------
  double AverageWavelength(const std::vector<double>& allwavelengths,
                    const int& idxexclude)
  // Average list of wavelengths
  // if idxexclude >= 0, exclude entry with this index
  {
    int nc = allwavelengths.size();
    if (idxexclude >= 0) nc--;
    if (nc <= 0) return 0.0;
    else if (nc == 1) return allwavelengths[0];
    // We have 2 or more,average
    double sumwavelength = 0.0;
    nc = 0;
    for (size_t k=0; k<allwavelengths.size(); k++) {
      if (idxexclude < 0 || int(k) != idxexclude) {
        if (allwavelengths[k] > 0.0) {
          nc++;
          sumwavelength += allwavelengths[k];
        }
      }
    }
    return double(sumwavelength/double(nc));
  }
  //--------------------------------------------------------------
  double AverageDsetWavelength(const std::vector<Dataset>& datasets)
  // Average wavelength over all datasets & store average
  // On entry:
  //  datasets     list of datasets
  // Returns:   average wavelength
  {
    double averagewvl;
    int ndatasets = datasets.size();
    double Sumwvl = 0.0;
    int n = 0;
    for (int k=0; k<ndatasets; k++) {
      std::vector<double> allwavelengths = datasets[k].AllWavelengths();
      for (size_t j=0; j<allwavelengths.size(); j++) {
        Sumwvl += allwavelengths[j];
        n++;
      }
    }
    averagewvl = Sumwvl/double(n);
    return averagewvl;
  }
  // ------------------------------------------------------------
  bool CompareIFpair(const std::pair<int,float>& p1,const std::pair<int,float>& p2)
  {return (p1.second < p2.second);}
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  MeanSD::MeanSD(const std::vector<float>& list)
    : sum_sc(0.0), sum_sc2(0.0), count(0)
  {
    for (size_t i=0;i<list.size();i++)  {Add(list[i]);}
  }
  //--------------------------------------------------------------
  MeanSD::MeanSD(const std::vector<double>& list) : sum_sc(0.0), sum_sc2(0.0), count(0)
  {
    for (size_t i=0;i<list.size();i++)  {Add(list[i]);}
  }
  //--------------------------------------------------------------
  void MeanSD::Add(const float& v)
  {
    sum_sc += v;
    sum_sc2 += v*v;
    count++;
  }
  //--------------------------------------------------------------
  void MeanSD::Add(const double& v)
  {
    sum_sc += v;
    sum_sc2 += v*v;
    count++;
  }
  //--------------------------------------------------------------
  double MeanSD::Mean() const
  {
    return (count > 0) ? sum_sc/double(count) : 0.0;
  }
  //--------------------------------------------------------------
  double MeanSD::SD() const
  {
    if (count <= 1) return 0.0;
    double a = Max(0.0, sum_sc2 - sum_sc*sum_sc/double(count));
    return sqrt(a/double(count-1));
  }
  //--------------------------------------------------------------
  double MeanSD::Var() const
  {
    return (count > 1) ?
      (sum_sc2 - sum_sc*sum_sc/double(count))/double(count-1) : 0.0;
  }
  //--------------------------------------------------------------
  std::string MeanSD::format() const
  {
    std::string s = "MeanSD: sumsc, count "+clipper::String(sum_sc)+
      +" "+clipper::String(count)+", "+
      "Mean "+clipper::String(Mean())+"  sd "+clipper::String(SD());
    return s;
    //    return "Mean "+clipper::String(Mean())+"  sd "+clipper::String(SD());
  }
  //--------------------------------------------------------------
  MeanSD& MeanSD::operator +=(const MeanSD& other)
  {
    sum_sc += other.sum_sc;
    sum_sc2 += other.sum_sc2;
    count += other.count;

    return *this;
  }
  //--------------------------------------------------------------
  MeanSD operator +
  (const MeanSD& a, const MeanSD& b)
  {
    MeanSD c = a;
    return c += b;
  }
  //--------------------------------------------------------------
  bool MeanSD::MeanSDsmallerSD(const MeanSD& a, const MeanSD& b)
  {
    return a.SD() < b.SD();
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  MeanValue::MeanValue(const std::vector<float>& list)
    : sum_sc(0.0), sum_w(0.0), count(0)
  {
    for (size_t i=0;i<list.size();i++)  {Add(list[i]);}
  }
  //--------------------------------------------------------------
  MeanValue::MeanValue(const std::vector<double>& list)
    : sum_sc(0.0), sum_w(0.0), count(0)
  {
    for (size_t i=0;i<list.size();i++)  {Add(list[i]);}
  }
  //--------------------------------------------------------------
  void MeanValue::Add(const float& v, const float& w)
  {
    sum_sc += w * v;
    sum_w += w;
    count++;
  }
  //--------------------------------------------------------------
  void MeanValue::Add(const double& v, const double& w)
  {
    sum_sc += w * v;
    sum_w += w;
    count++;
  }
  //--------------------------------------------------------------
  void MeanValue::Subtract(const double& v, const double& w)
  {
    sum_sc -= w * v;
    sum_w -= w;
    count--;
  }
  //--------------------------------------------------------------
  double MeanValue::Mean() const
  {
    return (sum_w > 0) ? sum_sc/sum_w : 0.0;
  }
  //--------------------------------------------------------------
  std::string MeanValue::format() const
  {
    std::string s = "MeanValue: sumsc, sumw "+clipper::String(sum_sc)+
      +" "+clipper::String(sum_w)+", "
      "Mean "+clipper::String(Mean());
    return s;
    //^^    return "Mean "+clipper::String(Mean());
  }
  //--------------------------------------------------------------
  MeanValue& MeanValue::operator +=(const MeanValue& other)
  {
    sum_sc += other.sum_sc;
    sum_w += other.sum_w;
    count += other.count;

    return *this;
  }
  //--------------------------------------------------------------
  MeanValue operator +
  (const MeanValue& a, const MeanValue& b)
  {
    MeanValue c = a;
    return c += b;
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  MeanVariance::MeanVariance(const std::vector<double>& list)
    : sum_sc(0.0), sum_w(0.0), sum_w2(0.0), sum_sc2(0.0), count(0)
  {
    for (size_t i=0;i<list.size();i++)  {Add(list[i]);}
  }
  //--------------------------------------------------------------
  MeanVariance::MeanVariance(const std::vector<float>& list)
    : sum_sc(0.0), sum_w(0.0), sum_w2(0.0), sum_sc2(0.0), count(0)
  {
    for (size_t i=0;i<list.size();i++)  {Add(list[i]);}
  }
  //--------------------------------------------------------------
  void MeanVariance::clear()
  {sum_sc=0.0; sum_w=0.0; sum_w2=0.0; sum_sc2 = 0.0; count=0;}
  //--------------------------------------------------------------
  void MeanVariance::Add(const float& v, const float& w)
  {
    sum_sc += w * v;
    sum_sc2 += w * v * v;
    sum_w += w;
    sum_w2 += w * w;
    count++;
  }
  //--------------------------------------------------------------
  void MeanVariance::Add(const double& v, const double& w)
  {
    sum_sc += w * v;
    sum_sc2 += w * v * v;
    sum_w += w;
    sum_w2 += w * w;
    count++;
  }
  //--------------------------------------------------------------
  void MeanVariance::Subtract(const double& v, const double& w)
  {
    sum_sc -= w * v;
    sum_sc2 -= w * v * v;
    sum_w -= w;
    sum_w2 -= w * w;
    count--;
  }
  //--------------------------------------------------------------
  double MeanVariance::Mean() const
  {
    return (sum_w > 0) ? sum_sc/sum_w : 0.0;
  }
  //--------------------------------------------------------------
  double MeanVariance::VarianceofMeanFromWeights() const
  // variance of mean from the weights ie 1/Sum(weights)
  {
    double var = 0.0;
    if (count > 0) {
      var = 1.0 / sum_w;
    }
    return var;
  }
  //--------------------------------------------------------------
  double MeanVariance::SDofMeanfromWeights() const
  // SD of mean from the weights ie sqrt(1/Sum(weights))
  {
    return sqrt(VarianceofMeanFromWeights());
  }
  //--------------------------------------------------------------
  double MeanVariance::VarianceofMean() const
  // variance of sample
  {
    double var = 0.0;
    if (count > 1) {
      // Variance of mean is just 1/N sampleVariance
      //  assuming the weights are proportional to the true variances
      var = SampleVariance()/double(count);
    }
  return var;
  }
  //--------------------------------------------------------------
  double MeanVariance::SDofMean() const
  // SD of mean
  {
    return sqrt(VarianceofMean());
  }
  //--------------------------------------------------------------
  double MeanVariance::SampleVariance() const
  // Variance of sample from the distribution,= 0.0 if < 2 values
  {
    double var = 0.0;
    if (count > 1) {
      // fac = Sum(w)/[(Sum(w))^2 - Sum(w^2)] to allow for bias
      //     = 1/(sum_w-(sumw^2/sum_w))
      //   equivalent to n/(n-1) correction
      //     for unit weights, fac = 1/(n-1)
      //     see Wikipedia Weighted Mean,  1/(V1 - V2/V1)
      double fac = 1.0/(sum_w - sum_w2/sum_w);
      var = Max(0.0, sum_sc2 - sum_sc*sum_sc/sum_w) * fac;
      //^
      //      std::cout <<"SampleVariance "<<count<<" "<<1.0/(count-1.0)<<" "
      //                <<var<<" "<<fac<<" "<<sum_sc2*fac<<" "<<
      //        sum_sc/sum_w <<" "<<(sum_sc/sum_w)*(sum_sc/sum_w)
      //                <<" "<<sum_sc2/(count-1.0) - (sum_sc/sum_w)*(sum_sc/sum_w)
      //        <<"\n";
    }
    return var;
  }
  //--------------------------------------------------------------
  double MeanVariance::SampleVarianceOmit1(const double& v,
                                     const double& w) const
  // Variance of mean, omitting one observation (v, weight w)
  /*  from add
    sum_sc += w * v;
    sum_sc2 += w * v * v;
    sum_w += w;
    sum_w2 += w * w;
    count++;
  */
  {
    double var = 0.0;
    double sum_w_less1 = sum_w - w;
    if (sum_w_less1 <= 0.0) {return var;}
    if (count > 2) {
      // fac = Sum(w)/[(Sum(w))^2 - Sum(w^2)] to allow for bias
      //     = 1/(sum_w-(sumw^2/sum(w)))
      //   equivalent to n/(n-1) correction
      //     see Wikipedia Weighted Mean,  1/(V1 - V2/V1)
      double fac =
        1.0/(sum_w_less1 - (sum_w2 - w*w)/sum_w_less1);
      double sum_sc_less1 = sum_sc - w * v;
      var = Max(0.0,
                (sum_sc2 - w*v*v) -
                sum_sc_less1*sum_sc_less1/sum_w_less1) * fac;
      //      std::cout <<"Omit1 "<<v<<" "<<w<<" "<<var<<" "<<sum_sc<<" "<<
      //        sum_sc_less1<<" "<<sum_sc2-w*v*v<<
      //        sum_w_less1<<" "<<sum_w2-w*w<<std::endl; //^-
    } else if (count == 2) {
      var = 1.0/sum_w_less1;
    }
    return var;
  }
  //--------------------------------------------------------------
  double MeanVariance::VarianceofMeanOmit1(const double& v,
                                     const double& w) const
  {
    double var = 0.0;
    if (count > 2) {
      // Variance of mean is just 1/N sampleVariance
      //  assuming the weights are proportional to the true variances
      var = SampleVarianceOmit1(v,w)/double(count-1);
    }
    return var;
  }
  //--------------------------------------------------------------
  double MeanVariance::SampleSD() const
  {
    return sqrt(SampleVariance());
  }
  //--------------------------------------------------------------
  std::string MeanVariance::format() const
  {
    std::string s = "Mean "+clipper::String(Mean()) +
      ", SD "+ clipper::String(SampleSD())+
      ", SD(weights) " + clipper::String(sqrt(count)*SDofMeanfromWeights());
    return s;
  }
  //--------------------------------------------------------------
  MeanVariance& MeanVariance::operator +=(const MeanVariance& other)
  {
    sum_sc += other.sum_sc;
    sum_w += other.sum_w;
    sum_w2 += other.sum_w2;
    sum_sc2 += other.sum_sc2;
    count += other.count;

    return *this;
  }
  //--------------------------------------------------------------
  MeanVariance operator +
  (const MeanVariance& a, const MeanVariance& b)
  {
    MeanVariance c = a;
    return c += b;
  }
  //--------------------------------------------------------------
  MeanVariance& MeanVariance::operator -=(const MeanVariance& other)
  {
    sum_sc -= other.sum_sc;
    sum_w -= other.sum_w;
    sum_w2 -= other.sum_w2;
    sum_sc2 -= other.sum_sc2;
    count -= other.count;

    return *this;
  }
  //--------------------------------------------------------------
  MeanVariance operator -
  (const MeanVariance& a, const MeanVariance& b)
  {
    MeanVariance c = a;
    return c -= b;
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
}
