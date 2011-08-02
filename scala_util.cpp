// scala_util.cpp
#include <cstdlib>
#include "scala_util.hh"
#include "string_util.hh"

#define ASSERT assert
#include <assert.h>

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
    for (int i=0;i<3;i++) line += clipper::String(cell[i],w,p);
    line += "   ";
    for (int i=0;i<3;i++) line += clipper::String(cell[i+3],w,p);
    return line;
  } 
  //--------------------------------------------------------------
  std::string FormatCell(const Scell& cell, const int w, const int p)
  {
    clipper::String line;
    for (int i=0;i<3;i++) line += clipper::String(cell[i],w,p);
    line += "   ";
    for (int i=0;i<3;i++) line += clipper::String(cell[i+3],w,p);
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
  Scell AverageDsetCell(const std::vector<Xdataset>& datasets)
  // Average unit cells over all datasets & store average
  // On entry:
  //  datasets     list of datasets
  // Returns:   average cell
  {
    Scell averagecell;
    int ndatasets = datasets.size();
    if (ndatasets <= 1) {
      averagecell = datasets[0].cell();
    } else {
      std::vector<Scell> allcells;
      for (int k=0; k<ndatasets; k++) {
	allcells.push_back(datasets[k].cell());
      }
      averagecell = UnitCellSet(allcells).AverageCell();
    }
    return averagecell;
  }
  //--------------------------------------------------------------
  std::vector<Scell> AverageBatchCell(const std::vector<Batch>& batches,
				      const int& ndatasets,
				      std::vector<float>& averageMosaicity,
				      std::vector<float>& averageWavelength)
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
      int idx = batches[k].index(); // dataset index
      allcells[idx].push_back(batches[k].cell()); // add batch cell
      n[idx]++;
      averageMosaicity[idx] += batches[k].Mosaicity(); 
      averageWavelength[idx] += batches[k].Wavelength(); 
    }
    for (int j=0;j<ndatasets;j++) {
      if (n[j] > 0) {
	averagecell[j] = UnitCellSet(allcells[j]).AverageCell();
	averageMosaicity[j] /= float(n[j]);
	averageWavelength[j]  /= float(n[j]);
      }
    }
    return averagecell;
  }
  //--------------------------------------------------------------
  float AverageWavelength(const std::vector<float>& allwavelengths,
		    const int& idxexclude)
  // Average list of wavelengths
  // if idxexclude >= 0, exclude entry with this index
  {
    int nc = allwavelengths.size();
    if (idxexclude >= 0) nc--;
    if (nc <= 0) return 0.0;
    else if (nc == 1) return allwavelengths[0];
    // We have 2 or more,average
    float sumwavelength = 0.0;
    for (size_t k=0; k<allwavelengths.size(); k++) {
      if (idxexclude < 0 || int(k) != idxexclude) {
	sumwavelength += allwavelengths[k];
      }
    }
    return sumwavelength/nc;
  }
  //--------------------------------------------------------------
  float AverageDsetWavelength(const std::vector<Xdataset>& datasets)
  // Average wavelength over all datasets & store average
  // On entry:
  //  datasets     list of datasets
  // Returns:   average wavelength
  {
    float averagewvl;
    int ndatasets = datasets.size();
    if (ndatasets <= 1) {
      averagewvl = datasets[0].wavelength();
    } else {
      double Sumwvl = 0.0;
      for (int k=0; k<ndatasets; k++) {
	Sumwvl += datasets[k].wavelength();
      }
      averagewvl = Sumwvl/float(ndatasets);
    }
    return averagewvl;
  }
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
  void MeanSD::initExclude(const std::vector<float>& list,
			   const int& idxexclude) 
  // if idxexclude >= 0, exclude entry with this index
  {
    for (size_t i=0;i<list.size();i++)  {
      if (idxexclude < 0 || int(i) != idxexclude) Add(list[i]);
    }
  }
  //--------------------------------------------------------------
  void MeanSD::initExclude(const std::vector<double>& list,
			   const int& idxexclude) 
  // if idxexclude >= 0, exclude entry with this index
  {
    for (size_t i=0;i<list.size();i++)  {
      if (idxexclude < 0 || int(i) != idxexclude) Add(list[i]);
    }
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
    return (count > 1) ?
      sqrt((sum_sc2 - sum_sc*sum_sc/double(count))/double(count-1)) : 0.0;
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
    return "Mean "+clipper::String(Mean())+"  sd "+clipper::String(SD()); 
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
  MeanSD& operator +
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
    : sum_sc(0.0), count(0)
  {
    for (size_t i=0;i<list.size();i++)  {Add(list[i]);}
  }
  //--------------------------------------------------------------
  MeanValue::MeanValue(const std::vector<double>& list) : sum_sc(0.0), count(0)
  {
    for (size_t i=0;i<list.size();i++)  {Add(list[i]);}
  }
  //--------------------------------------------------------------
  void MeanValue::Add(const float& v)
  {
    sum_sc += v;
    count++;
  }
  //--------------------------------------------------------------
  void MeanValue::Add(const double& v)
  {
    sum_sc += v;
    count++;
  }
  //--------------------------------------------------------------
  double MeanValue::Mean() const
  {
    return (count > 0) ? sum_sc/double(count) : 0.0;
  }
  //--------------------------------------------------------------
  std::string MeanValue::format() const
  {
    return "Mean "+clipper::String(Mean()); 
  }
  //--------------------------------------------------------------
  MeanValue& MeanValue::operator +=(const MeanValue& other)
  {
    sum_sc += other.sum_sc;
    count += other.count;
  
    return *this;
  }
  //--------------------------------------------------------------
  MeanValue& operator +
  (const MeanValue& a, const MeanValue& b)
  {
    MeanValue c = a;
    return c += b;
  }
  //--------------------------------------------------------------
}


