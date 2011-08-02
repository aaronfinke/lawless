// printing.hh

#ifndef PRINTING_HEADER
#define PRINTING_HEADER

// Routines for printing stuff

#include "aimless.hh"
#include "scalemodel.hh"
#include "score_datatypes.hh"
#include "halfdataset.hh"
#include "intensitybin.hh"
#include "summarystatistics.hh"

using namespace scala;

//--------------------------------------------------------------
void PrintTitle(phaser_io::Output& output);
//--------------------------------------------------------------
void PrintFileInfoToXML(const std::string& StreamName,
			const std::string& FileName,
			const Scell& cell,
			const std::string& SpaceGroupName,
			phaser_io::Output& output);
//--------------------------------------------------------------
void PrintOutlierSettings(const all_controls& controls, phaser_io::Output& output);
//--------------------------------------------------------------
// Print scale factors
void PrintScales(const ScaleModel& AllScales, phaser_io::Output& output);
//--------------------------------------------------------------
void PrintScalesByBatch(const PxdName& dataset_pxd,
			const std::vector<Batch>& batches, const std::vector<Run>& RunList,
			const int& datasetIndex,
			const std::vector<float>& scale0batch,
			const std::vector<float>& bfacbatch,
			const std::vector<MeanSD>& scalebatch,
			phaser_io::Output& output);
//--------------------------------------------------------------
void PrintDeviationsByBatch(const PxdName& dataset_pxd,
			    const std::vector<Batch>& batches,
			    const int& datasetIndex,
			    const std::vector<MeanSD>& imeanbatch,
			    const std::vector<MeanSD>& rmsDbatch,
			    const std::vector<Rfactor>& rmergebatch,
			    const std::vector<Rfactor>& rmergebatchsmoothed,
			    const std::vector<int>& rejectedbatch,
			    const std::vector<float>& batchcompleteness,
			    const std::vector<float>& batchanomcompleteness,
			    const std::vector<double>& maxresbatch,
			    const std::vector<double>& maxresbatchsmoothed,
			    const double& MinimumIoverSigma,
			    const int& nbatchsmooth,
			    const ResoRange& ResRange,
			    phaser_io::Output& output);
//--------------------------------------------------------------
void PrintDeviationsByResolution(const PxdName& dataset_pxd,
				 const ResoRange& ResRange,
				 const std::vector<Rfactor>& rmergeRes,
				 const std::vector<Rfactor>& rmeasRes,
				 const std::vector<Rfactor>& rpimRes,
				 const std::vector<MeanSD>&  imeanRes,
				 const std::vector<MeanSD>&  rmsDRes,
				 const std::vector<MeanSD>&  avSdRes,
				 const std::vector<MeanSD>&  mnIsdRes,
				 const std::vector<MeanSD>&  biasRes,
				 const std::vector<MeanSD>&  biasIRes,
				 const double& MinimumIoverSigma,
				 SummaryStatistics& summarystatistics,
				 phaser_io::Output& output);
//--------------------------------------------------------------
// Statistics against overall mean I+- (only if ANOMALOUS ON)
void PrintDeviationsByResolutionOv(const PxdName& dataset_pxd,
				   const ResoRange& ResRange,
				   const std::vector<Rfactor>& rmergeRes,
				   const std::vector<Rfactor>& rmeasRes,
				   const std::vector<Rfactor>& rpimRes,
				   const std::vector<Rfactor>& rmergeResOv,
				   const std::vector<Rfactor>& rmeasResOv,
				   const std::vector<Rfactor>& rpimResOv,
				   SummaryStatistics& summarystatistics,
				   phaser_io::Output& output);
//--------------------------------------------------------------
void PrintDeviationsByIntensity(const PxdName& dataset_pxd,
				 const IntensityBin& Irange,
				 const std::vector<Rfactor>& rmergeInt,
				 const std::vector<Rfactor>& rmeasInt,
				 const std::vector<Rfactor>& rpimInt,
				 const std::vector<MeanSD>&  imeanInt,
				 const std::vector<MeanSD>&  rmsDInt,
				 const std::vector<MeanSD>&  avSdInt,
				 const std::vector<MeanSD>&  mnIsdInt,
				 const std::vector<MeanSD>&  biasInt,
				 const std::vector<MeanSD>&  biasIInt,
				phaser_io::Output& output);
//--------------------------------------------------------------
void PrintCompletenessMultiplicity(const PxdName& dataset_pxd,
				   const ResoRange& ResRange,
				   const hkl_symmetry& symmetry,
				   const Scell& cell,
				   std::vector<int>& NumRef,
				   std::vector<int>& NumObs,
				   std::vector<int>& NumRefSphere,
				   std::vector<int>& NumCentric,
				   std::vector<int>& NumACentric,
				   std::vector<int>& NumAnom,
				   std::vector<int>& NumAnomSphere,
				   std::vector<double>& SNumAnomPairs,
				   SummaryStatistics& summarystatistics,
				   phaser_io::Output& output);
//--------------------------------------------------------------
void PrintHalfDatasetCorrelations(const PxdName& dataset_pxd,
				  const ResoRange& ResRange,
				  const HalfDataset& halfDatasetScores,
				  SummaryStatistics& summarystatistics,
				  phaser_io::Output& output);
//--------------------------------------------------------------
void PrintAnisotropyAnalysis(const PxdName& dataset_pxd,
			     const ResoRange& ResRange,
			     const HalfDataset& halfDatasetScores,
			     const std::vector<std::vector<MeanSD> >& mnIsdResCone,
			     const double& coneangledegrees,
			     const double& MinimumIoverSigma,
			     SummaryStatistics& summarystatistics,
			     phaser_io::Output& output);
//--------------------------------------------------------------
void PrintUnmergedHeaderStuff(const scala::hkl_unmerge_list& hkl_list,
			      phaser_io::Output& output,
			      const int& verbose);
//--------------------------------------------------------------
#endif
