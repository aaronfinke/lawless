// Clipper
#include <clipper/clipper.h>
using clipper::Message;
using clipper::Message_fatal;
using clipper::Message_warn;

#include "dataset.hh"
#include "string_util.hh"

namespace scala {
  //--------------------------------------------------------------
  //! construct from one Xdataset
  Dataset::Dataset(const Xdataset& xdataset)
  {
    xdatasets.clear();
    AddXdataset(xdataset);
  }
  //--------------------------------------------------------------
  //! Add in a new xdataset if dataset name is the same, return true if added
  bool Dataset::AddXdataset(const Xdataset& xdataset)
  {
    bool added = false;
    // Do we have this one already?
    if (xdatasets.size() > 0) {
      for (size_t i=0;i<xdatasets.size();++i) {
        if (xdataset.pxdname().dname() == xdatasets[i].pxdname().dname()) {
          added = true; // same Dname, so add it
          break;
        }
      }
    } else {
      added = true;  // add anyway if list is empty
    }
    if (added) {
      xdatasets.push_back(xdataset);
    }
    check(); // belt, braces and piece of string
    return added;
  }
  //--------------------------------------------------------------
  //! Add in all new xdatasets if only the dataset name is the same,
  //!  return true if added
  bool Dataset::AddDataset(const Dataset& other)
  {
    bool added = false;
    if (*this == other) {
      // don't do anything if they are identical, do it elsewhere
      return added;
    }
    if (Dname() == other.Dname()) {
      // same dataset, Project or Crystal different
      // add all Xdatasets from other
      for (size_t i=0;i<other.xdatasets.size();++i) {
        xdatasets.push_back(other.xdatasets[i]);
        added = true;
      }
    }
    check(); // belt, braces and piece of string
    return added;
  }
  //--------------------------------------------------------------
  //! Add batch number to list for this dataset ID
  void Dataset::add_batch(const int& setid, const int& batch_num)
  {
    int idx = XdatasetIndex(setid);
    if (idx < 0) { // not found
      Message::message(Message_fatal
                       ("Dataset::add_batch, setid not found "+clipper::String(setid)));
    }
    xdatasets[idx].add_batch(batch_num);  // add it
  }
  //--------------------------------------------------------------
  //! Add batch number to list for this PxdName
  void Dataset::add_batch(const PxdName& pxdname, const int& batch_num)
  {
    int idx = XdatasetIndex(pxdname);
    if (idx < 0) { // not found
      Message::message(Message_fatal
                       ("Dataset::add_batch, pxdname not found "+pxdname.format()));
    }
    xdatasets[idx].add_batch(batch_num);  // add it
  }
  //--------------------------------------------------------------
  //! return list of PXDnames
  std::vector<PxdName> Dataset::pxdnames() const //!< return all PXD names
  {
    std::vector<PxdName> vpxdnames;
    for (size_t i=0;i<xdatasets.size();++i) {
      vpxdnames.push_back(xdatasets[i].pxdname());
    }
    return vpxdnames;
  }
  //--------------------------------------------------------------
  //! return PXDname for given set ID
  PxdName Dataset::pxdname(const int& setid) const
  {
    int idx = XdatasetIndex(setid);
    if (idx < 0) { // not found
      Message::message(Message_fatal
                       ("Dataset::pxdname, setid not found "+clipper::String(setid)));
    }
    return xdatasets[idx].pxdname();
  }
  //--------------------------------------------------------------
  //! return consensus PXDname (just set Dname to "MultiCrystal")
  PxdName Dataset::pxdname() const
  {
    dieIfEmpty("pxdname");
    if (xdatasets.size() == 1) { //only one
      return xdatasets[0].pxdname();
    }
    // Pname and Dname should all be the same (Pname doesn't matter)
    std::string pname = xdatasets[0].pxdname().pname(); // use 1st Pname
    std::string dname = xdatasets[0].pxdname().dname(); // use 1st Dname
    //
    std::string xname = "MultiCrystal"; // can't see how to be clever!
    return PxdName(pname, xname, dname);
  }
  //--------------------------------------------------------------
  //! return project name (all Xdatasets have same project name)
  std::string Dataset::Pname() const
  {
    dieIfEmpty("Pname");
    return xdatasets[0].pxdname().pname();
  }
  //--------------------------------------------------------------
  //! return dataset name (all Xdatasets have same dataset name)
  std::string Dataset::Dname() const
  {
    dieIfEmpty("Dname");
    return xdatasets[0].pxdname().dname();
  }
  //--------------------------------------------------------------
  UnitCellSet Dataset::AllCellSet() const
  {
    UnitCellSet cellset;
    for (size_t k=0; k<xdatasets.size(); k++) {
      cellset.AddCellSet(xdatasets[k].AllCells());
    }
    return cellset;
  }
  //--------------------------------------------------------------
  std::vector<Scell> Dataset::AllCells() const
  {
    return AllCellSet().Cells();
  }
  //--------------------------------------------------------------
  //! return average (or sole) cell
  Scell Dataset::cell() const
  {
    dieIfEmpty("cell");
    if (xdatasets.size() == 1) { //only one
      return xdatasets[0].cell();
    }
    // Average cells
    return AllCellSet().AverageCell();
  }
  //--------------------------------------------------------------
  std::vector<float> Dataset::AllWavelengths() const
  {
    std::vector<float> allwavelengths;
    for (size_t k=0; k<xdatasets.size(); k++) {
      std::vector<float> wvl = xdatasets[k].AllWavelengths();
      allwavelengths.insert(allwavelengths.end(),
                            wvl.begin(), wvl.end());
    }
    return allwavelengths;
  }
  //--------------------------------------------------------------
  //! return average (or sole) wavelength
  float Dataset::wavelength() const
  {
    dieIfEmpty("wavelength");
    if (xdatasets.size() == 1) { //only one
      return xdatasets[0].wavelength();
    }
    // We have 2 or more,average
    double sumwavelength = 0.0;
    std::vector<float> allwavelengths = AllWavelengths();
    for (size_t k=0; k<allwavelengths.size(); k++) {
      sumwavelength += allwavelengths[k];
    }
    return float(sumwavelength/double(allwavelengths.size()));
  }
  //--------------------------------------------------------------
  //! return wavelength for named crystal
  float Dataset::wavelength(const std::string& xname) const
  {
    dieIfEmpty("wavelength named");
    int idx = XdatasetIndex(xname);
    if (idx < 0) {
      Message::message(Message_fatal
                       ("Dataset::wavelength, xname not found "+xname));
    }
    return xdatasets[idx].wavelength();
  }
  //--------------------------------------------------------------
  //! return wavelength range
  Range Dataset::wavelengthRange() const
  {
    dieIfEmpty("wavelength");
    if (xdatasets.size() == 1) { //only one
      return Range(xdatasets[0].wavelength(), xdatasets[0].wavelength());
    }
    // We have 2 or more,average
    Range wvlr;
    for (size_t k=0; k<xdatasets.size(); k++) {
      wvlr.update(xdatasets[k].wavelength());
    }
    return wvlr;
  }
  //--------------------------------------------------------------
  //!< return average mosaicity
  float Dataset::Mosaicity() const
  {
    dieIfEmpty("Mosaicity");
    if (xdatasets.size() == 1) { //only one
      return xdatasets[0].Mosaicity();
    }
    double summos = 0.0;
    for (size_t k=0; k<xdatasets.size(); k++) {
      summos += xdatasets[k].Mosaicity();
    }
    return float(summos/double(xdatasets.size()));
  }
  //--------------------------------------------------------------
  void Dataset::SetResRange(const ResoRange& resrange) //!< set resolution range
  // Set to minimum range of current range & requested range
  {
    dieIfEmpty("SetResRange");
    for (size_t k=0; k<xdatasets.size(); k++) {
      ResoRange resrangenow = xdatasets[k].ResRange();
      xdatasets[k].ResRange() = resrangenow.MinRange(resrange);
    }
  }
  //--------------------------------------------------------------
  ResoRange Dataset::ResRange() const   //!< return maximum resolution range
  {
    dieIfEmpty("ResRange");
    ResoRange maxresrange = xdatasets[0].ResRange();
    if (xdatasets.size() > 1) {
      for (size_t k=1; k<xdatasets.size(); k++) { // loop from 2nd
        maxresrange = maxresrange.MaxRange(xdatasets[k].ResRange());
      }
    }
    return maxresrange;
  }
  //--------------------------------------------------------------
  //! change basis: reindex to get new cell
  void Dataset::change_basis(const ReindexOp& reindex_op)
  {
    dieIfEmpty("change_basis");
    for (size_t k=0; k<xdatasets.size(); k++) {
      xdatasets[k].change_basis(reindex_op);
    }
  }
  //--------------------------------------------------------------
  //! Set unit cells for all Xdatasets
  void Dataset::SetCellWavelength(const Scell& cell, const float& wavel)
  ///  void Dataset::SetCell(const Scell& cell)
  {
    dieIfEmpty("SetCell");
    for (size_t k=0; k<xdatasets.size(); k++) {
      xdatasets[k].SetCellWavelength(cell, wavel);
    }
  }
  //--------------------------------------------------------------
  //! Set mosaicity for all Xdatasets
  void Dataset::SetMosaicity(const float& mosaicity)
  {
    dieIfEmpty("SetMosaicity");
    for (size_t k=0; k<xdatasets.size(); k++) {
      xdatasets[k].Mosaicity() = mosaicity;
    }
  }
  //--------------------------------------------------------------
  //! add in another cell and wavelength, put into 1st Xdataset
  void Dataset::AddCellWavelength(const Scell& newcell, const float& wavel)
  {
    dieIfEmpty("AddCellWavelength");
    xdatasets[0].AddCellWavelength(newcell, wavel);
  }
  //--------------------------------------------------------------
  //! add to run index list for given Xdataset
  void Dataset::AddRunIndex(const PxdName& pxdname, const int& RunIndex)
  {
    int idx = XdatasetIndex(pxdname);
    if (idx < 0) { // not found
      Message::message(Message_fatal
                       ("Dataset::AddRunIndex, pxdname not found "+pxdname.format()));
    }
    xdatasets[idx].AddRunIndex(RunIndex);  // add it
  }
  //--------------------------------------------------------------
  //! clear all runs
  void Dataset::ClearRunList()
  {
    dieIfEmpty("ClearRunList");
    for (size_t k=0; k<xdatasets.size(); k++) {
      xdatasets[k].ClearRunList();
    }
  }
  //--------------------------------------------------------------
  //!< return run index list for all xdatasets
  std::vector<int> Dataset::RunIndexList() const
  {
    dieIfEmpty("RunIndexList");
    if (xdatasets.size() == 1) { //only one
      return xdatasets[0].RunIndexList();
    }
    // We have 2 or more, concatenate
    std::vector<int> runindexlist;
    for (size_t k=0; k<xdatasets.size(); k++) {
      std::vector<int> ril = xdatasets[k].RunIndexList();
      runindexlist.insert(runindexlist.end(),
                          ril.begin(), ril.end());
    }
    return runindexlist;
  }
  //--------------------------------------------------------------
  //! return list of datasetIDs
  std::vector<int> Dataset::SetIDs() const
  {
    dieIfEmpty("SetIDs");
    std::vector<int> setids;
    for (size_t k=0; k<xdatasets.size(); k++) {
      setids.push_back(xdatasets[k].setid());
    }
    return setids;
  }
  //--------------------------------------------------------------
  //! return datasetID for given PXDname
  int Dataset::GetID(const PxdName& pxdname) const
  {
    dieIfEmpty("GetID");
    int idx = XdatasetIndex(pxdname);
    if (idx < 0) { // not found
      Message::message(Message_fatal
                       ("Dataset::GetID, pxdname not found "+pxdname.format()));
    }
    return xdatasets[idx].setid();
  }
  //--------------------------------------------------------------
  //! set datasetIDs to ID, ID+n etc, return last value used
  int Dataset::StoreSetID(const int& ID)
  {
    dieIfEmpty("StoreSetID");
    int id = ID - 1;
    for (size_t k=0; k<xdatasets.size(); k++) {
      xdatasets[k].setid() = ++id;
    }
    return id;
  }
  //--------------------------------------------------------------
  //! return SetID >=0 if xdataset with given PXD name is present here, else -1
  bool Dataset::IsPxdPresent(const PxdName& pxdname) const
  {
    int setid = XdatasetIndex(pxdname);
    if (setid >= 0) { // found,
      return true;
    }
    return false;
  }
  //--------------------------------------------------------------
  //! return true if xdataset with given ID is present here
  bool Dataset::IsSetidPresent(const int& setid) const
  {
    int idx = XdatasetIndex(setid);
    if (idx >= 0) { // found,
      return true;
    }
    return false;
  }
  //--------------------------------------------------------------
  void Dataset::ClearBatchList()   //!< clear batch lists
  {
    dieIfEmpty("ClearBatchList");
    for (size_t k=0; k<xdatasets.size(); k++) {
      xdatasets[k].ClearBatchList();
    }
  }
  //--------------------------------------------------------------
  int Dataset::num_batches() const //!< number of batches in all Xdatasets
  {
    dieIfEmpty("num_batches");
    int nb = 0;
    for (size_t k=0; k<xdatasets.size(); k++) {
      nb += xdatasets[k].num_batches();
    }
    return nb;
  }
  //--------------------------------------------------------------
  //! return maximum ID of datasets here or given ID (if >0)
  int Dataset::MaxID(const int& ID) const
  {
    dieIfEmpty("MaxID");
    int maxid = -1;
    for (size_t k=0; k<xdatasets.size(); k++) {
      maxid = Max(maxid, xdatasets[k].setid());
    }
    if (ID > 0) {
      maxid=  Max(maxid, ID);
    }
    return maxid;
  }
  //--------------------------------------------------------------
  //! format all PxdNames
  std::string Dataset::formatNames() const
  {
    dieIfEmpty("formatNames");
    std::string s = "";
    size_t ndts = xdatasets.size();
    for (size_t k=0; k<xdatasets.size(); k++) {
      s += xdatasets[k].pxdname().format();
      if (k < ndts-1) { // not last
        s += " and ";
      }
    }
    return s;
  }
  //--------------------------------------------------------------
  std::string Dataset::formatPrint() const //!< format
  {
    dieIfEmpty("formatPrint");
    std::string s = " * Dataset information *\n";
    size_t ndts = xdatasets.size();
    for (size_t k=0; k<xdatasets.size(); k++) {
      s += xdatasets[k].formatPrint();
      if (k < ndts-1) { // not last
        s += "\n";
      }
    }
    return s;
  }
  //--------------------------------------------------------------
  std::string Dataset::format() const //!< format
  {
    dieIfEmpty("format");
    std::string s;
    size_t ndts = xdatasets.size();
    for (size_t k=0; k<xdatasets.size(); k++) {
      bool first = (k==0);
      s += xdatasets[k].formatPrint(first);
      if (k < ndts-1) { // not last
        s += "\n";
      }
    }
    return s;
  }
  //--------------------------------------------------------------
  //! equality, just tests all pxdnames
  bool operator == (const Dataset& a,const Dataset& b)
  {
    int na = a.Number();
    int nb = b.Number();
    if (na != nb) return false;  // must have same number
    for (int i=0; i<na; i++) {
      if (!(a.pxdname() == b.pxdname())) return false;
    }
    return true; // all same
  }
  //--------------------------------------------------------------
  // Get index for xdataset setid, = -1 if absent
  int Dataset::XdatasetIndex(const int& setid) const
  {
    int idx = -1;
    if (xdatasets.size() <= 0) return idx;
    for (size_t i=0;i<xdatasets.size();++i) {
      if (setid == xdatasets[i].setid()) {
        return i;
      }}
    return idx;
  }
  //--------------------------------------------------------------
  // Get index for xdataset PxdName, = -1 if absent
  int Dataset::XdatasetIndex(const PxdName& pxdname) const
  {
    int idx = -1;
    if (xdatasets.size() <= 0) return idx;
    for (size_t i=0;i<xdatasets.size();++i) {
      if (pxdname == xdatasets[i].pxdname()) {
        return i;
      }}
    return idx;
  }
  //--------------------------------------------------------------
  // Get index for Xname, = -1 if absent
  int Dataset::XdatasetIndex(const std::string& xname) const
  {
    int idx = -1;
    if (xdatasets.size() <= 0) return idx;
    for (size_t i=0;i<xdatasets.size();++i) {
      if (xname == xdatasets[i].pxdname().xname()) {
        return i;
      }}
    return idx;
  }
  //--------------------------------------------------------------
  void Dataset::check() const
  // Sanity check, all xdatasets should have same Dname
  // Dies here if not
  {
    if (xdatasets.size() <= 1) return;
    for (size_t i=1;i<xdatasets.size();++i) { // loop from 2nd
      if (xdatasets[0].pxdname().dname() != xdatasets[i].pxdname().dname()) {
        // Die, die, die!
        std::string names = xdatasets[0].pxdname().dname()+" != "
          +xdatasets[i].pxdname().dname();
        Message::message(Message_fatal
                         ("Dataset::check failed "+names));
      }
    }
  }
  //--------------------------------------------------------------
  void Dataset::dieIfEmpty(const std::string& where) const
  // Die here if nothing in list
  {
    if (xdatasets.size() <= 0) {
      Message::message(Message_fatal
                       ("Dataset:: empty dataset list, function: "+where));
    }
  }
  //--------------------------------------------------------------
  //! return number of cells/wavelengths
  int Dataset::NumberofCells() const
  {
    int n = 0;
    for (size_t i=0; i<xdatasets.size(); i++) {
      n += xdatasets[i].NumberofCells();
    }
    return n;
  }
  //--------------------------------------------------------------
  std::string Dataset::formatAllCells() const //!< format cell & wavelength list if more than one
  //!< format cell & wavelength list if more than one
  {
    UnitCellSet allcells = AllCellSet();
    std::string s;
    if (allcells.Number() <= 1) return s;
    std::string blank(22,' ');
    double devmax = WorstDeviation();
    const double TOL = 0.001;
    if (devmax < TOL) {
      s += "         Files for dataset contain "+clipper::String(allcells.Number(), 3)+
        " near identical cells, "+
        " maximum deviation "+StringUtil::ftos(devmax, 6,4)+"\n";
      return s;
    }
    s += "\n"+std::string(10,' ')+
      "'deviation' is estimate of the difference of the cell (in A)\n";
    s += std::string(12,' ')+
      "from the mean of the others\n";
    s += std::string(10,' ')+
      "At worst this should be less the half the maximum resolution of the data\n\n";
    s += blank +
      "    a       b       c     alpha    beta   gamma  deviation lambda\n";
    if (allcells.Number() <= 1) return s;
    std::vector<Scell> cells = allcells.Cells();  // cells
    // differences from mean of others
    std::vector<double> dv = allcells.Deviations();
    std::vector<float> allwavelengths = AllWavelengths();
    ASSERT (allwavelengths.size() == cells.size());

    for (int i=0;i<allcells.Number();++i) {
      s += blank;
      s += cells[i].formatPrint(false)+
        StringUtil::ftos(dv[i], 7,2)+
        StringUtil::ftos(allwavelengths[i],9,4)+
        "\n";
    }
    s += std::string(7,' ')+"RMS deviation: ";
    dv = allcells.RmsD();
    for (int i=0;i<6;++i) {
      s += StringUtil::ftos(dv[i], 7,2)+" ";
    }
    return s;
  }
  //--------------------------------------------------------------
  //! return worst deviation (A), = 0 if only one
  double Dataset::WorstDeviation() const
  {
    double worst = -100000.;
    for (size_t i=0; i<xdatasets.size(); i++) {
      worst = Max(worst, xdatasets[i].WorstDeviation());
    }
    return worst;
  }
  //--------------------------------------------------------------
  //--------------------------------------------------------------
  bool in_datasets(const int& setid,
                   const std::vector<Dataset>& datasets,
                   int& idataset)
  // Return true if dataset setid is in datasets list
  //  & return dataset index idataset (-1 if not)
  // If setid == 0, assign to first dataset
  {
    ASSERT (datasets.size() > 0);
    if (setid <= 0) {
      idataset = 0;
      return true;
    }
    for (size_t k = 0; k < datasets.size(); ++k) {
      if (datasets[k].IsSetidPresent(setid)) {
        idataset = k;
        return true;
      }
    }
    idataset = -1;
    return false;
  }
  //--------------------------------------------------------------
  bool in_datasets(const PxdName& pxdname,
                   const std::vector<Dataset>& datasets,
                   int& idataset)
  // Return true if dataset pxdname is in datasets list
  //  & return dataset index idataset (-1 if not)
  // If pxdname is blank, assign to first dataset
  {
    ASSERT (datasets.size() > 0);
    if (pxdname.is_blank()) {
      idataset = 0;
      return true;
    }
    for (size_t k = 0; k < datasets.size(); ++k) {
      if (datasets[k].IsPxdPresent(pxdname)) {
        idataset = k;
        return true;
      }
    }
    idataset = -1;
    return false;
  }
}
