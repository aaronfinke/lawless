//  hkl_controls.hh
//
// Selection controls for hkl_unmerge

#ifndef HKL_CONTROLS_HEADER
#define HKL_CONTROLS_HEADER

#include "controls.hh"
#include "range.hh"
#include "InputAll.hh"


namespace scala
{
  //===================================================================
  class file_select
  //! Various things to control general selection from file

  //!     - dataset selection
  //!     - batch selection
  //!     - resolution limits (on s = 1/d**2)
  //!     - detector coordinate rejection ranges
  //!     - criteria for blank image
  //!
  //! Also store totals rejected for various reasons
  {
  public:
    file_select();
    //! constructor from input
    file_select(phaser_io::InputAll input, const int& NumFileSeries);

    void set_reslimits(const ResoRange& resrange); //!< set resolution limits
    ResoRange reslimits() const {return range_sel;} //!< return resolution limits

    //! test against limits, true if inside
    bool in_reslimits(const Rtype& s) const;

    bool accept_dataset(const std::vector<PxdName>& pxdname) const
    {return true;} //!< always true for now

    //! Return false if batch_number is in rejection lists
    bool accept_batch(const int& batch_number, const int& fileSeries) const;

    BatchSelection BatchExclude() const {return batchexclude;} //!< return list of excluded batches

    void incr_rej_mflag() {nrej_mflag++;} //!< Increment rejection counts
    void incr_rej_reso() {nrej_reso++;} //!< Increment rejection counts
    int Nrej_mflag() const {return nrej_mflag;} //!< return rejection counts
    int Nrej_reso() const {return nrej_reso;} //!< return rejection counts

    float  InputScale() const {return inputscale;}  //!< get scale to be applied on input
    float& InputScale() {return inputscale;}        //!< set scale to be applied on input

    //! Laue Lorentz correction applied on input, see the LORENTZ keyword.
    //! Off unless explicitly asked for: integration programs differ over
    //! whether they have already applied it, and applying it twice is
    //! as wrong as not applying it at all.
    //!   I' = I * sin^2(theta) / lambda^power
    //! power = 0 for a stationary crystal in a white beam (the geometric Laue
    //! factor, as LAUEGEN/LAUENORM and LSM use it), power = 4 for
    //! time-of-flight Laue, where the flight-time bin width adds lambda^4.
    bool  Lorentz() const {return lorentz;}
    void  SetLorentz(const bool& on, const double& power)
      {lorentz = on; lorentzpower = power;}
    double LorentzPower() const {return lorentzpower;}

    //! Test for blank batches: get fraction of maximum resolution to use
    float NullResFrac() const {return nullResolutionfraction;} // get
    //! Test for blank batches: set fraction of maximum resolution to use
    float& NullResFrac() {return nullResolutionfraction;}      // set
    //! Test for blank batches: get threshold on proportion of negative reflections
    float NullNegativeReject() const {return nullNegativeReject;}      // get
    //! Test for blank batches: set threshold on proportion of negative reflections
    float& NullNegativeReject() {return nullNegativeReject;}           // set

  private:
    void init();
    int nrej_mflag, nrej_reso;
    ResoRange range_sel;
    BatchSelection batchexclude;  // List/ranges of batches to exclude
    BatchSelection batchinclude;  // List/ranges of batches to include (from RUN)
    float inputscale;             // scale factor to apply on input
    bool lorentz;        // LORENTZ keyword given
    double lorentzpower; // power of lambda in the correction
    // fraction of maximum resolution to use in test for blank batches
    float nullResolutionfraction;
    // Threshold on proportion of negative reflections
    float nullNegativeReject;

  }; // file_select
}
#endif
