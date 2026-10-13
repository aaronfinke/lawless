// probe.hh
//
// Radiation type and instrument class, from the PROBE keyword
//
//   PROBE  NEUTRON | XRAY            (default XRAY)
//   PROBE  NEUTRON TOF | QUASILAUE   (default TOF)
//
// PROBE changes defaults, guards and reporting only.  It never applies a
// correction which cannot be seen in the log or overridden by another
// keyword, and everything it switches stays independently settable, so an
// X-ray Laue user can borrow a neutron term and a neutron user can turn one
// off.  Everything it changes is listed in the log by ReportProbeDefaults().
//
// This is a static class, like SelectI in hkl_unmerge.hh, so that reporting
// code deep in the program can ask what the radiation is without threading a
// control object through every signature.

#ifndef PROBE_HEADER
#define PROBE_HEADER

#include <string>

namespace scala {
  //--------------------------------------------------------------
  class Probe
  {
  public:
    enum Radiation {XRAY, NEUTRON};
    //! Time-of-flight separates diffraction orders in time, so there are no
    //! harmonics to deconvolute.  A quasi-Laue (reactor) instrument does not.
    enum Instrument {TOF, QUASILAUE};

    Probe(){}

    //! store the setting from the PROBE keyword
    static void Set(const Radiation& Rad, const Instrument& Inst);

    static bool IsNeutron() {return radiation == NEUTRON;}
    static bool IsXray()    {return radiation == XRAY;}
    //! true for a neutron instrument which does not separate orders
    static bool IsQuasiLaue() {return radiation == NEUTRON && instrument == QUASILAUE;}
    //! true if PROBE was given explicitly
    static bool Given() {return given;}
    //! "X-ray" or "neutron, time-of-flight", for the log
    static std::string Name();

  private:
    static Radiation radiation;
    static Instrument instrument;
    static bool given;
  };
}
#endif
