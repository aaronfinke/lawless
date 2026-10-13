// probe.cpp

#include "probe.hh"

namespace scala {
  Probe::Radiation  Probe::radiation  = Probe::XRAY;
  Probe::Instrument Probe::instrument = Probe::TOF;
  bool Probe::given = false;
  //--------------------------------------------------------------
  void Probe::Set(const Radiation& Rad, const Instrument& Inst)
  {
    radiation = Rad;
    instrument = Inst;
    given = true;
  }
  //--------------------------------------------------------------
  std::string Probe::Name()
  {
    if (radiation == XRAY) return "X-ray";
    if (instrument == QUASILAUE) return "neutron, quasi-Laue (reactor)";
    return "neutron, time-of-flight";
  }
}
