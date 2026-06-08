// anomalousstatus.hh

#ifndef ANOMALOUSSTATUS_HEADER
#define ANOMALOUSSTATUS_HEADER

namespace scala {
  class AnomalousStatus {
  public:
    // Status of anomalous data:
    //   ON    specified on input as ON
    //   OFF   specified on input as OFF
    //   else  unspecified
    //   FOUND  significant anomalous present
    //   ABSENT significant anomalous below thresholds
    //   NO_DATA  no information
    enum anomalousStatus {ANOMALOUS_ON_FOUND, ANOMALOUS_ON_ABSENT,
			  ANOMALOUS_OFF_FOUND, ANOMALOUS_OFF_ABSENT,
			  ANOMALOUS_FOUND, ANOMALOUS_ABSENT,
			  NO_ANOMALOUS_DATA};


    static bool isFound(const anomalousStatus& anomalousstatus) {
      if (anomalousstatus == ANOMALOUS_ON_FOUND ||
	  anomalousstatus == ANOMALOUS_OFF_FOUND ||
	  anomalousstatus == ANOMALOUS_FOUND) {
	return true;}
      return false;
    }    

  };
}

#endif
