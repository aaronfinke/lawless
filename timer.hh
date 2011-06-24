// timer.hh
//
// A timer class

#ifndef TIMER_HEADER
#define TIMER_HEADER

#include <ctime>

class Timer
{
public:
  Timer();      //!< construct and start
  void Start(); //!< start clock
  double Stop();  //!< reset clock & return time in seconds
  //! return time in seconds
  double Dtime() const;

private:
  std::clock_t t0;  // start time
  std::clock_t t1;  // end time
};

#endif
