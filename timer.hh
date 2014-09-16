// timer.hh
//
// A timer class

#ifndef TIMER_HEADER
#define TIMER_HEADER

#include <ctime>
#include <string>

class Timer
{
public:
  Timer();      //!< construct and start
  void Start(); //!< start clock
  double Stop();  //!< reset clock & return time in seconds
  //! return time in seconds
  double Dtime() const;
  //! return elapsed time
  double Etime() const;
  //! return formatted version, stop timer if stop true
  std::string format(const bool& stop=true);


private:
  std::clock_t t0;  // start time
  std::clock_t t1;  // end time
  std::time_t et0;  // start clock time
  double t;     // cpu time
  double det1;  // elapsed clock time
};

#endif
