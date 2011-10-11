// timer.cpp
//
//  A simple timer class

#include "timer.hh"
#include <iostream>
#include <sstream> 

//--------------------------------------------------------------
Timer::Timer()
{
  Start();
}
//--------------------------------------------------------------
void Timer::Start()
{
  et0 = std::time(NULL);
  t0 = std::clock();
}
//--------------------------------------------------------------
double Timer::Stop()
{
  t = Dtime();
  det1 = std::time(NULL) - et0;
  et0 = std::time(NULL);
  t0 = std::clock();
  return t;
}
//--------------------------------------------------------------
double Timer::Dtime() const
{
  return (std::clock()-t0)/double(CLOCKS_PER_SEC);
}
//--------------------------------------------------------------
double Timer::Etime() const
//! return elapsed time
{return det1;}
//--------------------------------------------------------------
std::string Timer::format(const bool& stop)
{
  if (stop) Stop();
  std::ostringstream oss;
  oss << "cpu time: " << t << " secs, elapsed time: " << det1 << " secs";
  return oss.str();
}
//--------------------------------------------------------------

