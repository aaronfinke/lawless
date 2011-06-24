// timer.cpp
//
//  A simple timer class

#include "timer.hh"
#include <iostream>

//--------------------------------------------------------------
Timer::Timer()
{
  Start();
}
//--------------------------------------------------------------
void Timer::Start()
{
  t0 = std::clock();
}
//--------------------------------------------------------------
double Timer::Stop()
{
  double t = Dtime();
  t0 = std::clock();
  return t;
}
//--------------------------------------------------------------
double Timer::Dtime() const
{
  return (std::clock()-t0)/double(CLOCKS_PER_SEC);
}
