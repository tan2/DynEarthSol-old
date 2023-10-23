#ifndef DYNEARTHSOL3D_UTILS_HPP
#define DYNEARTHSOL3D_UTILS_HPP

#include "parameters.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <cfloat>
#include <math.h>
#include <iomanip>
#if defined(_WIN32)
#include <windows.h>
#else
#include <time.h>
#endif

static void print(std::ostream& os, const double& x)
{
  os << x;
}


static void print(std::ostream& os, const int& x)
{
  os << x;
}


static void print(std::ostream& os, const std::size_t& x)
{
  os << x;
}


template <typename T1, typename T2>
void print(std::ostream& os, const std::pair<T1,T2>& x)
{
    os << x.first << ':' << x.second;
}


template <typename T>
void print(std::ostream& os, const T& A, std::size_t size)
{
  os << "[";
  for (std::size_t i = 0; i != size; ++i) {
    print(os, A[i]);
    if (i+1 != size)
      os << ", ";
  }
  os << "]";
}


template <typename Array>
void print(std::ostream& os, const Array& A)
{
  typename Array::const_iterator i;
  os << "[";
  for (i = A.begin(); i != A.end(); ++i) {
    print(os, *i);
    os << ", ";
  }
  os << "]";
}


static double trace(const double* s)
{
#ifdef THREED
    return s[0] + s[1] + s[2];
#else
    return s[0] + s[1];
#endif
}


static double second_invariant2(const double* t)
{
#ifdef THREED
    double a = (t[0] + t[1] + t[2]) / 3;
    return ( 0.5 * ((t[0]-a)*(t[0]-a) + (t[1]-a)*(t[1]-a) + (t[2]-a)*(t[2]-a))
             + t[3]*t[3] + t[4]*t[4] + t[5]*t[5] );
#else
    return 0.25*(t[0]-t[1])*(t[0]-t[1]) + t[2]*t[2];
#endif
}


static double second_invariant(const double* t)
{
    /* second invariant of the deviatoric part of tensor t
     * defined as: td = deviatoric(t); sqrt( td(i,j) * td(i,j) / 2)
     */
    return std::sqrt(second_invariant2(t));
}


static int findNearestNeighbourIndex( double x_new, const double_vec& x )
{
    /* find nearest neighbour index for interpolation
     * x vector only can be ascending
     */
    double dist = DBL_MAX;
    int idx = -1;
    for (size_t i = 0; i < x.size(); ++i ) {
        double newDist = x_new - x[i];
        if ( newDist >= 0 && newDist <= dist ) {
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}


static int64_t get_nanoseconds() {
    #if defined(_WIN32)
    LARGE_INTEGER frequency, counter;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&counter);
    return (int64_t)((double)counter.QuadPart / frequency.QuadPart * 1e9);
    #else
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (int64_t)ts.tv_sec * 1e9 + ts.tv_nsec;
    #endif
}

static void print_time_ns(const int64_t duration) {
    int hours = duration / (int64_t)3600000000000;
    int minutes = (duration % (int64_t)3600000000000) / (int64_t)60000000000;
    double seconds = (duration % (int64_t)60000000000) / 1e9;
    std::cout << std::setw(3) << std::setfill('0') << hours << ":"
    << std::setw(2) << std::setfill('0') << minutes << ":"
    << std::setw(9) << std::fixed << std::setprecision(6) << std::setfill('0') << seconds;
}

#endif
