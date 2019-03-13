#ifndef __GLOBALS_H__
#define __GLOBALS_H__ 

#include <cstdint>
#include <vector>
#include <numeric>
#include <cstdio>

#include "ns3/nstime.h"

#include "util.h"

// forward declarations:
std::vector<uint8_t> calculateR();

#define STRATEGY_INSPECT 0
#define STRATEGY_HNC     1
#define STRATEGY_RLNC    2
#define STRATEGY_RLNC1   3
#define STRATEGY_RLNC2   4
#define STRATEGY_RLNC3   5
#define STRATEGY_RLNC4   6
// ...

extern uint32_t g_strategy; // the algorithm used for generating linear combinations


#define OUTPUT_MODE_SINK_TIMING 1

extern uint32_t g_output_mode;
extern std::vector<uint8_t> r; // global layer layout vector (non cumulative)
extern uint32_t b; // global block size
extern uint32_t k_ahead; // maximum layer deviation
extern double   feedbackRatio;
extern double   g_dropoutFactor;

// below: some derived utility names that are in sync with paper notation

// generation size
uint8_t n();

// number of layers
uint8_t abs_r();

// a dummy dataset
extern std::vector<std::vector<uint8_t> > g_dummy_dataset;

//
// internal utilities:
//
template<typename T> std::vector<T> calculateR() {
  std::vector<T> newR(abs_r());
  for(size_t l=0;l<abs_r(); l++) {
    newR.at(l) = std::accumulate(r.begin(), r.begin()+l+1, 0);
  }
  return newR;
}

// cumulated layer vector
template<typename T = uint8_t> const std::vector<T>& R() {
  static std::vector<T> R = calculateR<T>();
  return R;
}

// lowest-priority decodable layer index
template<typename Ty> int const bestLayer(const std::vector<Ty>& T) {
  // // cumulate tau for Tau
  // std::vector T;
  // for(size_t l=0; l<abs_r(); l++) {
  //   T.at(l) = std::accumulate(r.begin(), r.begin()+l+1, 0);
  // }

  // find lowest-priority decodable layer (ie, greatest index)
  int lowest_priority_decodable = -1;
  for(size_t l=0; l<abs_r(); l++) {
    NS_ASSERT( T.at(l) <= R().at(l) );
    if( T.at(l) == R().at(l) ) {
      lowest_priority_decodable = l;
    }
  }
  
  return lowest_priority_decodable;
}


// for timing output results
extern uint32_t g_distance;
extern std::map<CycleId,ns3::Time> g_injectionTime;

#endif //__GLOBALS_H__
