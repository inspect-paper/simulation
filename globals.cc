#include "globals.h"

uint32_t g_output_mode = OUTPUT_MODE_SINK_TIMING;
std::vector<uint8_t> r; // global layer layout vector (non cumulative)
uint32_t b; // global block size
std::vector<std::vector<uint8_t> > g_dummy_dataset;
uint32_t k_ahead; // maximum layer deviation
double   feedbackRatio = 10.0; // feedback messages per data messages
double   g_dropoutFactor = 5.0; // maintenance period for removing neighbors 

// generation size
uint8_t n() {
  static uint8_t n = std::accumulate(r.begin(), r.end(), 0);
  return n;
}

// number of layers
uint8_t abs_r() {
  static uint8_t l = r.size();
  return l;
}

uint32_t g_strategy = STRATEGY_INSPECT;
uint32_t g_distance = 30.0;
std::map<CycleId,ns3::Time> g_injectionTime;
