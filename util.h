#ifndef __UTIL_H__
#define __UTIL_H__

#include <vector>
#include <cstdint>
#include <tuple>
#include <string>
#include <functional>

#include "ns3/network-module.h"
#include "ns3/random-variable-stream.h"

#include "network_coding.hxx"

typedef std::tuple<uint16_t, uint16_t, uint16_t> CycleId;

uint64_t MakeNumericCycleId(CycleId id);
std::string ShowCycleId(CycleId id);
std::string ShowCycleId(CycleId id, uint16_t offset);

void callFkt(std::function<void()> f);

ns3::Ipv4Address IdToIp(int id);

class Ns3GfRng : public ff::gf_rng
{
 public:
  Ns3GfRng();
  virtual ff::gf_elem roll();
  virtual ~Ns3GfRng() {};
 private:
  ns3::Ptr<ns3::UniformRandomVariable> rng;
};

template<typename T1,typename T2> std::vector<T1> rebase_vec(std::vector<T2> old) {
  std::vector<T1> res;
  res.insert(res.begin(), old.begin(), old.end());
  return res;
}

#endif //__UTIL_H__
