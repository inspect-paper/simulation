#include "util.h"

using namespace ns3;
using namespace ff;

uint64_t MakeNumericCycleId(CycleId id) {
  uint64_t res = 0;
  res += std::get<0>(id);
  res  = res << 16;
  res += std::get<1>(id);
  res  = res << 16;
  res += std::get<2>(id);
  res  = res << 16;
  return res;
}

std::string ShowCycleId(CycleId id) {
  std::string s = "XXXX/XXXX/XXXX ";
  snprintf(const_cast<char*>(s.c_str()), s.length(), "%04x/%04x/%04x", std::get<0>(id), std::get<1>(id), std::get<2>(id));
  return s;
}

std::string ShowCycleId(CycleId id, uint16_t offset) {
  std::string s = "XXXX/XXXX/XXXX ";
  snprintf(const_cast<char*>(s.c_str()), s.length(), "%04x/%04x/%04x", std::get<0>(id), std::get<1>(id), std::get<2>(id) + offset);
  return s;
}

void callFkt(std::function<void()> f) {
  f();
};

ns3::Ipv4Address IdToIp(int id) {
  NS_ASSERT(id > 0 && id < 256);
  std::stringstream ss;
  ss << "10.1.1.";
  ss << (int)id;
  return ns3::Ipv4Address(ss.str().c_str());
}

Ns3GfRng::Ns3GfRng()
{
  rng = CreateObject<UniformRandomVariable> ();
}

gf_elem Ns3GfRng::roll()
{
  return rng->GetInteger(0,255);
}
