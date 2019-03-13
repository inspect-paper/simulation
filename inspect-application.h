/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */

#include "ns3/core-module.h"
#include "ns3/application.h"
#include "ns3/data-rate.h"
#include "ns3/socket.h"

#include <memory>
#include <tuple>


namespace ns3 {

class InspectApplication : public Application
{
 public:
  static uint64_t tx;
  static TypeId GetTypeId (void);
  
  InspectApplication ();
  virtual ~InspectApplication(){}
  
 protected:
  virtual void StartApplication (void);
  virtual void StopApplication (void);

  virtual void PeriodicTx ();
  virtual void Rx (Ptr<Socket>);
  virtual void PeriodicFeedback ();
  virtual void PeriodicMaintenance ();

  virtual void Inject ();
  virtual bool IsSink ();
  uint16_t GetOwnMachineId () const; // determined from IP ending

  // attributes:
  DataRate m_dataRate;       // == maximum sending rate
  uint32_t m_numSensors;     // == number of parallel generations
  Time     m_cycleDuration;  // == 3x cycle duration equals wait time between sending (1x cycle duration + 2x cooldown)
  Time     m_injectionStart; // == time when the first injection starts
  uint32_t m_maxCycles;      // == maximum number of generations to transmit before stopping
  uint32_t m_cycleOffset;    // == only used for printing, useful to merge multiple simulations
  uint32_t m_blockSize;      // == size of blocks in generation

  // derived values:
  uint32_t GetMaxNumCoefficients() const;

 private:
  // implementation stuff:
  struct PrivT;
  std::shared_ptr<PrivT> m_priv;
};
  
}
