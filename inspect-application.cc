/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */

#include <map>
#include <set>
#include <utility>
#include <sstream>
#include <numeric>
#include <memory>
#include <iterator>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/applications-module.h"
#include "ns3/socket.h"
#include "ns3/udp-socket-factory.h"
#include "ns3/random-variable-stream.h"
#include "ns3/ipv4.h"

#include "inspect-application.h"
#include "lambda-callbacks.h"
#include "globals.h"
#include "packet-format.h"
#include "network_coding.hxx"

using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("InspectApplication");

NS_OBJECT_ENSURE_REGISTERED (InspectApplication);

struct InspectApplication::PrivT {
  Ptr<Socket> socket;
  Ptr<ExponentialRandomVariable> exp;
  Ptr<UniformRandomVariable> unif;
  
  // state for receiving, sending and injecting
  std::map<CycleId, std::shared_ptr<ff::nc_generator> > generatorStateSet;
  std::map<CycleId, std::shared_ptr<ff::naive_nc_absorber> > decoderStateSet;

  // Idea for improvement:
  // neighbors are at first "sporadic". When a second message is received within 8 feedback intervals, that neighbor is promoted to a "full" neighbor.
  // a regular timeout function demotes neighbors that did not receive feedback for at 4 intervals from "full" to "sporadic" and deletes neighbors after "8" intervals.
  // this gives neighbors with approximately >25% PDR and sporadic neighbors with approximately >12.5% PDR
  //std::map<uint16_t, ns3::Time> sporadicNeighbors;
  std::map<uint16_t, ns3::Time> fullNeighbors;
  std::shared_ptr<Ns3GfRng> gfRng;

  // feedback is sent only for incomplete generations, or (once) those for which a linear combination was received recently
  // this data structure tracks which generations are "active" with respect to feedback, i.e., for which generations such a
  // feedback message should be produced
  std::map<CycleId, bool> activeCycles;

  // state for injecting
  uint16_t nextCycleInjected;

  // state for sending
  
  // starts with 0, increases by one for each sent linear combination. 
  // used to implement approximate round robin send strategy for multiple generations.
  uint16_t nextGenStateOffset;
  uint16_t nextAbsStateOffset;

  // an optimization to stop nodes from trying to send things when there is nothing to send
  // will be set to true once a node determined that there is just nothing to do (all generations complete)
  // will be reset to false once a node (1) injects or (2) receives feedback.
  bool nothingToDo;
};

TypeId
InspectApplication::GetTypeId (void) {
  static TypeId tid = TypeId ("ns3::InspectApplication")
    .SetParent<Application> ()
    .SetGroupName("Applications")
    .AddConstructor<InspectApplication> ()
    .AddAttribute ("DataRate", "The data rate when sending",
                   DataRateValue (DataRate ("1Mbps")),
                   MakeDataRateAccessor (&InspectApplication::m_dataRate),
                   MakeDataRateChecker ())
    .AddAttribute ("NumSensors", "Number of sensors",
                   UintegerValue (4),
                   MakeUintegerAccessor (&InspectApplication::m_numSensors),
                   MakeUintegerChecker<uint32_t> (0,16))
    .AddAttribute ("CycleDuration", "Duration of the cycle",
                   TimeValue (Seconds(25)),
                   MakeTimeAccessor (&InspectApplication::m_cycleDuration),
                   MakeTimeChecker ())
    .AddAttribute ("InjectionStart", "Start delay of the first cycle",
                   TimeValue (Seconds(1)),
                   MakeTimeAccessor (&InspectApplication::m_injectionStart),
                   MakeTimeChecker ())
    .AddAttribute ("MaxCycles", "Number of cycles to send, 0 means no limit",
                   UintegerValue (0),
                   MakeUintegerAccessor (&InspectApplication::m_maxCycles),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("CycleOffset", "An offset that is added to the cycle counter",
                   UintegerValue (0),
                   MakeUintegerAccessor (&InspectApplication::m_cycleOffset),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("BlockSize", "Size of all blocks in bytes",
                   UintegerValue (1024),
                   MakeUintegerAccessor (&InspectApplication::m_blockSize),
                   MakeUintegerChecker<uint32_t> (1,4096))
    ;
  return tid;
}

InspectApplication::InspectApplication ()
  :m_priv(new PrivT())
{
  m_priv->exp = CreateObject<ExponentialRandomVariable> ();
  m_priv->unif = CreateObject<UniformRandomVariable> ();
  m_priv->nextCycleInjected = 0;
  m_priv->gfRng = std::make_shared<Ns3GfRng>();
  m_priv->nextGenStateOffset = 0;
  m_priv->nextAbsStateOffset = 0;
  m_priv->nothingToDo = false;
}

void InspectApplication::StartApplication() {
  NS_LOG_DEBUG(Simulator::Now().GetSeconds() << ": " << GetOwnMachineId() << " started");

  m_priv->socket = Socket::CreateSocket(GetNode(), UdpSocketFactory::GetTypeId ());
  m_priv->socket->SetAllowBroadcast(true);
  m_priv->socket->Bind(InetSocketAddress(Ipv4Address::GetAny(), 9090));
  
  m_priv->socket->SetRecvCallback(MakeCallback(&InspectApplication::Rx, this));

  if(!IsSink()) {
    PeriodicTx();
  }

  PeriodicFeedback();

  PeriodicMaintenance();

  if(!IsSink()) {
    Simulator::Schedule(m_injectionStart, &InspectApplication::Inject, this);
  }
}

void InspectApplication::StopApplication() {
  
}

void InspectApplication::PeriodicTx () {
  const size_t genNumStates = m_priv->generatorStateSet.size();
  const size_t absNumStates = m_priv->decoderStateSet.size();

  // check whether incomplete generation queue is empty and wait 10ms if it is is empty
  if( genNumStates + absNumStates == 0 || m_priv->nothingToDo == true ) {
    Simulator::Schedule(MilliSeconds(10), &InspectApplication::PeriodicTx, this);
    return;
  }

  CycleId cycleId; // will be set later
  InspectHeader            header;
  header.messageType = INSPECT_LINEAR_COMBINATION;
  // ... rest of fields will be set once cycle ID is known
  InspectLinearCombination lcHeader;
  // ... will be set once cycle ID is known

  // figure out generation to use for sending
  // if available, use own generation, else, use absorbers in a round-robin fashion
  std::shared_ptr<ff::nc_general> generationState;
  
  if(genNumStates) {
    size_t initialOffset = (m_priv->nextGenStateOffset) % genNumStates;
    size_t roundRobinIdx = initialOffset;
    do {
      // case (1) generate LC from own generation
      auto it = m_priv->generatorStateSet.begin();
      std::advance(it, roundRobinIdx);
      cycleId         = (*it).first;
      generationState = (*it).second;
      NS_ASSERT(generationState);
      if(generationState->completed) {
        // skip fully sent generations
        // idea: remove state after a timeout to improve performance
        NS_LOG_DEBUG(Simulator::Now().GetSeconds() << ": " << GetOwnMachineId() << " SKIPS generation " << ShowCycleId(cycleId)
                     << " total gen num states " << genNumStates);
        // NS_LOG_DEBUG("... rrIdx: " << roundRobinIdx);
        // NS_LOG_DEBUG("... iiOff: " << initialOffset);
      }
      roundRobinIdx = (++m_priv->nextGenStateOffset) % genNumStates;
    } while((generationState->completed) && roundRobinIdx != initialOffset);
  }

  if(absNumStates && (!generationState || generationState->completed)) {
    // try absorbers instead if no decoder available:
    size_t initialOffset = (m_priv->nextAbsStateOffset) % absNumStates;
    size_t roundRobinIdx = initialOffset;
    do {
      // case (2) generate LC from other node's generation
      auto it = m_priv->decoderStateSet.begin();
      std::advance(it, roundRobinIdx);
      cycleId         = (*it).first;
      generationState = (*it).second;
      NS_ASSERT(generationState);
      if(generationState->completed || 0 == generationState->get_current_degree()) {
        // skip fully sent generations
        // idea for improvement: remove state after a timeout to improve simulation performance
        NS_LOG_DEBUG(Simulator::Now().GetSeconds() << ": " << GetOwnMachineId() << " SKIPS generation " << ShowCycleId(cycleId)
                     << " total abs num states " << absNumStates);
        // NS_LOG_DEBUG("... rrIdx: " << roundRobinIdx);
        // NS_LOG_DEBUG("... iiOff: " << initialOffset);
      }
      roundRobinIdx = (++m_priv->nextAbsStateOffset) % absNumStates;
    } while((generationState->completed || 0 == generationState->get_current_degree()) && roundRobinIdx != initialOffset);
  }

  if(generationState->completed || 0 == generationState->get_current_degree()) {
    // neither generator nor decoder state useful, nothing to do here. continue
    m_priv->nothingToDo = true;
    Simulator::Schedule(MilliSeconds(10), &InspectApplication::PeriodicTx, this);
    return;
  }

  // execute in any case for assumed feedback updates
  size_t layerChoice = generationState->classChoice(k_ahead);
  if (g_strategy == STRATEGY_INSPECT) {
    ; // use layer choice from above
  } else if (g_strategy >= STRATEGY_RLNC) {
    ; // this is addressed in simulation.cc by modifying r to reflect the strategy
  } else if (g_strategy == STRATEGY_HNC) {
    layerChoice = m_priv->unif->GetInteger(0,abs_r() - 1);
  } else {
    NS_ASSERT(false && "invalid strategy choice");
  }
  NS_LOG_DEBUG(Simulator::Now().GetSeconds() << ": " << GetOwnMachineId() << " chooses generation " << ShowCycleId(cycleId) << " layer " << layerChoice);
    
  ff::linear_combination lc;
  if(std::dynamic_pointer_cast<ff::nc_generator>(generationState)) {
    lc = std::dynamic_pointer_cast<ff::nc_generator>(generationState)\
      ->new_zerotail_linear_combination(R<size_t>().at(layerChoice));
    NS_ASSERT(lc.combination.size() == n());
  } else {
    lc = std::dynamic_pointer_cast<ff::naive_nc_absorber>(generationState)\
      ->new_zerotail_linear_combination(layerChoice, n(), b);
    NS_ASSERT(lc.combination.size() == n());
    // if( ! (lc.combination.size() == n()) ) {
    //   printf("relevant node: %u\n", GetOwnMachineId());
    //   printf("relevant combination: %s\n", ShowCycleId(cycleId).c_str());
    //   printf("choice: %zu, size: %zu, n(): %u\n", layerChoice, lc.combination.size(), n());
    //   printf("degs: ");
    //   for ( auto d : std::dynamic_pointer_cast<ff::naive_nc_absorber>(generationState)->get_current_degrees() ) {
    //     printf("%u ", d);
    //   }
    //   printf("\n");
    //   NS_ASSERT(false);
    // }
  }

    
  lcHeader = InspectLinearCombination::FromFfLc (lc);


  header.machineId = std::get<0>(cycleId);
  header.sensorId = std::get<1>(cycleId);
  header.cycleId = std::get<2>(cycleId);

  Ptr<Packet> pkg = Create<Packet>();
  pkg->AddHeader (lcHeader);
  pkg->AddHeader (header);

  // actually transmit packet
  m_priv->socket->SendTo(pkg,
                         0,
                         InetSocketAddress(Ipv4Address::GetBroadcast(), 9090)
                         );
  NS_LOG_DEBUG(Simulator::Now().GetSeconds() << ": " << GetOwnMachineId() << " sent packet ");


  // reschedule sending event after exponentially distributed send time
  Time txTime = m_dataRate.CalculateBytesTxTime (pkg->GetSize ());
  Time txDelay = MicroSeconds (m_priv->exp->GetValue (txTime.GetMicroSeconds () ,0));
  Simulator::Schedule(txDelay, &InspectApplication::PeriodicTx, this);
}

void InspectApplication::PeriodicMaintenance () {
  Time dataTxTime = m_dataRate.CalculateBytesTxTime (20 + InspectHeader().GetSerializedSize() + InspectLinearCombination().GetSerializedSize());
  Time feedbackPeriod    = dataTxTime * feedbackRatio;
  Time maintenancePeriod = feedbackPeriod * g_dropoutFactor;
  Time now = Simulator::Now();

  std::set<uint16_t> demoted;
  for( auto full : m_priv->fullNeighbors ) {
    // when x+5 is in the past, we have not received any feedback from neighbor recently
    if(full.second + maintenancePeriod < now) {
      demoted.insert(full.first);
      NS_LOG_DEBUG(Simulator::Now().GetSeconds() << ": " << GetOwnMachineId() << " DEMOTED " << full.first);
    }
  }
  for( auto d : demoted ) {
    m_priv->fullNeighbors.erase(d);

    for( auto gen : m_priv->generatorStateSet ) {
      gen.second->delFB( d );
    }
    for( auto gen : m_priv->decoderStateSet ) {
      gen.second->delFB( d );
    }
  }
  
  Time maintenancePeriod1 = maintenancePeriod/2;
  Time maintenancePeriod2 = MilliSeconds (m_priv->exp->GetValue ((maintenancePeriod/2).GetMilliSeconds () ,0));
  Simulator::Schedule(maintenancePeriod1 + maintenancePeriod2, &InspectApplication::PeriodicMaintenance, this);
}

void InspectApplication::Rx (Ptr<Socket> sock) {
  uint32_t ret;
  Address fromAddress;
  Ptr<Packet> pkg = sock->RecvFrom (fromAddress);
  InetSocketAddress fromIpv4 = InetSocketAddress::ConvertFrom(fromAddress);
  uint16_t fromId = fromIpv4.GetIpv4().Get () % 256;

  if(fromId == GetOwnMachineId()) {
    // our own packet!
    return;
  }
  
  NS_LOG_DEBUG(Simulator::Now().GetSeconds() << ": " << GetOwnMachineId() << " received a packet");
  InspectHeader header;
  ret = pkg->RemoveHeader(header);
  NS_ASSERT(ret && "could not read inspect header");
  // check whether we got feedback or a new linear combination
  if( header.messageType == INSPECT_LINEAR_COMBINATION ) {
    InspectLinearCombination lc;
    ret = pkg->RemoveHeader(lc);
    // ignore if we are the generator ourself
    if( m_priv->generatorStateSet.find(header.GetCycleId()) != m_priv->generatorStateSet.end() ) {
      return;
    }
    // create new absorber state if it does not exist already
    if( m_priv->decoderStateSet.find(header.GetCycleId()) == m_priv->decoderStateSet.end() ) {
      auto absorber = std::make_shared<ff::naive_nc_absorber>(m_priv->gfRng.get(), R<size_t>());
      m_priv->decoderStateSet[header.GetCycleId()] = absorber;
    }
    auto state = m_priv->decoderStateSet.at(header.GetCycleId());

    auto old_degs = state->get_current_degrees();
    int old_layer = bestLayer(old_degs);

    m_priv->activeCycles[header.GetCycleId()] = true;
    
    if(old_degs.at(abs_r()-1) == n()) {
      // already fully decodable
      return;
    }
    state->saveLinComb( lc.GetFfLc() );
    NS_LOG_DEBUG(Simulator::Now().GetSeconds() << ": " << GetOwnMachineId() << " saved linear combination of gen: " << ShowCycleId(header.GetCycleId()));

    auto new_degs = state->get_current_degrees();
    int new_layer = bestLayer(new_degs);
    bool innovative = false;
    for(int i=old_layer; i<new_layer; i++) {
      NS_ASSERT(i >= -1);
      NS_ASSERT(g_injectionTime.find(header.GetCycleId()) != g_injectionTime.end());
      printf("%u, %u, %u, %2.2f, %u, %d, %s, %2.3f\n",
             g_strategy,
             g_distance,
             (unsigned int)feedbackRatio,
             Simulator::Now().GetSeconds(),
             GetOwnMachineId(),
             i+1,
             ShowCycleId(header.GetCycleId(), m_cycleOffset).c_str(),
             (Simulator::Now() - g_injectionTime.at(header.GetCycleId())).GetSeconds());
      innovative = true;
    }
    
    auto max = state->get_current_degree();
    if(max == n() && innovative) {
      NS_LOG_DEBUG(Simulator::Now().GetSeconds() << ": " << GetOwnMachineId() << " COMPLETED linear combination of gen: " << ShowCycleId(header.GetCycleId()));
    }
  } else if (header.messageType == INSPECT_FEEDBACK) {
    // register neighbor first
    m_priv->fullNeighbors[fromId] = Simulator::Now();
    
    // inspect feedback header
    InspectFeedback fb;
    ret = pkg->RemoveHeader(fb);
    NS_ASSERT(ret && "could not parse feedback message");

    // we need either generator state or decoder (=absorber) state to make sense of this feedback message
    std::shared_ptr<ff::nc_general> state;
    if(m_priv->generatorStateSet.find(header.GetCycleId()) != m_priv->generatorStateSet.end()) {
      // we have generator state for this feedback message's cycle
      state = m_priv->generatorStateSet[header.GetCycleId()];
    }  else if(m_priv->decoderStateSet.find(header.GetCycleId()) != m_priv->decoderStateSet.end()) {
      // we have decoder (=absorber) state for this feedback message's cycle
      state = m_priv->decoderStateSet[header.GetCycleId()];
    } else {
      // create new (empty) decoder state, since we just learned of a new generation
      auto absorber = std::make_shared<ff::naive_nc_absorber>(m_priv->gfRng.get(), R<size_t>());
      state = (m_priv->decoderStateSet[header.GetCycleId()] = absorber);
    }
    if(state) {
      // found either generator state or decoder state
      for(auto i :rebase_vec<uint16_t>(fb.feedbackVector)) {
        NS_LOG_DEBUG("...feedback: " << (int)i);
      }
      state->saveFB(fromId, rebase_vec<uint16_t>(fb.feedbackVector));
      state->checkComplete();
      if(!state->completed) {
        NS_LOG_DEBUG("...no longer nothing to do");
        // this feedback messages caused us to become active again!
        m_priv->nothingToDo = false;
      }
    }
  } else {
    NS_ASSERT(false && "invalid message type");
  }
}
void InspectApplication::PeriodicFeedback () {
  // send feedback for all generations
  for(auto p : m_priv->decoderStateSet) {
    auto cycleId = p.first;
    auto generationState = p.second;

    if (generationState->get_current_degree() == n()) {
      // this generation is fully received, only send if necessary
      if(!m_priv->activeCycles[cycleId]) {
        m_priv->activeCycles.erase(cycleId);
        continue;
      } else {
        // note: the if statement actually creates entry (with default entry: false), so the following erase is needed in this branch here as well!
        m_priv->activeCycles.erase(cycleId);
      }
    }
    
    InspectHeader    header;
    InspectFeedback  fbHeader;

    header.messageType = INSPECT_FEEDBACK;
    header.machineId = std::get<0>(cycleId);
    header.sensorId  = std::get<1>(cycleId);
    header.cycleId   = std::get<2>(cycleId);

    auto degs = generationState->get_current_degrees();
    NS_LOG_DEBUG(Simulator::Now().GetSeconds() << ": " << GetOwnMachineId() << " generating feedback for: " << ShowCycleId(cycleId));
    for (auto d : degs) {
      NS_LOG_DEBUG("..." << (int)d);
    }
    fbHeader.feedbackVector.insert(fbHeader.feedbackVector.begin(), degs.begin(), degs.end());
    NS_ASSERT(fbHeader.feedbackVector.size() == abs_r());

    Ptr<Packet> pkg = Create<Packet>();
    pkg->AddHeader(fbHeader);
    pkg->AddHeader(header);

    m_priv->socket->SendTo(pkg,
                           0,
                           InetSocketAddress(Ipv4Address::GetBroadcast(), 9090)
                           );
  }
  
  // feedback send period has 50% constant and 50% exponential parts
  Time dataTxTime = m_dataRate.CalculateBytesTxTime (20 + InspectHeader().GetSerializedSize() + InspectLinearCombination().GetSerializedSize());
  Time feedbackPeriod = dataTxTime * feedbackRatio;
  Time feedbackTxDelay1 = feedbackPeriod/2;
  Time feedbackTxDelay2 = MicroSeconds (m_priv->exp->GetValue ((feedbackPeriod/2).GetMicroSeconds () ,0));
  Simulator::Schedule(feedbackTxDelay1 + feedbackTxDelay2, &InspectApplication::PeriodicFeedback, this);
}

void InspectApplication::Inject ()
{
  for(uint16_t sensor=1; sensor <= m_numSensors; sensor++) {
    CycleId cycleId(GetOwnMachineId(),sensor,m_priv->nextCycleInjected);

    g_injectionTime[cycleId] = Simulator::Now();
    // PDB printf("%2.2f: ... %u injected %s\n", Simulator::Now().GetSeconds(), GetOwnMachineId(), ShowCycleId(cycleId).c_str());

    
    NS_ASSERT( m_priv->generatorStateSet.find(cycleId) == m_priv->generatorStateSet.end() );
    NS_ASSERT( m_priv->decoderStateSet.find(cycleId) == m_priv->decoderStateSet.end() );

    auto g = std::make_shared<ff::nc_generator>(m_priv->gfRng.get(), &g_dummy_dataset, R<size_t>());
    m_priv->generatorStateSet[cycleId] = g;

    NS_LOG_DEBUG(Simulator::Now().GetSeconds() << ": " << GetOwnMachineId() << " injected generation " << ShowCycleId(cycleId));

    m_priv->nothingToDo = false;
  }
  
  ++m_priv->nextCycleInjected;
  if(m_priv->nextCycleInjected < m_maxCycles || m_maxCycles == 0) {
    Simulator::Schedule(1 * m_cycleDuration, &InspectApplication::Inject, this);
  }
}
bool InspectApplication::IsSink () {
  return GetOwnMachineId() == 1;
}
// determines from last IP octet
uint16_t InspectApplication::GetOwnMachineId () const
{
  Ptr<Ipv4> ipv4 = GetNode()->GetObject<Ipv4> ();
  Ipv4InterfaceAddress iaddr = ipv4->GetAddress (1,0); 
  Ipv4Address ipAddr = iaddr.GetLocal ();
  return ipAddr.Get() % 256;
}
