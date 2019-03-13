/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/applications-module.h"
#include "ns3/internet-module.h"
#include "ns3/mobility-module.h"
#include "ns3/wifi-module.h"

#include "lambda-callbacks.h"
#include "globals.h"
#include "inspect-application.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <string>

using namespace ns3;


NS_LOG_COMPONENT_DEFINE ("inspect");

std::vector<uint8_t> split1(const std::string& str)
{
  std::vector<uint8_t> cont;
  std::istringstream ss(str);
  std::string token;
  while (std::getline(ss, token, ',')) {
    cont.push_back(std::stol(token));
  }
  return cont;
}

void PrintPositions(NodeContainer*);

int main (int argc, char *argv[])
{
  // LogComponentEnable("InspectApplication", LOG_LEVEL_DEBUG);

  std::string phyMode ("ErpOfdmRate54Mbps");
  // disable fragmentation
  Config::SetDefault ("ns3::WifiRemoteStationManager::FragmentationThreshold", StringValue ("2200"));
  Config::SetDefault ("ns3::WifiRemoteStationManager::RtsCtsThreshold", StringValue ("2200"));
  // Fix non-unicast data rate to be the same as that of unicast
  Config::SetDefault ("ns3::WifiRemoteStationManager::NonUnicastMode",
    StringValue (phyMode));
  
  uint32_t rngRun = 0;
  uint32_t numSensors = 4;
  uint32_t cycleDuration = 25;
  uint32_t maxCycles = 8;
  double pathLossExponent = 3.0; // default exponent
  double referenceLoss = 40.1;    // reference loss for 2.4 GHz at 1m, by friis (fixed slope)
  double referenceDistance = 1.0;
  double nakagamiM0 = 1.0; // rayleigh
  Time warmup = Seconds(100.0);
  Time duration = Seconds(600.0);
  Time cooldown = Seconds(120.0);
  uint32_t cycleOffset = 0;
  uint32_t gridWidth = 3;
  uint32_t gridLength = 5;
  bool pcap = false;
  bool killOneNode = false;
  bool onlyOneSource = false;
  uint32_t onlyNSources = 0;
  bool disc = false;
  double discRadius;
  DataRate rate = DataRate("1Mbps");
  bool ppos = false;

  // global system parameters
  b = 1024;
  std::string r_string = "2,2,3,3";

  CommandLine cmd;
  cmd.AddValue("rngRun", "run-# of the PRNG", rngRun);
  cmd.AddValue("maxCycles", "maximum number of cycles", maxCycles);
  cmd.AddValue("distance", "h/v distance between two nodes", g_distance);
  cmd.AddValue("cycleDuration", "duration of each cycle", cycleDuration);
  cmd.AddValue("ple", "sets the path loss exponent", pathLossExponent);
  cmd.AddValue("rl", "sets the reference loss", referenceLoss);
  cmd.AddValue("rd", "sets the reference distance", referenceDistance);
  cmd.AddValue("nm", "sets the nakagami m0,1,2 value", nakagamiM0);
  cmd.AddValue("cycleOffset", "an offset to the cycle counter", cycleOffset);
  cmd.AddValue("gridWidth", "width of the grid", gridWidth);
  cmd.AddValue("gridLength", "length of the grid", gridLength);
  cmd.AddValue("outputMode", "stdout mode (1:timing, ...)", g_output_mode);
  cmd.AddValue("pcap", "enable pcap generation", pcap);
  cmd.AddValue("r", "layer configuration", r_string);
  cmd.AddValue("feedbackRatio", "feedbackRatio", feedbackRatio);
  cmd.AddValue("killOneNode", "remove one node to make grid topology deliberately irregular", killOneNode);
  cmd.AddValue("onlyOneSource", "have only one (the highest-index) node send data", onlyOneSource);
  cmd.AddValue("onlyNSources", "have only one (the highest-index) node send data", onlyNSources);
  cmd.AddValue("numSensors", "number of sensors", numSensors);
  cmd.AddValue("strategy", "strategy to use (0:inspect, 1:hnc, 2:rlnc, 3+:baseline rlnc for (n-2)-th layer)", g_strategy);
  cmd.AddValue("rate", "data rate for sending", rate);
  cmd.AddValue("disc", "use disc instead of grid", disc);
  cmd.AddValue("radiusPerNode", "used to calculate total disc area for nodes", g_distance);
  cmd.AddValue("dropoutFactor", "determines the maintenance period as a multiple of the feedback ratio", g_dropoutFactor);
  cmd.AddValue("ppos", "print positions, then exit", ppos);

  cmd.Parse (argc, argv);

  // option derived values

  NS_ASSERT(!(onlyNSources && onlyOneSource) && "You cannot specifiy both options.");

  uint32_t numNodes = gridWidth * gridLength;
  if(killOneNode) {
    --numNodes;
  }
 
  discRadius = std::sqrt (g_distance * g_distance * numNodes);

  // layered mode r configuration
  r = split1(r_string); 
  if(g_strategy == STRATEGY_RLNC) {
    // use one layer only: the lowest-priority layer
    r = { (uint8_t)std::accumulate(r.begin(), r.end(), 0) };
  } else if(g_strategy > STRATEGY_RLNC) {
    // use one layer only: the base-layer
    r = { (uint8_t)std::accumulate(r.begin(), r.begin()+(g_strategy - STRATEGY_RLNC), 0) };
  }
  
  k_ahead = n();
  for (int i=0; i<n(); i++) {
    std::vector<uint8_t> dummyMessage(b);
    // fill with sequential dummy data, could be useful for consistency checks later
    for( int j=0; j<(int)b; j++ ) {
      dummyMessage.at(j) = (i*b+j) & 256;
    }
    g_dummy_dataset.push_back(dummyMessage);
  }

  if(cycleOffset == 0) {
    cycleOffset = rngRun*numNodes*maxCycles;
  }
  
  RngSeedManager::SetRun(rngRun);
  
  WifiHelper wifi;
  YansWifiPhyHelper wifiPhy = YansWifiPhyHelper::Default ();
  YansWifiChannelHelper wifiChannel;
  WifiMacHelper wifiMac;
  
  wifi.SetStandard (WIFI_PHY_STANDARD_80211g);
  
  wifiChannel.SetPropagationDelay ("ns3::ConstantSpeedPropagationDelayModel");
  wifiChannel.AddPropagationLoss ("ns3::LogDistancePropagationLossModel",
                                  "Exponent", DoubleValue (pathLossExponent),
                                  "ReferenceLoss", DoubleValue(referenceLoss),
                                  "ReferenceDistance", DoubleValue(referenceDistance));
  wifiChannel.AddPropagationLoss ("ns3::NakagamiPropagationLossModel",
                                  "m0", DoubleValue (nakagamiM0),
                                  "m1", DoubleValue (nakagamiM0),
                                  "m2", DoubleValue (nakagamiM0));
  wifiPhy.SetChannel (wifiChannel.Create ());
  wifiPhy.SetPcapDataLinkType (YansWifiPhyHelper::DLT_IEEE802_11_RADIO);

  // Add a non-QoS upper mac, and disable rate control
  wifi.SetRemoteStationManager ("ns3::ConstantRateWifiManager",
                                "DataMode",StringValue (phyMode),
                                "ControlMode",StringValue (phyMode));

  Ssid ssid = Ssid ("preview");
  wifiMac.SetType ("ns3::AdhocWifiMac",
                   "Ssid", SsidValue (ssid));
  
  NodeContainer nodes;
  nodes.Create(numNodes);
  MobilityHelper mobility;
  if (disc) {
    mobility.SetPositionAllocator(
                                  "ns3::UniformDiscPositionAllocator",
                                  "rho", DoubleValue (discRadius)
                                  );
  } else {
    mobility.SetPositionAllocator("ns3::GridPositionAllocator",
                                  "GridWidth", UintegerValue (gridWidth),
                                  "DeltaX", DoubleValue (g_distance),
                                  "DeltaY", DoubleValue (g_distance));
  }
  mobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  mobility.Install(nodes);

  NetDeviceContainer devices = wifi.Install (wifiPhy, wifiMac, nodes);

  InternetStackHelper stack;
  stack.Install (nodes);
  Ipv4AddressHelper address;
  address.SetBase ("10.1.1.0", "255.255.255.0");
  Ipv4InterfaceContainer interfaces = address.Assign(devices);

  if(g_output_mode > 0) {
    // {
    //   Ptr<Application> sinkApp = CreateObject<InspectApplication>();
    //   sinkApp->SetStartTime(Seconds(10.0));
    //   nodes.Get(0)->AddApplication(sinkApp);
    // }
    auto rng = CreateObject<UniformRandomVariable> ();
    for(size_t i=0; i<numNodes; i++)
      {
        Ptr<Application> app = CreateObject<InspectApplication>();
        Time delay = Seconds(rng->GetValue (0,cycleDuration));
        app->SetStartTime(warmup);
        app->SetAttribute("InjectionStart", TimeValue(delay));
        app->SetAttribute("DataRate", DataRateValue (rate));
        app->SetAttribute("MaxCycles", UintegerValue (maxCycles));
        app->SetAttribute("CycleDuration", TimeValue (Seconds(cycleDuration)));
        // app->SetAttribute("NumBlocks", UintegerValue(51));
        app->SetAttribute("MaxCycles", UintegerValue(maxCycles));
        if ((!onlyOneSource && !onlyNSources) || (onlyOneSource && i==numNodes-1) || (onlyNSources && i<onlyNSources)) {
          app->SetAttribute("NumSensors", UintegerValue(numSensors));
        } else {
          app->SetAttribute("NumSensors", UintegerValue(0));
        }
        nodes.Get(i)->AddApplication(app);
      }
    // cycleOffset is useful for data aggregation over multiple simulations, only affects sink since it prints output
    Config::Set ("/NodeList/*/ApplicationList/*/$ns3::InspectApplication/CycleOffset",
                 UintegerValue(cycleOffset));
    // sink should not inject cycles or send generations
    // Config::Set ("/NodeList/0/ApplicationList/*/$ns3::InspectApplication/NumBlocks",
    //              UintegerValue(0));
  }

  if(pcap) {
    wifiPhy.EnablePcap ("inspect-wlan", devices, true);
  }

  if(ppos) {
    Simulator::Schedule(Seconds(10.0), &PrintPositions, &nodes);
    Simulator::Stop(Seconds(11.0));
  }
  
  Simulator::Stop(warmup+duration+cooldown);
  Simulator::Run();
  Simulator::Destroy();
  
  return 0;
}

void PrintPositions(NodeContainer* nodes)
{
  uint32_t nodeNum = 0;
  std::cerr << "Positions at time t=" << Simulator::Now().GetMilliSeconds()/1000 << "s:\n";
  for(NodeContainer::Iterator it = nodes->Begin();
      it != nodes->End();
      it++) {
    Vector pos = (*it)->GetObject<MobilityModel>()->GetPosition();
    std::cout << "    Node #" << (++nodeNum) << " has pos " << pos << std::endl;
  }
}
