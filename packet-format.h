#include <ns3/packet.h>
#include <ns3/header.h>

#include <vector>

#include "globals.h"
#include "util.h"
#include "network_coding.hxx" 

#define INSPECT_LINEAR_COMBINATION 0
#define INSPECT_FEEDBACK           1

struct InspectHeader : public ns3::Header {
  static ns3::TypeId GetTypeId (void) {
    static ns3::TypeId tid = ns3::TypeId ("ns3::InspectHeader")
      .SetParent<ns3::Tag> ()
      .AddConstructor<InspectHeader> ()
      ;
    return tid;
  };
  virtual ns3::TypeId GetInstanceTypeId (void) const {
    return GetTypeId ();
  };
  
  virtual uint32_t GetSerializedSize (void) const {
    return 1 // message type
      + 3*2; // meta data fields
  };
  virtual void Serialize (ns3::Buffer::Iterator i) const {
    i.WriteU8 (messageType);
    i.WriteU16(machineId);
    i.WriteU16(sensorId);
    i.WriteU16(cycleId);
  };
  virtual uint32_t Deserialize (ns3::Buffer::Iterator i) {
    messageType = i.ReadU8();
    machineId   = i.ReadU16();
    sensorId    = i.ReadU16();
    cycleId     = i.ReadU16();
    return GetSerializedSize();
  };
  virtual void Print (std::ostream &os) const {
    os << "InspectHeader{"
       << (int)messageType << ", "
       << machineId << ", "
       << sensorId << ", "
       << cycleId << "}"
      ;
  }

  uint8_t  messageType; // type of message (byte is smallest unit
  uint16_t machineId;   // id. to ORIGIN in original paper
  uint16_t sensorId;    // ..|
  uint16_t cycleId;     // ..|-- together id. to CYCLE ID in original paper
  
  CycleId GetCycleId() const {
    return std::make_tuple(machineId,sensorId,cycleId);
  }
};

struct InspectFeedback : public ns3::Header {
  static ns3::TypeId GetTypeId (void) {
    static ns3::TypeId tid = ns3::TypeId ("ns3::InspectFeedback")
      .SetParent<ns3::Tag> ()
      .AddConstructor<InspectFeedback> ()
      ;
    return tid;
  };
  virtual ns3::TypeId GetInstanceTypeId (void) const {
    return GetTypeId ();
  };
  
  virtual uint32_t GetSerializedSize (void) const {
    return 1 * abs_r(); // feedback vector (1 byte per entry)
  };
  virtual void Serialize (ns3::Buffer::Iterator i) const {
    for( auto fb : feedbackVector ) {
      i.WriteU8(fb);
    }
  };
  virtual uint32_t Deserialize (ns3::Buffer::Iterator i) {
    feedbackVector.resize( abs_r() );
    for(size_t j=0; j< abs_r() ; ++j) {
      feedbackVector.at(j) = i.ReadU8();
    }
    return GetSerializedSize();
  };
  virtual void Print (std::ostream &os) const {
    os << "InspectFeedback";
  }

  std::vector<uint8_t> feedbackVector;
};

struct InspectLinearCombination : public ns3::Header {
  static ns3::TypeId GetTypeId (void) {
    static ns3::TypeId tid = ns3::TypeId ("ns3::InspectLinearCombination")
      .SetParent<ns3::Tag> ()
      .AddConstructor<InspectLinearCombination> ()
      ;
    return tid;
  };
  virtual ns3::TypeId GetInstanceTypeId (void) const {
    return GetTypeId ();
  };
  
  virtual uint32_t GetSerializedSize (void) const {
    return n() + b;
  };
  virtual void Serialize (ns3::Buffer::Iterator i) const {
    NS_ASSERT(encodingVector.size() == n());
    NS_ASSERT(informationVector.size() == b);
    for( auto v : encodingVector ) {
      i.WriteU8(v._poly);
    }
    for( auto v : informationVector ) {
      i.WriteU8(v._poly);
    }
  };
  virtual uint32_t Deserialize (ns3::Buffer::Iterator i) {
    encodingVector.resize(    n() );
    informationVector.resize( b );
    for(size_t j=0; j<n() ; ++j) {
      encodingVector.at(j) = i.ReadU8();
    }
    for(size_t j=0; j<b ; ++j) {
      informationVector.at(j) = i.ReadU8();
    }
    return GetSerializedSize();
  };
  virtual void Print (std::ostream &os) const {
    os << "InspectLinearCombination";
  }

  std::vector<ff::gf_elem> encodingVector;    // size N
  std::vector<ff::gf_elem> informationVector; // size K

  ff::linear_combination GetFfLc () const {
    return {encodingVector, informationVector};
  }
  static InspectLinearCombination FromFfLc(const ff::linear_combination& lc) {
    InspectLinearCombination res;
    res.encodingVector = lc.combination;
    res.informationVector = lc.data;
    return res;
  }
};
