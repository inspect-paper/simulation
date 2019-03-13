#ifndef __NETWORK_CODING_HXX
#define __NETWORK_CODING_HXX

#include "linear_algebra.hxx"
#include <map>
#include <random>

namespace ff {

struct linear_combination
{
  std::vector<gf_elem> combination;
  std::vector<gf_elem> data;
};

class gf_rng
{
public:
  virtual gf_elem roll() = 0;
  virtual ~gf_rng(){};
};

class Feedback
{
public: 
	uint32_t node;
	std::vector<uint16_t> fb;
	std::vector<uint16_t> addedfb;
	std::vector<uint16_t> assumedaddedfb;
};

class nc_general
{
public:
  virtual ~nc_general() {}
  nc_general()
    :completed(false)
    ,counter(0) 
  {}
  virtual void saveFB(uint32_t, const std::vector<uint16_t>&);
  virtual void delFB(uint32_t);
  virtual int classChoice(int);
  virtual void checkComplete() = 0;
  virtual uint16_t get_current_degree() const = 0;
  bool completed;
  int counter;
protected:
  std::vector<Feedback> _feedbackList;
  std::vector<std::vector<uint8_t> >* _datasets;
  std::vector<size_t> _classes;
};
class nc_generator : public nc_general
{
public:
  nc_generator(gf_rng* rng, std::vector<std::vector<uint8_t> >* datasets, const std::vector<size_t>& R);
  linear_combination new_linear_combination() const;
  linear_combination new_zerotail_linear_combination(size_t first_zero_idx) const;
  virtual uint16_t get_current_degree() const;
  virtual void checkComplete();
  int generationNumber;
private:
  gf_rng* _rng;
  std::vector<std::vector<uint8_t> >* _datasets;
};

class abstract_nc_absorber : public nc_general 
{
public:
  virtual ~abstract_nc_absorber(){};

  virtual bool inject(const linear_combination& lc) = 0;

};

class nc_absorber : public abstract_nc_absorber
{
	
public:
  nc_absorber(gf_rng* rng, uint16_t size_of_generation);
  virtual ~nc_absorber(){}

  virtual bool inject(const linear_combination& lc);
  virtual linear_combination new_zerotail_linear_combination(size_t total, size_t data) const;
  virtual std::vector<std::vector<uint8_t> > solve();
  virtual std::vector<std::vector<uint8_t> > partial_solve(uint16_t sub_rank);

  virtual uint16_t get_current_degree() const {return _degree;}
  virtual void     checkComplete() {assert(false && "not implemented");}

  uint16_t get_max_degree() const {return _P.dim1();};

  
  std::vector<linear_combination>* _receivedLC;
  
private:
  gf_rng* _rng;
  uint16_t _degree;
  Numeric_lib::Matrix<gf_elem,2> _P;
  Numeric_lib::Matrix<gf_elem,2> _L;
  Numeric_lib::Matrix<gf_elem,2> _U;
  std::vector<std::vector<gf_elem> > _bs;

  
  void perform_swap(size_t idx_a, size_t idx_b);
  uint16_t lu_decomposition();
};

// pseudo-random rng with fixed seed
class simple_rng : public gf_rng {
public:
  simple_rng()
    :engine(1337)
    ,dist(1,255)
  {}

  virtual gf_elem roll() {
    return dist(engine);
  }
private:
  std::mt19937 engine;
  std::uniform_int_distribution<short> dist;
};

class naive_nc_absorber : public abstract_nc_absorber {
public:
  naive_nc_absorber(gf_rng* rng, const std::vector<size_t>& R);
  virtual ~naive_nc_absorber(){};

  virtual bool inject(const linear_combination& lc);

  virtual uint16_t get_current_degree() const;
  virtual std::vector<uint16_t>  get_current_degrees() const;
  linear_combination new_zerotail_linear_combination(int classindex, size_t total, size_t data);
  virtual void checkComplete();
  void saveLinComb(const linear_combination& lc);

private:
  std::vector<nc_absorber> absorbers;
  uint16_t maximum_solvable;
 };

};

#endif //__NETWORK_CODING_HXX
