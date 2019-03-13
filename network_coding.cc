/* -*- Mode: C++; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*- */
#include "network_coding.hxx"
#include <cassert>
#include <algorithm>
#include <numeric>
#include <cstdio>

#define UNIMPLEMENTED() assert(false && "Not yet implemented.")

using namespace ff;

nc_generator::nc_generator(gf_rng* rng, std::vector<std::vector<uint8_t>>* datasets, const std::vector<size_t>& R)
  :_rng(rng)
  ,_datasets(datasets)
{
	_classes = R;
	std::vector<Feedback> _feedbackList;
	completed = false;
#ifndef NDEBUG	
  size_t n = _datasets->at(0).size();
  for(auto v : *_datasets) {
    assert(v.size() == n);
  }
#endif
}

linear_combination nc_generator::new_linear_combination() const
{
  size_t size_of_generation = _datasets->size();
  size_t size_of_dataset    = _datasets->at(0).size();

  linear_combination lc;
  lc.combination.resize(size_of_generation);
  lc.data.resize(size_of_dataset);

  for(size_t i=0; i<size_of_generation; i++) {
    lc.combination.at(i) = _rng->roll();
    for( size_t j=0; j<size_of_dataset; j++) {
      lc.data.at(j) += lc.combination.at(i) * _datasets->at(i).at(j);
    }
  }

  return lc;
}

uint16_t nc_generator::get_current_degree() const
{
  return _datasets->size();
}
linear_combination nc_generator::new_zerotail_linear_combination(size_t first_zero_idx) const
{
  size_t size_of_generation = _datasets->size();
  size_t size_of_dataset    = _datasets->at(0).size();
  assert(first_zero_idx <= size_of_generation);
  assert(first_zero_idx > 0);

  linear_combination lc;
  lc.combination.resize(size_of_generation);
  lc.data.resize(size_of_dataset);

  for(size_t i=0; i<first_zero_idx; i++) {
    lc.combination.at(i) = _rng->roll();
    for(size_t j=0; j<size_of_dataset; j++) {
      lc.data.at(j) += lc.combination.at(i) * _datasets->at(i).at(j);
    }
  }
	  
  return lc;
}
/*
 * check if all neighbors have completely received generation
 * if so increase generation counter 
 */ 

void nc_generator::checkComplete(){
	//check if all feedback messages are up-to-date
	bool complete = true;
	//some feedback must have been received
	if(!((int)_feedbackList.size() == 0)){
		for(size_t j = 0; j < _feedbackList.size(); j++){
			if(!(_feedbackList.at(j).addedfb.at(_classes.size()-1) == _classes.at(_classes.size()-1))){
				complete = false;
				break;
			}
		}
		//case complete 
		if(complete){
			completed = true;
		}
	}	
}
/*
 * node - ID of node 
 * addedfb - added feedback
 * 
 * update in feedback list or add new entry  
 */
void nc_general::saveFB(uint32_t node, const std::vector<uint16_t>& addedfb){
	std::vector<uint16_t> fb;	
	for(size_t k =0; k < _classes.size(); k++){	
		//compute fb vector
		if(k==0){
			fb.push_back(addedfb.at(k));
		}
		else{
			fb.push_back(addedfb.at(k)-addedfb.at(k-1));
		}
	}
	
	bool found = false;	
	for(size_t i=0; i<_feedbackList.size(); i++){
		//update feedback
		if(_feedbackList.at(i).node == node){
			found= true;			
			_feedbackList.at(i).fb= fb;
			_feedbackList.at(i).assumedaddedfb= addedfb;
			_feedbackList.at(i).addedfb= addedfb;
			break;
		}
	}
	//new feedback entry
	if(!found){
		Feedback newentry;
		newentry.node = node;
		newentry.fb = fb;
		newentry.addedfb= addedfb;
		newentry.assumedaddedfb = addedfb;
		_feedbackList.push_back(newentry);
	}
}

/*
 * delete feedback from node with ID given in input: node  
 */ 
void nc_general::delFB(uint32_t node){
	for(size_t i =0; i < _feedbackList.size(); i++){
		if(_feedbackList.at(i).node == node){
			_feedbackList.erase(_feedbackList.begin()+i);
			break;
		}
	}
}

/*
 * m_kende - range how many layers are candidates
 * returns computed layer
 * realizes algorithm iNsPECt
 * RTord - RT_{Ord}, 
 * RTsc - RT_{SL}, 
 * DCord - DC_{Ord}, 
 * DCsc - DC_{SL}, 
 * QRT - Q_{rt},
 * old.. value that is still required from one of the last iterations
 * ...X - values that only compute somethin for the feedback of a single neighbor
 * choice - layer
 */ 
int nc_general::classChoice(int m_kende)
{		
	//no feedback yet, return first
	if(_feedbackList.size() == 0){
		return 0;
	}
	int startvar = _classes.size()-1;
	for(size_t i=0; i <_feedbackList.size(); i++){
		
		int candidate = startvar;
		int j = startvar;
		while(_feedbackList.at(i).assumedaddedfb.at(j) < _classes.at(j)){		
			candidate = j;
			if(j == 0){
				break;
			}				
			j--; 

		}
		if(candidate < startvar){
			startvar= candidate;
		}
	}	
	int ende = _classes.size();
	if( (startvar+ m_kende) < (int) _classes.size()){
		ende = startvar+m_kende;
	}
	
	int RTsc=0;
	int RTscX=0;
	int oldRTscX=0;
	int RTord =0;
	int forRTord=0;
	
	int DCscX =0;
	int DCsc =0;
	int DCord=0;
	
	int QRT;
	int oldQRT=0;
	int choice = startvar;
	
	//save single DCordX
	std::vector<int> pta;
	for(int x=0; x < (int) _feedbackList.size(); x++ ){
		pta.push_back(0);
	}
	
	for(int i = 0; i< ende; i++){
		DCsc = 0;		
		for(int x =0; x < (int) _feedbackList.size(); x++){
			for(int j=i; j < (int) _classes.size(); j++){
				//RTscX
				if(j ==i){
					RTscX = _classes.at(j) - _feedbackList.at(x).assumedaddedfb.at(j);
					//compare old and new value
					oldRTscX = 0;
					if(j >0){
						oldRTscX = _classes.at(j-1) - _feedbackList.at(x).assumedaddedfb.at(j-1);
						if(RTscX < oldRTscX){
							oldRTscX = RTscX;
						}
					}
				}
				else{
					int zwisch = _classes.at(j) - _feedbackList.at(x).assumedaddedfb.at(j);
					if(zwisch < RTscX){
						RTscX = zwisch;
					}
					if(zwisch < oldRTscX){
						oldRTscX = zwisch;
					}
				}
			}
			int multiple = 0;
			for(int j=i; j >= 0; j--){
				if(_classes.at(j) == _feedbackList.at(x).assumedaddedfb.at(j)){
					break;
				}
				else{
					multiple ++;
				} 
			}
			
			// RTsc + for RTord
			if(x ==0){
				RTsc = RTscX;
				forRTord = RTscX -oldRTscX;
			}
			else{
				if(RTsc < RTscX){
					RTsc = RTscX;
				}
				if(forRTord < (RTscX -oldRTscX)){
					forRTord = RTscX -oldRTscX;
				}
			}
		    DCscX = (multiple) * RTscX;
		    
		    //from here on: sum up 
		    if(i==0){
				pta.at(x) = RTscX;
			}
			else{
				if(RTscX != oldRTscX){
					pta.at(x) = RTord + RTscX - oldRTscX;
				}
			}
			DCord += pta.at(x);
			DCsc += DCscX;
		}
		RTord = RTord + forRTord;
		QRT= RTord - RTsc; 
		if(i == startvar){
			oldQRT = QRT;
		}
		else{
			if((i>startvar) && (QRT > oldQRT) && ((DCord - DCsc)>=0)){
				choice =i;
				oldQRT = QRT;
				
			}
		}
	}
	
	for(size_t i=0; i <_feedbackList.size(); i++){
		for(int j = choice; j < (int )_classes.size(); j++){
			_feedbackList.at(i).assumedaddedfb.at(j)++;
		}
	}
	return choice;
}

nc_absorber::nc_absorber(gf_rng* rng, uint16_t size_of_generation)
  :_rng(rng)		 
  ,_degree(0)
  ,_P(size_of_generation,size_of_generation)
  ,_L(size_of_generation,size_of_generation)
  ,_U(size_of_generation,size_of_generation)
  ,_bs(size_of_generation)
{
  for(int i=0; i<size_of_generation; i++) {
    _P(i,i) = 1;
    _L(i,i) = 1;
  }
}

// warning: only works with finite fields!!! (no abs used at maximum-search)
uint16_t nc_absorber::lu_decomposition()
{
  int n = _U.dim1();

  // pivotization
  for (int i=0; i<n; i++) {
#ifdef NC_DEBUG
    printf("pivotization for i=%d\n",i);
#endif
    // search for non-zero element in this column
    int nzRow = -1;
    for (int k=i; k<n; k++) {
#ifdef NC_DEBUG
      printf("  checking k=%d, _U(k,i)=%u\n",k,_U(k,i)._poly);
#endif
      if (_U(k,i) != gf_elem::zero) {
	nzRow = k;
#ifdef NC_DEBUG
      printf("  breaking (found != zero)\n");
#endif
	break;
      }

    }
    if(nzRow == -1) {
      // no non-zero row element found
      // we cannot continue to solve the system
      // there may be other non-zero columns later (which increases the degree)
      for(int column_offset=0; i+column_offset<n; ++column_offset) {
	if(_U(i,i+column_offset) != gf_elem::zero) {
	  return i+1;
	}
      }
      return i;
    }
    
    if(nzRow != i) {
#ifdef NC_DEBUG
      printf("...swapping rows %u & %u\n",i,nzRow);
#endif
      _U.swap_rows(nzRow,i);
      _P.swap_rows(nzRow,i);
      for(int k=0;k<i;k++) {
	gf_elem tmp;
	std::swap(_L(i,k), _L(nzRow,k));
      }
    }
    
    // reduce all cells below current cell in this column to zero
    for (int k=i+1; k<n; k++) {
      gf_elem c = -_U(k,i)/_U(i,i);
      if(c == gf_elem::zero) {
	continue;
      }
      mtx_add_rows(_U, k, i, c);
      // fill in corresponding _L entry
      _L(k,i) = c;
    }
  }

  return n;
}

inline void counter_eliminate(Numeric_lib::Matrix<gf_elem,2>& U,
			      Numeric_lib::Matrix<gf_elem,2>& L,
			      size_t row,
			      size_t column)
{
  gf_elem c = L(row,column);
#ifdef NC_DEBUG
  gf_elem old = U(row,column);
#endif
  if(c != gf_elem::zero) {
    mtx_add_rows(U,row,column,c);
    L(row,column) = 0;
  }
#ifdef NC_DEBUG
  printf("counter eliminating (%lu,%lu): %u -> %u with %u*%u\n",
	 row,column,old._poly,U(row,column)._poly,c._poly,U(column,column)._poly);
#endif
}

void nc_absorber::perform_swap(size_t idx_a, size_t idx_b)
{
  // sanity check prep work
#ifndef NDEBUG
  auto check_mtx = transpose(_P) * _L * _U;
#endif
#ifdef NC_DEBUG
  printf(" --------------<> SWAP OP\n");
  printf("check mtx:\n");
  print_matrix(check_mtx);
  printf("U_x:\n");
  print_matrix(_U);
  printf("L_x:\n");
  print_matrix(_L);
  printf("P_x:\n");
  print_matrix(_P);
  printf("\n");
#endif
  size_t greater = std::max(idx_a,idx_b);
  size_t lower = std::min(idx_a,idx_b);
  
  // printf("ztr: %u, degree: %u, swap_row: %zu\n", zerotail_rank, _degree, swap_row);
  if(greater == lower) {
    return;
  }

  for(int row = _L.dim1() - 1; row > (int)greater; row--) {
    counter_eliminate(_U,_L,row,lower);
   }
  
  // swap in upper matrix
  _U.swap_rows(greater, lower);
  // swap in pivot matrix
  _P.swap_rows(greater, lower);
  // swap in lower matrix
  for(size_t k=0;k<lower;k++) {
    std::swap(_L(lower,k), _L(greater,k));
  }

  
#ifdef NC_DEBUG
  printf("P^T * L * U:\n");
  print_matrix(transpose(_P) * _L * _U);
  printf("U_x:\n");
  print_matrix(_U);
  printf("L_x:\n");
  print_matrix(_L);
  printf("P_x:\n");
  print_matrix(_P);
  printf("\n");
#endif
  assert( check_mtx == transpose(_P)*_L*_U );
#ifdef NC_DEBUG
  printf(" <>---------------- END SWAP OP\n");
#endif
}

template<typename T>
std::vector<T> gen_maxvec(const std::vector<T>& v) {
  std::vector<T> maxvec(v.size());
  maxvec.at(0) = v.at(0);
  for(size_t i=1; i<v.size(); i++) {
    maxvec.at(i) = std::max(v.at(i),maxvec.at(i-1));
  }
  return maxvec;
}

bool nc_absorber::inject(const linear_combination& lc)
{
#ifdef NC_DEBUG
  printf("INJECTING: \n");
  print_vector(lc.combination);
#endif
  if(_degree == _U.dim1()) {
    // already fully determined
    return true;
  }

  std::vector<gf_elem> v( _U.dim1() );
  std::iota(v.begin(), v.end(), 0);
  uint8_t idx = (_P * v).at(_degree)._poly;
  _bs.at(idx) = lc.data;
  
  // write coefficients into upper echolon matrix
  for(size_t i = 0; i < _U.dim2(); ++i) {
    _U(_degree,i) = lc.combination.at(i);
  }
  for(size_t i = 0; i < _degree; ++i) {
    _L(_degree,i) = gf_elem::zero;
  }

  

#ifdef NC_DEBUG
  printf("P^T * L * U:\n");
  print_matrix(transpose(_P) * _L * _U);
  printf("U1:\n");
  print_matrix(_U);
  printf("\n");
  printf("L1:\n");
  print_matrix(_L);
  printf("P1:\n");
  print_matrix(_P);
  printf("\n");
#endif
  
  // check if combination was linear dependent
#ifndef NDEBUG
  auto check_mtx = transpose(_P) * _L * _U;
#endif  
  uint16_t new_degree = this->lu_decomposition();
#ifdef NC_DEBUG
  printf("P^T * L * U:\n");
  print_matrix(transpose(_P) * _L * _U);
  printf("U2:\n");
  print_matrix(_U);
  printf("L2:\n");
  print_matrix(_L);
  printf("P2:\n");
  print_matrix(_P);
  printf("\n");
#endif
  assert( check_mtx == transpose(_P) * _L * _U);
  if(new_degree > _degree)
    {
    }
  else {
    // remove obsolete entries in L
    for(size_t i=0; i<_degree; i++) {
      _L(_degree,i) = gf_elem::zero;
    }

  }
  _degree = new_degree;
  
  if(_degree == _U.dim1()) {
    // found solution
    return true;
  }
  else {
    return false;
  }
}

linear_combination nc_absorber::new_zerotail_linear_combination(size_t total, size_t data) const
{ 
  linear_combination out;
  out.combination.resize(total);
  out.data.resize(data);
  for(size_t i =0; i < total; i++){
	  out.combination.at(i)=0;
  }
  for(size_t i =0; i < data; i++){
	  out.data.at(i)=0;
  }  
  for(size_t i=0; i< _receivedLC->size(); i++) {
	  gf_elem combi;
		while((combi = _rng->roll()) == gf_elem::zero )
			; // a zero factor is fatal here when only one combination is in the list!
		for(size_t j=0; j< _receivedLC->at(i).combination.size(); j++){
			//data
			//combination
			out.combination.at(j) += combi * _receivedLC->at(i).combination.at(j);
	  }
		for(size_t j=0; j< data; j++){
			out.data.at(j)+= combi * _receivedLC->at(i).data.at(j);
		} 
  }
  return out;
}

std::vector<std::vector<uint8_t> > nc_absorber::solve() {
  return partial_solve(_U.dim1());
}

std::vector<std::vector<uint8_t> > nc_absorber::partial_solve(uint16_t sub_rank) {
  assert(sub_rank <= _degree);
  
  size_t data_len = _bs.at(0).size();
  std::vector<std::vector<uint8_t> > solutions(sub_rank);
  for(size_t j=0; j<sub_rank; j++) {
    solutions.at(j).reserve(data_len);
  }
  for(size_t d=0; d<data_len; d++) {
    std::vector<uint8_t> solution;
    std::vector<gf_elem> b(_U.dim1());
    for(size_t i=0; i<_U.dim1(); i++) {
      if(_bs.at(i).begin() == _bs.at(i).end()) {
	break;
      }
      b.at(i) = _bs.at(i).at(d);
    }
    size_t j=0;
    
    for( gf_elem s : lu_solve(_P,_L,_U,b,sub_rank)) {
      solutions.at(j++).push_back(s._poly);
    }
  }
  return solutions;
}


naive_nc_absorber::naive_nc_absorber(gf_rng* rng, const std::vector<size_t>& R)
  :maximum_solvable(0)
{
  _classes=R;
  std::vector<Feedback> _feedbackList;	
  for(size_t n : R) {
    if(n > 1) {
      absorbers.push_back(nc_absorber(rng, n));
    }
  }
  for(size_t n = 0; n < absorbers.size(); n++){
	  absorbers.at(n)._receivedLC = new std::vector<linear_combination>();
	  
  }
  completed=false;
}

linear_combination naive_nc_absorber::new_zerotail_linear_combination(int classindex, size_t total, size_t data){
	linear_combination lc;
	for(size_t i = classindex; i< absorbers.size(); i++){
		if(absorbers.at(i)._receivedLC->size() != 0){
			lc = absorbers.at(i).new_zerotail_linear_combination(total, data);
			return lc;
		}
	}
	assert(false && "could not generate linear combination of given layer");
	return lc;
}

bool naive_nc_absorber::inject(const linear_combination& lc)
{
  // detect zerotail-rank
  uint32_t zerotail_rank = lc.combination.size();
  for(int i = lc.combination.size()-1; i >= 0; --i){
	  if(lc.combination.at(i) != gf_elem::zero){
		  zerotail_rank = i+1;
		  break;
	  }
  }
  assert(zerotail_rank > 0);

  bool absfin = false;
  for(size_t i=0; i<absorbers.size(); i++) {
    if(zerotail_rank <= absorbers.at(i).get_max_degree()) {
      uint16_t degalt = absorbers.at(i).get_current_degree();
      bool fin = absorbers.at(i).inject(lc);
      //add lc
      uint16_t new_deg = absorbers.at(i).get_current_degree();
      if(degalt < new_deg){
				absorbers.at(i)._receivedLC->push_back(lc);
			}
      if(fin) {
				maximum_solvable = std::max(maximum_solvable,new_deg);
				absfin = true;
      }
      else {
				absfin = false;
      }
    }
  }
  return absfin;
}

uint16_t naive_nc_absorber::get_current_degree() const {
  uint16_t max_cur_deg = 0;
  for(size_t i=0; i<absorbers.size(); i++) {
    max_cur_deg = std::max(max_cur_deg, absorbers.at(i).get_current_degree());
  }
  return max_cur_deg;
};

std::vector<uint16_t> naive_nc_absorber::get_current_degrees() const{
  std::vector<uint16_t> max_cur_deg;
  for(size_t i=0; i<absorbers.size(); i++) {
    max_cur_deg.push_back(absorbers.at(i).get_current_degree());
  }
  return max_cur_deg;
};


void naive_nc_absorber::saveLinComb(const linear_combination& lc){
	//inject
	inject(lc);

};

/*
 * check in received feedback if neighbors completed transmission
 * if so set completed value to true
 */
void naive_nc_absorber::checkComplete(){
	bool complete = true;
	//feedbacklist must not be empty!
	if(!(_feedbackList.empty())){
		for(size_t j = 0; j < _feedbackList.size(); j++){
			if(!(_feedbackList.at(j).addedfb.at(_classes.size()-1) == _classes.at(_classes.size()-1))){
				complete = false;
				break;
			}
		}
		//completed!
		if(complete && (get_current_degrees().at(_classes.size()-1) ==_classes.at(_classes.size()-1))){
			completed = true;
		} else {
			completed = false;
		}
	}
}
