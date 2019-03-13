#include "linear_algebra.hxx"
#include "Matrix.h"

#include <cassert>
#include <cstdio>
#include <iostream>
#include <boost/format.hpp>

/////////////////////////////////////////////////////
//          FINITE FIELD OPERATIONS                //
/////////////////////////////////////////////////////

using ff::gf_elem;

ff::gf_elem::gf_elem(int poly)
{
  assert(poly >= 0);
  assert(poly < 256);
  _poly = poly;
}

bool ff::gf_elem::operator==(const gf_elem& rhs) const
{
  return _poly == rhs._poly;
}
bool ff::gf_elem::operator!=(const gf_elem& rhs) const
{
  return _poly != rhs._poly;
}
gf_elem ff::gf_elem::operator+ (const gf_elem& rhs) const
{
  return _poly ^ rhs._poly;
};
gf_elem ff::gf_elem::operator+=(const gf_elem& rhs)
{
  return _poly ^= rhs._poly;
};
gf_elem ff::gf_elem::operator- () const
{
  return gf_elem::zero - *this;
};
gf_elem ff::gf_elem::operator- (const gf_elem& rhs) const
{
  return _poly ^ rhs._poly;
};
gf_elem ff::gf_elem::operator-=(const gf_elem& rhs)
{
  return _poly ^= rhs._poly;
};
gf_elem ff::gf_elem::operator*(const gf_elem& rhs) const
{
  if(_poly == 0 || rhs._poly == 0)
    return 0;
  else
    return _exp_table[ (_log_table[_poly] + _log_table[rhs._poly]) % 255 ];
}
gf_elem ff::gf_elem::operator*=(const gf_elem& rhs)
{
  if(_poly == 0 || rhs._poly == 0)
    return _poly = 0;
  else
    return _poly = _exp_table[ (_log_table[_poly] + _log_table[rhs._poly]) % 255 ];
}
gf_elem ff::gf_elem::operator/(const gf_elem& rhs) const
{
  assert(rhs._poly > 0 && "error: division by zero in finite field");
  if(_poly == 0) {
    return 0;
  }
  else {
    uint8_t log_lhs = _log_table[_poly];
    uint8_t log_rhs = _log_table[rhs._poly];
    if(log_lhs >= log_rhs) {
      return _exp_table[ log_lhs - log_rhs ];
    }
    else {
      return _exp_table[ 255 - log_rhs + log_lhs ];
    }
  }
}
gf_elem ff::gf_elem::operator/=(const gf_elem& rhs) 
{
  assert(rhs._poly > 0 && "error: division by zero in finite field");
  if(_poly == 0) {
    return 0;
  }
  else {
    uint8_t log_lhs = _log_table[_poly];
    uint8_t log_rhs = _log_table[rhs._poly];
    if(log_lhs >= log_rhs) {
      return _poly = _exp_table[ log_lhs - log_rhs ];
    }
    else {
      return _poly = _exp_table[ 255 - log_rhs + log_lhs ];
    }
  }
}

ff::gf_elem::_init::_init() {
  uint8_t x = 1;
  uint8_t prim = 3;
  for(int i=0;i<255;i++) {
    // invariant here: prim^i = x
    _log_table[x] = i;
    _exp_table[i] = x;
    //printf("%3d^%3u = %3u\n",3,i,x);
	
    x = gmul(x,prim);
  }
  _exp_table[255] = 1;
}

/* Multiply two numbers in the GF(2^8) finite field defined 
 * by the polynomial x^8 + x^4 + x^3 + x + 1 = 0
 * using the Russian Peasant Multiplication algorithm
 * (the other way being to do carry-less multiplication followed by a modular reduction)
 */
uint8_t ff::gf_elem::gmul(uint8_t a, uint8_t b) {
  uint8_t p = 0; /* the product of the multiplication */
  while (b) {
    if (b & 1) /* if b is odd, then add the corresponding a to p (final product = sum of all a's corresponding to odd b's) */
      p ^= a; /* since we're in GF(2^m), addition is an XOR */

    if (a & 0x80) /* GF modulo: if a >= 128, then it will overflow when shifted left, so reduce */
      a = (a << 1) ^ 0x11b; /* XOR with the primitive polynomial x^8 + x^4 + x^3 + x + 1 -- you can change it but it must be irreducible */
    else
      a <<= 1; /* equivalent to a*2 */
    b >>= 1; /* equivalent to b // 2 */
  }
  return p;
}

uint8_t ff::gf_elem::_log_table[256];
uint8_t ff::gf_elem::_exp_table[256];
gf_elem::_init ff::gf_elem::_initializer;
const gf_elem ff::gf_elem::zero(0);
const gf_elem ff::gf_elem::one(1);


std::ostream& ff::operator<< (std::ostream& stream, const gf_elem& m)
{
  stream << boost::format(" %3u") % (unsigned)m._poly;
  return stream;
}

/////////////////////////////////////////////////////
//          MATRIX OPERATIONS                      //
/////////////////////////////////////////////////////


template<>
void ff::print_vector(const std::vector<gf_elem>& v) {
  for(size_t i = 0; i < v.size(); ++i) {
    printf("| ");
    printf("%3u",v.at(i)._poly);
    printf(" |\n");
  }
};
void ff::print_matrix(const Numeric_lib::Matrix<gf_elem,2>& m) {
  for(size_t i = 0; i < m.dim1(); ++i) {
    printf("| ");
    for(size_t j = 0; j < m.dim2(); ++j) {
      printf("%3u ",m(i,j)._poly);
    }
    printf("|\n");
  }
}

void ff::mtx_add_rows(Numeric_lib::Matrix<gf_elem,2>& m, size_t dst, size_t src, gf_elem factor) {
  for(size_t j = 0; j < m.dim2(); ++j) {
    m(dst,j) += m(src,j) * factor;
  }
}
void ff::mtx_multiply(Numeric_lib::Matrix<gf_elem,2>& m, size_t row, gf_elem factor) {
  for(size_t j = 0; j < m.dim2(); ++j) {
    m(row,j) *= factor;
  }
}

void ff::mtx_swap_columns(Numeric_lib::Matrix<gf_elem,2>& m, size_t dst, size_t src)
{
  gf_elem tmp;
  for(size_t i=0; i<m.dim1(); i++) {
    tmp = m(i,dst);
    m(i,dst) = m(i,src);
    m(i,src) = tmp;
  }
}

// warning: only works with finite fields!!! (no abs used at maximum-search)
uint16_t ff::gauss(Numeric_lib::Matrix<gf_elem,2>& A)
{
  // printf("original: \n"); print_matrix(A); printf("\n");

  int n = A.dim1();

  for (int i=0; i<n; i++) {
    // Search for maximum in this column
    int nzRow = -1;
    for (int k=i; k<n; k++) {
      if (A(k,i) != gf_elem::zero) {
	nzRow = k;
	break;
      }
    }
    if(nzRow == -1) {
      // no non-zero row element found
      // we cannot solve the system
      return i;
    }

    A.swap_rows(nzRow,i);

    // Make all rows below this one 0 in current column
    for (int k=i+1; k<n; k++) {
      gf_elem c = -A(k,i)/A(i,i);
      mtx_add_rows(A, k, i, c);
    }
    // printf("partial ecolon form: \n"); print_matrix(A); printf("\n");
  }

  // Solve equation Ax=b for an upper triangular matrix A
  // printf("ecolon form: \n"); print_matrix(A); printf("\n");
  // printf("...starting solve: \n");

  std::vector<gf_elem> x(n);
  for (int i=n-1; i>=0; i--) {
    mtx_multiply(A,i,gf_elem::one/A(i,i));
    //x.at(i) = A(i,n);
    // printf("step A: \n"); print_matrix(A); printf("\n");
    for(int k=i-1;k>=0; k--) {
      mtx_add_rows(A,k,i,-A(k,i));
    }
    // printf("step B: \n"); print_matrix(A); printf("\n");
  }
  return n;
}


// warning: only works with finite fields!!! (no abs used at maximum-search)
uint16_t ff::lu_decomposition(Numeric_lib::Matrix<gf_elem,2>& P, Numeric_lib::Matrix<gf_elem,2>& L, Numeric_lib::Matrix<gf_elem,2>& U)
{
  // printf("original: \n"); print_matrix(U); printf("\n");
  //uint32_t degree = 0;

  int n = U.dim1();

  // pivotization
  for (int i=0; i<n; i++) {
    // search for non-zero element in this column
    int nzRow = -1;
    for (int k=i; k<n; k++) {
      if (U(k,i) != gf_elem::zero) {
	nzRow = k;
	break;
      }
    }
    if(nzRow == -1) {
      // no non-zero row element found
      // we cannot continue to solve the system
      // there may be other non-zero columns later (which increases the degree)
      for(int column_offset=0; i+column_offset<n; ++column_offset) {
	if(U(i,i+column_offset) != gf_elem::zero) {
	  return i+1;
	}
      }
      return i;
    }
  
    // printf("...swapping rows %u & %u\n",i,nzRow);
    U.swap_rows(nzRow,i);
    P.swap_rows(nzRow,i);
    for(int k=0;k<i;k++) {
      gf_elem tmp;
      // printf("swapping L(%u,%u) <-> L(%u,%u) [%3u <-> %3u]",i,k,nzRow,i,L(i,k)._poly,L(nzRow,i)._poly);
      tmp = L(i,k);
      L(i,k) = L(nzRow,k);
      L(nzRow,k) = tmp;
    }
    
    // Make all rows below this one 0 in current column
    for (int k=i+1; k<n; k++) {
      gf_elem c = -U(k,i)/U(i,i);
      if(c == gf_elem::zero) {
	continue;
      }
      mtx_add_rows(U, k, i, c);
      // fill in corresponding L entry
      // printf("...writing %3u at L(%u,%u) [== %u]\n", c._poly, k, i, L(k,i)._poly);
      L(k,i) = c;
    }
    // printf("partial echolon U: \n"); print_matrix(U); printf("\n");
    // printf("partial echolon L: \n"); print_matrix(L); printf("\n");
  }

  return n;
}

std::vector<gf_elem> ff::lu_solve(const Numeric_lib::Matrix<gf_elem,2>& P,
				  const Numeric_lib::Matrix<gf_elem,2>& L,
				  const Numeric_lib::Matrix<gf_elem,2>& U,
				  const std::vector<gf_elem>& b,
				  uint16_t sub_rank)
{
  int n = P.dim1();
  assert(P.dim1() == P.dim2());
  assert(L.dim1() == L.dim2());
  assert(U.dim1() == U.dim2());

  assert(L.dim1() == (size_t)n);
  assert(U.dim1() == (size_t)n);
  
  // Solve equation Ux=b for an upper triangular matrix U
  // printf("P: \n"); print_matrix(P); printf("\n");
  // printf("L: \n"); print_matrix(L); printf("\n");
  // printf("U: \n"); print_matrix(U); printf("\n");
  // printf("b: ");
  // for(int i=0;i<b.size(); ++i) {
  //   printf("%3u,",b.at(i)._poly);
  // }
  // printf("\n");
  // printf("...starting solve: \n");

  if(!sub_rank) {
    sub_rank = n;
  }
  
  std::vector<gf_elem> b_ = P*b;
  b_.resize(sub_rank);
  std::vector<gf_elem> y(sub_rank);
  std::vector<gf_elem> x(sub_rank);
  // solve L*y = b for y  where y := U*x
  for (int i=0; i<sub_rank; i++) {
    y.at(i) = b_.at(i); // since L(i,i) == 1 forall i in 0..n-1
    for(int k=0;k<i; k++) {
      // printf("y(%u) -= %3u*%3u;\n",i,L(i,k)._poly,y.at(k)._poly);
      y.at(i) -= L(i,k)*y.at(k);
    }
    // printf("y: ");
    // for(int i=0;i<y.size(); ++i) {
    //   // printf("%3u,",y.at(i)._poly);
    // }
    // printf("\n");
  }
  for (int i=sub_rank-1; i>=0; i--) {
    x.at(i) = y.at(i);
    for(int k=i+1;k<sub_rank; k++) {
      // printf("x(%u) -= %3u*%3u;\n",i,U(i,k)._poly,x.at(k)._poly);
      x.at(i) -= U(i,k)*x.at(k);
    }
    x.at(i) /=  U(i,i);
    // printf("x: ");
    // for(int i=0;i<x.size(); ++i) {
    //   printf("%3u,",x.at(i)._poly);
    // }
    // printf("\n");
  }
  return x;
}
