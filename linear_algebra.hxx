#ifndef __LINEAR_ALGEBRA_HXX
#define __LINEAR_ALGEBRA_HXX

/**
 * Implements linear algebra over the finite field GF(2^8)
 * modulus the primitive polynom x^8+x^4+x^3+x+1 (^= 283).
 * (The same polynomial is used by AES)
 * 
 * The speed of multiplication over polynomials is improved with a
 * logarithm and exponent lookup table for the generator x+1 (^= 3).
 * Both lookup tables together need 512 bytes of memory.
 * The tables are initialized automatically upon program start using
 * C++ static member semantics.
 *
 * Supported operations to date are:
 *   - all standard finite field operations
 *   - gaussian elemination
 */

#include "Matrix.h"

#include <stdint.h>
#include <vector>
#include <cassert>
#include <iostream>

namespace ff
{

class gf_elem {
public:
  uint8_t _poly;

  gf_elem():_poly(0){};
  gf_elem(int poly);
  
  bool operator==(const gf_elem& rhs)    const;
  bool operator!=(const gf_elem& rhs)    const;
  
  gf_elem operator+ (const gf_elem& rhs) const;
  gf_elem operator+=(const gf_elem& rhs);
  gf_elem operator- ()                   const;
  gf_elem operator- (const gf_elem& rhs) const;
  gf_elem operator-=(const gf_elem& rhs);
  gf_elem operator*(const gf_elem& rhs)  const;
  gf_elem operator*=(const gf_elem& rhs);
  gf_elem operator/(const gf_elem& rhs)  const;
  gf_elem operator/=(const gf_elem& rhs);

  static const gf_elem zero;
  static const gf_elem one;
private:
  static struct _init { _init(); } _initializer;
  static uint8_t _log_table[256];
  static uint8_t _exp_table[256];
  
  uint8_t static gmul(uint8_t a, uint8_t b);
};

std::ostream& operator<< (std::ostream& stream, const gf_elem& m);

uint16_t gauss(Numeric_lib::Matrix<gf_elem,2>& A);
uint16_t lu_decomposition(Numeric_lib::Matrix<gf_elem,2>& P,
			  Numeric_lib::Matrix<gf_elem,2>& L,
			  Numeric_lib::Matrix<gf_elem,2>& U);

std::vector<gf_elem> lu_solve(const Numeric_lib::Matrix<gf_elem,2>& P,
			      const Numeric_lib::Matrix<gf_elem,2>& L,
			      const Numeric_lib::Matrix<gf_elem,2>& U,
			      const std::vector<gf_elem>& b,
			      uint16_t sub_rank = 0);

void mtx_add_rows(Numeric_lib::Matrix<gf_elem,2>& m, size_t dst, size_t src, gf_elem factor);
void mtx_swap_columns(Numeric_lib::Matrix<gf_elem,2>& m, size_t dst, size_t src);
void mtx_multiply(Numeric_lib::Matrix<gf_elem,2>& m, size_t row, gf_elem factor);
void print_matrix(const Numeric_lib::Matrix<gf_elem,2>& m);

template<typename T>
void print_vector(const std::vector<T>& v) {
  for(size_t i = 0; i < v.size(); ++i) {    printf("| ");
    printf("%3u",v.at(i));
    printf(" |\n");
  }
};
template<>
void print_vector(const std::vector<gf_elem>& v);

template <typename T>
Numeric_lib::Matrix<T,2> operator*(const Numeric_lib::Matrix<T,2>& A, const Numeric_lib::Matrix<T,2>& B)
{
  assert(A.dim2() == B.dim1());
  Numeric_lib::Matrix<T,2> X(A.dim1(), B.dim2());
  for(size_t row = 0; row < A.dim1(); ++row) {
    for(size_t column = 0; column < B.dim2(); ++column) {
      for(size_t inner = 0; inner < A.dim2(); inner++) {
      X(row,column) += A(row,inner) * B(inner,column);
      }
    }
  }
  return X;
}

template <typename T>
bool operator==(const Numeric_lib::Matrix<T,2>& A, const Numeric_lib::Matrix<T,2>& B)
{
  if(A.dim1() != B.dim1() || A.dim2() != B.dim2()) {
    return false;
  }

  for(size_t row = 0; row < A.dim1(); ++row ) {
    for(size_t column = 0; column < A.dim2(); ++column ) {
      if(A(row,column) != B(row,column))
	return false;
    }
  }
  return true;
}

template<typename T>
std::ostream& operator<< (std::ostream& stream, const Numeric_lib::Matrix<T,2>& m)
{
  stream << std::endl;
  for(size_t i = 0; i < m.dim1(); ++i) {
    stream << "| ";
    for(size_t j = 0; j < m.dim2(); ++j) {
      stream << m(i,j);
    }
    stream << "|\n";
  }
  return stream;
}

template<typename T>
std::vector<T> operator*(const Numeric_lib::Matrix<T,2>& A, std::vector<T> b) {
  assert(A.dim2() == b.size());

  std::vector<T> x(A.dim1());
  
  for(size_t row=0; row<A.dim1(); row++) {
    for(size_t column=0; column<A.dim2(); column++) {
      x.at(row) += A(row,column) * b.at(column);
    }
  }
  return x;
}

template<typename T>
Numeric_lib::Matrix<T,2> transpose(const Numeric_lib::Matrix<T,2>& A) {
  Numeric_lib::Matrix<T,2> X(A.dim2(), A.dim2());
  for(size_t i=0; i<A.dim1(); i++) {
    for(size_t j=0; j<A.dim2(); j++) {
      X(j,i) = A(i,j);
    }
  }
  return X;
}
  
};
  
#endif //__LINEAR_ALGEBRA_HXX
