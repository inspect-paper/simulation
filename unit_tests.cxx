#include "linear_algebra.hxx"

using namespace ff;

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE( "Some finite element checks for multiplication and division", "[gf_some]")
{
  REQUIRE( gf_elem(1)*1    == 1   );
  REQUIRE( gf_elem(1)/1    == 1   );
  REQUIRE( gf_elem(7)*3    == 9   );
  REQUIRE( gf_elem(9)/3    == 7   );
  REQUIRE( gf_elem(9)/7    == 3   );
  REQUIRE( gf_elem(7)*21   == 107 );
  REQUIRE( gf_elem(107)/21 == 7   );
  REQUIRE( gf_elem(107)/7  == 21  );
  REQUIRE( gf_elem(0)*1    == 0   );
  REQUIRE( gf_elem(1)*1    == 1   );
  REQUIRE( gf_elem(0)+1    == 1   );
}

TEST_CASE( "Finite elements: associativity of multiplication", "[gf_assoc][!hide]" )
{
  for(int i=0; i<256; ++i) {
    for(int j=0; j<256; ++j) {
      for(int k=0; k<256; ++k) {
	// CAPTURE( i );
	// CAPTURE( j );
	// CAPTURE( k );
	REQUIRE( (gf_elem(i)*gf_elem(j))*gf_elem(k) == gf_elem(i)*(gf_elem(j)*gf_elem(k)) );
      }
    }
  }
}
//~ 
TEST_CASE( "Finite elements: distributivity", "[gf_dist]" )
{
  for(int i=0; i<256; ++i) {
    for(int j=0; j<256; ++j) {
      for(int k=0; k<256; ++k) {
	// CAPTURE( i );
	// CAPTURE( j );
	// CAPTURE( k );
	REQUIRE( (gf_elem(i)+gf_elem(j))*gf_elem(k) == gf_elem(i)*gf_elem(k)+gf_elem(j)*gf_elem(k) );
      }
    }
  }
}

TEST_CASE( "Matrix-vector multiplication", "[gf_mat_vec_mul]" )
{
  Numeric_lib::Matrix<gf_elem,2> A(4,3);
  A(0,0) = 1;
  A(1,1) = 1;
  A(2,2) = 1;
  A(3,0) = 1;
  A(3,2) = 1;

  std::vector<gf_elem> v(3);
  v[0] = 1;
  v[1] = 2;
  v[2] = 3;

  std::vector<gf_elem> v2 = A*v;

  REQUIRE( v2.at(0) == v.at(0) );
  REQUIRE( v2.at(1) == v.at(1) );
  REQUIRE( v2.at(2) == v.at(2) );
}

TEST_CASE( "FF-Gauss, 4x4 matrix", "[gf_gauss_4x4]" )
{
  using namespace Numeric_lib;
  Matrix<gf_elem,2> system(4,5);
  uint16_t degree;
  
  degree = gauss(system);

  REQUIRE(degree == 0);

  system(0,0) = gf_elem(1);
  system(0,1) = gf_elem(1);
  system(0,2) = gf_elem(1);
  system(0,3) = gf_elem(1);
  system(0,4) = gf_elem(2);
  degree = gauss(system);

  REQUIRE(degree == 1);

  system(1,0) = gf_elem(1);
  system(1,1) = gf_elem(2);
  system(1,2) = gf_elem(2);
  system(1,3) = gf_elem(1);
  system(1,4) = gf_elem(9);
  degree = gauss(system);

  REQUIRE(degree == 2);

  system(2,0) = gf_elem(3);
  system(2,1) = gf_elem(2);
  system(2,2) = gf_elem(4);
  system(2,3) = gf_elem(2);
  system(2,4) = gf_elem(3);
  degree = gauss(system);

  REQUIRE(degree == 3);

  system(3,0) = gf_elem(3);
  system(3,1) = gf_elem(2);
  system(3,2) = gf_elem(4);
  system(3,3) = gf_elem(1);
  system(3,4) = gf_elem(9);
  degree = gauss(system);

  REQUIRE(degree == 4);

  REQUIRE( system(0,4) == 244 );
  REQUIRE( system(1,4) ==  85 );
  REQUIRE( system(2,4) == 165 );
  REQUIRE( system(3,4) ==   6 );
}

TEST_CASE( "FF-LU-Decomposition 3x3 matrix", "[gf_lu_decomp_3x3]" )
{
  using namespace Numeric_lib;
  Matrix<gf_elem,2> A(3,3);
  Matrix<gf_elem,2> TMP(3,3);
  Matrix<gf_elem,2> U(3,3);
  Matrix<gf_elem,2> L(3,3);
  Matrix<gf_elem,2> P(3,3);
  for(size_t i=0; i<3; ++i) {
    L(i,i) = 1;
    P(i,i) = 1;
  }

  uint16_t degree;

  degree = lu_decomposition(P,L,U);
  REQUIRE(degree == 0);
  TMP = L*U;
  REQUIRE( TMP == P*A );

  U(0,0) = gf_elem(2);
  U(0,1) = gf_elem(4);
  U(0,2) = gf_elem(7);
  A=U;

  degree = lu_decomposition(P,L,U);
  REQUIRE(degree == 1);
  TMP = L*U;
  REQUIRE( TMP == P*A );

  U(1,0) = gf_elem(13);
  A(1,0) = gf_elem(13);;
  U(1,1) = gf_elem(5);
  A(1,1) = gf_elem(5);
  U(1,2) = gf_elem(99);
  A(1,2) = gf_elem(99);

  degree = lu_decomposition(P,L,U);
  REQUIRE(degree == 2);
  TMP = L*U;
  REQUIRE( TMP == P*A );

  U(2,0) = gf_elem(1);
  A(2,0) = gf_elem(1);
  U(2,1) = gf_elem(17);
  A(2,1) = gf_elem(17);
  U(2,2) = gf_elem(240);
  A(2,2) = gf_elem(240);

  degree = lu_decomposition(P,L,U);
  REQUIRE(degree == 3);
  TMP = L*U;
  REQUIRE( TMP == P*A );

  std::vector<gf_elem> x(3);
  x.at(0) = gf_elem(3);
  x.at(1) = gf_elem(4);
  x.at(2) = gf_elem(5);

  // solve for x
  std::vector<gf_elem> b(3);
  b.at(0) =  x[0]* 2 + x[1]* 4 +  x[2]*  7; 
  b.at(1) =  x[0]*13 + x[1]* 5 +  x[2]* 99;
  b.at(2) =  x[0]* 1 + x[1]*17 +  x[2]*240;
  std::vector<gf_elem> x_ = lu_solve(P,L,U,b);

  REQUIRE( x_.at(2) == x.at(2) ); // ?= 5
  REQUIRE( x_.at(1) == x.at(1) ); // ?= 4
  REQUIRE( x_.at(0) == x.at(0) ); // ?= 3
}

TEST_CASE( "FF-LU-Decomposition 4x4 matrix", "[gf_lu_decomp_4x4]" )
{
  using namespace Numeric_lib;
  Matrix<gf_elem,2> A(4,4);
  Matrix<gf_elem,2> TMP(4,4);
  Matrix<gf_elem,2> U(4,4);
  Matrix<gf_elem,2> L(4,4);
  Matrix<gf_elem,2> P(4,4);
  for(size_t i=0; i<4; ++i) {
    L(i,i) = 1;
    P(i,i) = 1;
  }

  uint16_t degree;
  degree = lu_decomposition(P,L,U);

  REQUIRE(degree == 0);

  U(0,0) = gf_elem(1);
  U(0,1) = gf_elem(1);
  U(0,2) = gf_elem(1);
  U(0,3) = gf_elem(1);
  A=U;

  degree = lu_decomposition(P,L,U);
  REQUIRE( degree == 1 );

  TMP = L*U;
  REQUIRE( TMP == P*A );

  U(1,0) = gf_elem(1);
  U(1,1) = gf_elem(2);
  U(1,2) = gf_elem(2);
  U(1,3) = gf_elem(1);
  A(1,0) = gf_elem(1);
  A(1,1) = gf_elem(2);
  A(1,2) = gf_elem(2);
  A(1,3) = gf_elem(1);
  
  degree = lu_decomposition(P,L,U);
  REQUIRE(degree == 2);

  TMP = L*U;
  REQUIRE( TMP == P*A );

  U(2,0) = gf_elem(3);
  U(2,1) = gf_elem(2);
  U(2,2) = gf_elem(4);
  U(2,3) = gf_elem(2);
  A(2,0) = gf_elem(3);
  A(2,1) = gf_elem(2);
  A(2,2) = gf_elem(4);
  A(2,3) = gf_elem(2);

  degree = lu_decomposition(P,L,U);
  REQUIRE(degree == 3);

  TMP = L*U;
  REQUIRE( TMP == P*A );
  
  U(3,0) = gf_elem(3);
  U(3,1) = gf_elem(2);
  U(3,2) = gf_elem(4);
  U(3,3) = gf_elem(1);
  A(3,0) = gf_elem(3);
  A(3,1) = gf_elem(2);
  A(3,2) = gf_elem(4);
  A(3,3) = gf_elem(1);


  degree = lu_decomposition(P,L,U);
  REQUIRE(degree == 4);

  TMP = L*U;
  REQUIRE( TMP == P*A );
  
  // solve for x
  std::vector<gf_elem> b(4);
  b.at(0) = gf_elem(2);
  b.at(1) = gf_elem(9);
  b.at(2) = gf_elem(3);
  b.at(3) = gf_elem(9);
  print_matrix(P);
  print_matrix(L);
  print_matrix(U);
  std::vector<gf_elem> x_ = lu_solve(P,L,U,b);

  REQUIRE( x_.at(3) ==   6 );
  REQUIRE( x_.at(2) == 165 );
  REQUIRE( x_.at(1) ==  85 );
  REQUIRE( x_.at(0) == 244 );
}

TEST_CASE( "FF-LU-Decomposition 3x3 matrix with pivoting", "[gf_lu_decomp_3x3_pivoting]" )
{
  using namespace Numeric_lib;
  Matrix<gf_elem,2> A(3,3);
  Matrix<gf_elem,2> TMP(3,3);
  Matrix<gf_elem,2> U(3,3);
  Matrix<gf_elem,2> L(3,3);
  Matrix<gf_elem,2> P(3,3);
  for(size_t i=0; i<3; ++i) {
    L(i,i) = 1;
    P(i,i) = 1;
  }

  uint16_t degree;

  // degree = lu_decomposition(P,L,U);
  // printf("P0:\n"); print_matrix(P); printf("\n");
  // printf("A0:\n"); print_matrix(A); printf("\n");
  // REQUIRE(degree == 0);
  // TMP = L*U;
  // REQUIRE( TMP == P*A );

  U(0,0) = gf_elem(0);
  U(0,1) = gf_elem(0);
  U(0,2) = gf_elem(7);
  A=U;

  // degree = lu_decomposition(P,L,U);
  // printf("P1:\n"); print_matrix(P); printf("\n");
  // printf("A1:\n"); print_matrix(A); printf("\n");
  // // REQUIRE(degree == 1);
  // TMP = L*U;
  // REQUIRE( TMP == P*A );

  U(1,0) = gf_elem(13);
  A(1,0) = gf_elem(13);;
  U(1,1) = gf_elem(5);
  A(1,1) = gf_elem(5);
  U(1,2) = gf_elem(99);
  A(1,2) = gf_elem(99);

  // degree = lu_decomposition(P,L,U);
  // printf("P2:\n"); print_matrix(P); printf("\n");
  // printf("A2:\n"); print_matrix(A); printf("\n");
  // // REQUIRE(degree == 2);
  // TMP = L*U;
  // REQUIRE( TMP == P*A );

  U(2,0) = gf_elem(1);
  A(2,0) = gf_elem(1);
  U(2,1) = gf_elem(17);
  A(2,1) = gf_elem(17);
  U(2,2) = gf_elem(240);
  A(2,2) = gf_elem(240);

  degree = lu_decomposition(P,L,U);
  printf("P3:\n"); print_matrix(P); printf("\n");
  printf("A3:\n"); print_matrix(A); printf("\n");
  REQUIRE(degree == 3);
  TMP = L*U;
  REQUIRE( TMP == P*A );

  std::vector<gf_elem> x(3);
  x.at(0) = gf_elem(3);
  x.at(1) = gf_elem(4);
  x.at(2) = gf_elem(5);

  // solve for x
  std::vector<gf_elem> b(3);
  b.at(0) =  x[0]* 0 + x[1]* 0 +  x[2]*  7;
  b.at(1) =  x[0]*13 + x[1]* 5 +  x[2]* 99;
  b.at(2) =  x[0]* 1 + x[1]*17 +  x[2]*240;
  std::vector<gf_elem> x_ = lu_solve(P,L,U,b);

  REQUIRE( x_.at(2) == x.at(2) ); // ?= 5
  REQUIRE( x_.at(1) == x.at(1) ); // ?= 4
  REQUIRE( x_.at(0) == x.at(0) ); // ?= 3
}

TEST_CASE( "FF-LU-Decomposition 3x3 matrix with incremental pivoting", "[gf_lu_decomp_3x3_inc_pivoting]" )
{
  using namespace Numeric_lib;
  Matrix<gf_elem,2> A(3,3);
  Matrix<gf_elem,2> TMP(3,3);
  Matrix<gf_elem,2> U(3,3);
  Matrix<gf_elem,2> L(3,3);
  Matrix<gf_elem,2> P(3,3);
  for(size_t i=0; i<3; ++i) {
    L(i,i) = 1;
    P(i,i) = 1;
  }

  uint16_t degree;

  degree = lu_decomposition(P,L,U);
  // printf("P0:\n"); print_matrix(P); printf("\n");
  // printf("A0:\n"); print_matrix(A); printf("\n");
  REQUIRE(degree == 0);
  TMP = L*U;
  REQUIRE( TMP == P*A );

  U(0,0) = gf_elem(0);
  U(0,1) = gf_elem(0);
  U(0,2) = gf_elem(7);
  A=U;

  degree = lu_decomposition(P,L,U);
  // printf("P1:\n"); print_matrix(P); printf("\n");
  // printf("A1:\n"); print_matrix(A); printf("\n");
  REQUIRE(degree == 1);
  TMP = L*U;
  REQUIRE( TMP == P*A );

  U(1,0) = gf_elem(13);
  A(1,0) = gf_elem(13);;
  U(1,1) = gf_elem(5);
  A(1,1) = gf_elem(5);
  U(1,2) = gf_elem(99);
  A(1,2) = gf_elem(99);

  degree = lu_decomposition(P,L,U);
  // printf("P2:\n"); print_matrix(P); printf("\n");
  // printf("A2:\n"); print_matrix(A); printf("\n");
  REQUIRE(degree == 2);
  TMP = L*U;
  REQUIRE( TMP == P*A );

  U(2,0) = gf_elem(1);
  A(2,0) = gf_elem(1);
  U(2,1) = gf_elem(17);
  A(2,1) = gf_elem(17);
  U(2,2) = gf_elem(240);
  A(2,2) = gf_elem(240);

  degree = lu_decomposition(P,L,U);
  // printf("P3:\n"); print_matrix(P); printf("\n");
  // printf("A3:\n"); print_matrix(A); printf("\n");
  REQUIRE(degree == 3);
  TMP = L*U;
  REQUIRE( TMP == P*A );

  std::vector<gf_elem> x(3);
  x.at(0) = gf_elem(3);
  x.at(1) = gf_elem(4);
  x.at(2) = gf_elem(5);

  // solve for x
  std::vector<gf_elem> b(3);
  b.at(0) =  x[0]* 0 + x[1]* 0 +  x[2]*  7;
  b.at(1) =  x[0]*13 + x[1]* 5 +  x[2]* 99;
  b.at(2) =  x[0]* 1 + x[1]*17 +  x[2]*240;
  // b = P*b;
  std::vector<gf_elem> x_ = lu_solve(P,L,U,b);

  REQUIRE( x_.at(2) == x.at(2) ); // ?= 5
  REQUIRE( x_.at(1) == x.at(1) ); // ?= 4
  REQUIRE( x_.at(0) == x.at(0) ); // ?= 3
}
